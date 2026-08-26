package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"

	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rewrite/protocol"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/mat"
)

const (
	partyCount       = 3
	auxiliaryPartyID = 0
	cohortAPartyID   = 1
	cohortBPartyID   = 2
	mpcFieldBits     = 256
	divSqrtChunkSize = 64
)

type runOptions struct {
	Party          int
	Lane           int
	Input          string
	Output         string
	TimingOutput   string
	PortBase       int
	SharedKeys     string
	CKKS           string
	DataBits       int
	FractionalBits int
	Probes         int
	Seed           int64
}

type secureResult struct {
	GeneIndex      int
	PhenotypeIndex int
	BurdenP        float64
	SKATWHP        float64
}

func runSecure(options runOptions) (runErr error) {
	recorder := newTimingRecorder(options)
	runSpan := recorder.startPhase("run", "secure_process_total", "")
	defer func() {
		panicValue := recover()
		status := "success"
		if runErr != nil || panicValue != nil {
			status = "failure"
		}
		runSpan.finish(status)
		if err := recorder.write(options.TimingOutput); err != nil && runErr == nil {
			runErr = fmt.Errorf("write timing output: %w", err)
		}
		if panicValue != nil {
			panic(panicValue)
		}
	}()

	// 1. Validate the process-local runtime options.
	if options.Party < 0 || options.Party >= partyCount {
		return fmt.Errorf("party must be 0, 1, or 2")
	}
	if options.Input == "" || options.Output == "" {
		return fmt.Errorf("input and output are required")
	}
	if options.SharedKeys == "" {
		return fmt.Errorf("shared-keys is required")
	}
	if options.DataBits <= options.FractionalBits {
		return fmt.Errorf("data-bits must be greater than frac-bits")
	}
	if options.PortBase < 1 || options.PortBase > 65532 {
		return fmt.Errorf("invalid port-base %d", options.PortBase)
	}

	// 2. Load public metadata and build the natural HE batch schedule.
	inputSpan := recorder.startPhase("run", "load_public_input", "secure_process_total")
	input, err := loadPreprocessedInput(options.Input)
	if err != nil {
		inputSpan.finish("failure")
		return err
	}
	inputSpan.finish("success")
	recorder.setInput(input)

	parameterSpan := recorder.startPhase("run", "resolve_ckks_parameters", "secure_process_total")
	literal, err := securecrypto.ResolveCKKSParametersLiteral(options.CKKS)
	if err != nil {
		parameterSpan.finish("failure")
		return err
	}
	heParameters, err := ckks.NewParametersFromLiteral(literal)
	if err != nil {
		parameterSpan.finish("failure")
		return err
	}
	parameterSpan.finish("success")

	scheduleSpan := recorder.startPhase("run", "build_batch_schedule", "secure_process_total")
	fullCryptoParams, err := protocol.BuildCryptoParams(
		input.DataParams, options.Probes, heParameters.MaxSlots(),
	)
	if err != nil {
		scheduleSpan.finish("failure")
		return err
	}
	scheduleSpan.finish("success")

	networkSpan := recorder.startPhase("run", "network_init", "secure_process_total")
	networks := mpc.InitCommunication(
		"127.0.0.1",
		localhostServers(options.PortBase),
		options.Party,
		partyCount,
		1,
		options.SharedKeys,
	)
	networkSpan.finish("success")
	defer func() {
		for _, network := range networks {
			network.CloseAll()
		}
	}()

	collectiveSpan := recorder.startPhase("run", "collective_init", "secure_process_total")
	heContext := mpc.ParallelNetworks(networks).CollectiveInit(
		&heParameters,
		mpcFieldBits,
		true,
		galoisElements(heParameters, fullCryptoParams),
	)
	collectiveSpan.finish("success")

	mpcSpan := recorder.startPhase("run", "mpc_env_init", "secure_process_total")
	mpcObject := mpc.InitParallelMPCEnv(
		networks,
		mpc_core.LElem256Zero,
		options.DataBits,
		options.FractionalBits,
	)[0]
	mpcObject.SetDivSqrtMaxLen(divSqrtChunkSize)
	mpcSpan.finish("success")

	// 3. Load A/B phenotype inputs and fit the shared null model once.
	nullInputSpan := recorder.startPhase("run", "load_null_inputs", "secure_process_total")
	x, y, err := loadPartyPhenotypeInput(input, options.Party)
	if err != nil {
		nullInputSpan.finish("failure")
		return err
	}
	nullInputSpan.finish("success")

	nullSpan := recorder.startPhase("run", "null_model", "secure_process_total")
	beta, xtxInv, rss := protocol.SetupNull(mpcObject, input.DataParams, x, y)
	nullSpan.finish("success")

	// 4. Load and run one natural HE batch at a time.
	results := make([]secureResult, 0, len(input.Genes)*input.DataParams.PhenotypeCount)
	kernelSpan := recorder.startPhase("run", "secure_kernel_total", "secure_process_total")
	for batchIndex, batch := range fullCryptoParams.Batches {
		batchResults, err := runGeneBatch(
			mpcObject,
			heContext,
			input,
			batch,
			x,
			y,
			beta,
			xtxInv,
			rss,
			options,
			batchIndex,
			recorder,
		)
		if err != nil {
			kernelSpan.finish("failure")
			return err
		}
		results = append(results, batchResults...)
	}
	kernelSpan.finish("success")

	// 5. Only cohort A writes the released chromosome results.
	if options.Party == cohortAPartyID {
		writeSpan := recorder.startPhase("run", "write_secure_results", "secure_process_total")
		if err := writeSecureResults(options.Output, input, results); err != nil {
			writeSpan.finish("failure")
			return err
		}
		writeSpan.finish("success")
	}
	return nil
}

func runGeneBatch(
	mpcObject *mpc.MPC,
	heContext *securecrypto.CryptoParams,
	input preprocessedInput,
	batch protocol.GeneBatch,
	x, y *mat.Dense,
	beta, xtxInv mpc_core.RMat,
	rss mpc_core.RVec,
	options runOptions,
	batchIndex int,
	recorder *timingRecorder,
) ([]secureResult, error) {
	batchSpan := recorder.startBatchPhase(
		input, batchIndex, batch, nil,
		"batch_total", "secure_kernel_total", protocol.NoTimingPhenotype,
	)

	// 1. Convert the global batch indices to one batch-local DataParams.
	scheduleSpan := recorder.startBatchPhase(
		input, batchIndex, batch, nil,
		"build_batch_context", "batch_total", protocol.NoTimingPhenotype,
	)
	genes := make([]protocol.Gene, len(batch.GeneIndices))
	localGeneIndices := make([]int, len(batch.GeneIndices))
	for position, globalIndex := range batch.GeneIndices {
		genes[position] = input.DataParams.Genes[globalIndex]
		localGeneIndices[position] = position
	}
	batchDataParams := input.DataParams
	batchDataParams.Genes = genes
	batchCryptoParams := protocol.CryptoParams{
		Slots: heContext.GetSlots(),
		R:     options.Probes,
		Batches: []protocol.GeneBatch{{
			W:           batch.W,
			GeneIndices: localGeneIndices,
		}},
	}
	scheduleSpan.finish("success")

	// 2. Load only this party's genotype blocks for the current batch.
	loadSpan := recorder.startBatchPhase(
		input, batchIndex, batch, nil,
		"load_gene_batch", "batch_total", protocol.NoTimingPhenotype,
	)
	gp, gv, err := loadPartyGeneBatch(input, options.Party, batch.GeneIndices)
	if err != nil {
		loadSpan.finish("failure")
		batchSpan.finish("failure")
		return nil, err
	}
	privateVariantCounts := loadedPrivateVariantCounts(options.Party, gv)
	loadSpan.privateVariantCounts = privateVariantCounts
	loadSpan.finish("success")
	batchSpan.privateVariantCounts = privateVariantCounts

	// 3. Compute weights, packed statistics, final pivots, and released p-values.
	weightSpan := recorder.startBatchPhase(
		input, batchIndex, batch, privateVariantCounts,
		"compute_weights", "batch_total", protocol.NoTimingPhenotype,
	)
	weight, signedWeight := protocol.ComputeWeights(
		mpcObject, heContext, batchDataParams, gp,
	)
	weightSpan.finish("success")

	statisticsSpan := recorder.startBatchPhase(
		input, batchIndex, batch, privateVariantCounts,
		"packed_statistics", "batch_total", protocol.NoTimingPhenotype,
	)
	gpQ, gpL, gvQ, gvL, geneV, geneS1, geneS2, geneS3 :=
		protocol.ComputePackedStatistics(
			mpcObject,
			heContext,
			batchDataParams,
			batchCryptoParams,
			gp,
			gv,
			x,
			y,
			beta,
			xtxInv,
			weight,
			signedWeight,
			options.Seed,
			recorder.batchObserver(
				input, batchIndex, batch, privateVariantCounts,
			),
		)
	statisticsSpan.finish("success")

	finalizeSpan := recorder.startBatchPhase(
		input, batchIndex, batch, privateVariantCounts,
		"finalization", "batch_total", protocol.NoTimingPhenotype,
	)
	b, z := protocol.Finalize(
		mpcObject,
		batchDataParams,
		gpQ,
		gpL,
		gvQ,
		gvL,
		rss,
		geneV,
		geneS1,
		geneS2,
		geneS3,
	)
	finalizeSpan.finish("success")

	releaseSpan := recorder.startBatchPhase(
		input, batchIndex, batch, privateVariantCounts,
		"release", "batch_total", protocol.NoTimingPhenotype,
	)
	burdenP, skatP := protocol.Release(
		mpcObject, heContext, batchDataParams, b, z,
	)
	releaseSpan.finish("success")
	if options.Party != cohortAPartyID {
		batchSpan.finish("success")
		return nil, nil
	}

	// 4. Map cohort A's batch-local results back to global gene indices.
	mappingSpan := recorder.startBatchPhase(
		input, batchIndex, batch, privateVariantCounts,
		"map_released_results", "batch_total", protocol.NoTimingPhenotype,
	)
	phenotypeCount := input.DataParams.PhenotypeCount
	results := make([]secureResult, 0, len(batch.GeneIndices)*phenotypeCount)
	for position, globalIndex := range batch.GeneIndices {
		for phenotype := 0; phenotype < phenotypeCount; phenotype++ {
			localIndex := position*phenotypeCount + phenotype
			results = append(results, secureResult{
				GeneIndex:      globalIndex,
				PhenotypeIndex: phenotype,
				BurdenP:        burdenP[localIndex],
				SKATWHP:        skatP[localIndex],
			})
		}
	}
	mappingSpan.finish("success")
	batchSpan.finish("success")
	return results, nil
}

func loadedPrivateVariantCounts(party int, gv []*mat.Dense) []int {
	if party != cohortBPartyID {
		return nil
	}

	counts := make([]int, len(gv))
	for index, genotype := range gv {
		if genotype != nil {
			_, counts[index] = genotype.Dims()
		}
	}
	return counts
}

func galoisElements(
	heParameters ckks.Parameters,
	cryptoParams protocol.CryptoParams,
) []uint64 {
	unique := make(map[uint64]struct{})
	for _, batch := range cryptoParams.Batches {
		nu := cryptoParams.Slots / batch.W
		for _, element := range heParameters.GaloisElementsForInnerSum(nu, batch.W) {
			unique[element] = struct{}{}
		}
	}
	for _, element := range protocol.GtGTransformGaloisElements(
		heParameters, cryptoParams,
	) {
		unique[element] = struct{}{}
	}

	elements := make([]uint64, 0, len(unique))
	for element := range unique {
		elements = append(elements, element)
	}
	sort.Slice(elements, func(left, right int) bool {
		return elements[left] < elements[right]
	})
	return elements
}

func writeSecureResults(
	path string,
	input preprocessedInput,
	results []secureResult,
) error {
	sort.Slice(results, func(left, right int) bool {
		if results[left].GeneIndex != results[right].GeneIndex {
			return results[left].GeneIndex < results[right].GeneIndex
		}
		return results[left].PhenotypeIndex < results[right].PhenotypeIndex
	})

	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	temporary := path + ".tmp"
	file, err := os.Create(temporary)
	if err != nil {
		return err
	}
	writer := csv.NewWriter(file)
	if err := writer.Write([]string{
		"chromosome",
		"gene_index",
		"gene_id",
		"gene_symbol",
		"gene_order",
		"phenotype_index",
		"phenotype_name",
		"secure_burden_p",
		"secure_skat_wh_p",
	}); err != nil {
		file.Close()
		return err
	}
	for _, result := range results {
		gene := input.Genes[result.GeneIndex]
		if err := writer.Write([]string{
			gene.Chromosome,
			strconv.Itoa(result.GeneIndex),
			gene.ID,
			gene.Symbol,
			strconv.Itoa(gene.Order),
			strconv.Itoa(result.PhenotypeIndex),
			input.Metadata.PhenotypeColumns[result.PhenotypeIndex],
			strconv.FormatFloat(result.BurdenP, 'g', 17, 64),
			strconv.FormatFloat(result.SKATWHP, 'g', 17, 64),
		}); err != nil {
			file.Close()
			return err
		}
	}
	writer.Flush()
	if err := writer.Error(); err != nil {
		file.Close()
		return err
	}
	if err := file.Close(); err != nil {
		return err
	}
	return os.Rename(temporary, path)
}

func localhostServers(portBase int) map[string]mpc.Server {
	return map[string]mpc.Server{
		"party0": {
			IpAddr: "127.0.0.1",
			Ports: map[string]string{
				"party1": strconv.Itoa(portBase),
				"party2": strconv.Itoa(portBase + 1),
			},
		},
		"party1": {
			IpAddr: "127.0.0.1",
			Ports: map[string]string{
				"party2": strconv.Itoa(portBase + 2),
			},
		},
		"party2": {
			IpAddr: "127.0.0.1",
			Ports:  map[string]string{},
		},
	}
}
