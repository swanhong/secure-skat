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
	Input          string
	Output         string
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

func runSecure(options runOptions) error {
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
	input, err := loadPreprocessedInput(options.Input)
	if err != nil {
		return err
	}
	literal, err := securecrypto.ResolveCKKSParametersLiteral(options.CKKS)
	if err != nil {
		return err
	}
	heParameters, err := ckks.NewParametersFromLiteral(literal)
	if err != nil {
		return err
	}
	fullCryptoParams, err := protocol.BuildCryptoParams(
		input.DataParams, options.Probes, heParameters.MaxSlots(),
	)
	if err != nil {
		return err
	}

	networks := mpc.InitCommunication(
		"127.0.0.1",
		localhostServers(options.PortBase),
		options.Party,
		partyCount,
		1,
		options.SharedKeys,
	)
	defer func() {
		for _, network := range networks {
			network.CloseAll()
		}
	}()

	heContext := mpc.ParallelNetworks(networks).CollectiveInit(
		&heParameters,
		mpcFieldBits,
		true,
		galoisElements(heParameters, fullCryptoParams),
	)
	mpcObject := mpc.InitParallelMPCEnv(
		networks,
		mpc_core.LElem256Zero,
		options.DataBits,
		options.FractionalBits,
	)[0]
	mpcObject.SetDivSqrtMaxLen(divSqrtChunkSize)

	// 3. Load A/B phenotype inputs and fit the shared null model once.
	x, y, err := loadPartyPhenotypeInput(input, options.Party)
	if err != nil {
		return err
	}
	beta, xtxInv, rss := protocol.SetupNull(mpcObject, input.DataParams, x, y)

	// 4. Load and run one natural HE batch at a time.
	results := make([]secureResult, 0, len(input.Genes)*input.DataParams.PhenotypeCount)
	for _, batch := range fullCryptoParams.Batches {
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
		)
		if err != nil {
			return err
		}
		results = append(results, batchResults...)
	}

	// 5. Only cohort A writes the released chromosome results.
	if options.Party == cohortAPartyID {
		return writeSecureResults(options.Output, input, results)
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
) ([]secureResult, error) {
	// 1. Convert the global batch indices to one batch-local DataParams.
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

	// 2. Load only this party's genotype blocks for the current batch.
	gp, gv, err := loadPartyGeneBatch(input, options.Party, batch.GeneIndices)
	if err != nil {
		return nil, err
	}

	// 3. Compute weights, packed statistics, final pivots, and released p-values.
	weight, signedWeight := protocol.ComputeWeights(
		mpcObject, heContext, batchDataParams, gp,
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
	burdenP, skatP := protocol.Release(
		mpcObject, heContext, batchDataParams, b, z,
	)
	if options.Party != cohortAPartyID {
		return nil, nil
	}

	// 4. Map cohort A's batch-local results back to global gene indices.
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
	return results, nil
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
