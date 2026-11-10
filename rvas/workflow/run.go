package workflow

import (
	"crypto/hmac"
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rvas/protocol"
	"gonum.org/v1/gonum/mat"
)

type chromosomeResult struct {
	Chromosome int
	Genes      []protocol.Gene
	BurdenP    []float64
	SKATWHP    []float64
}

type nullModel struct {
	beta     mpc_core.RMat
	omega    mpc_core.RMat
	rssOverN mpc_core.RVec
}

func runChromosome(
	secureSession *session,
	model nullModel,
	input *PartyInput,
	chromosomeIndex int,
	seed int64,
	globalMaxMAF float64,
) (chromosomeResult, error) {
	chromosome := input.Chromosomes[chromosomeIndex]
	dataParams := secureSession.dataParams[chromosomeIndex]
	networks := mpc.ParallelNetworks(secureSession.networks)
	laneObservers := make([]func(stage string) func(), len(networks))
	for lane, network := range networks {
		laneObservers[lane] = secureSession.metrics.observe(
			chromosome.Chromosome,
			mpc.ParallelNetworks{network},
		)
	}
	observe := laneObservers[0]
	observeBatch := secureSession.metrics.observePackedWidth(
		chromosome.Chromosome,
		chromosome.CryptoParams.R,
	)
	mpcObject := secureSession.mpcObjects[0]

	var localDosage *mat.Dense
	if chromosome.GenotypeDirectory != "" {
		var err error
		localDosage, err = readPublicDosageSums(
			chromosome.GenotypeDirectory,
			chromosome.Genes,
			input.SampleCount,
		)
		if err != nil {
			return chromosomeResult{}, err
		}
	}

	filtered := globalMaxMAF < 0.5
	var cohortCounts *mat.Dense
	if filtered && localDosage != nil {
		var flat []float64
		for _, counts := range chromosome.CohortCounts {
			flat = append(flat, counts...)
		}
		cohortCounts = mat.NewDense(1, len(flat), flat)
	}

	done := secureSession.metrics.start(
		"compute_weights",
		chromosome.Chromosome,
		networks,
	)
	weight, signedWeight := protocol.ComputeWeights(
		mpcObject,
		secureSession.heContext,
		dataParams,
		localDosage,
		cohortCounts,
		globalMaxMAF,
	)
	done()

	loadGenotypes := func(
		batch protocol.GeneBatch,
	) ([]*mat.Dense, []*mat.Dense, error) {
		public, private, err := readGeneBatch(
			chromosome.GenotypeDirectory,
			chromosome.PrivateGenotypeDirectory,
			chromosome.Genes,
			batch,
			input.SampleCount,
		)
		if err != nil || !filtered || private == nil {
			return public, private, err
		}
		// B-private variants are filtered locally: a zero column equals a zero weight,
		// and the column count (hence what B reveals) is unchanged.
		for _, geneIndex := range batch.GeneIndices {
			genotypes := private[geneIndex]
			if genotypes == nil {
				continue
			}
			rows, columns := genotypes.Dims()
			counts := chromosome.PrivateCohortCounts[geneIndex]
			if len(counts) != columns {
				return nil, nil, fmt.Errorf("gene %s: %d private allele counts for %d private variants",
					chromosome.Genes[geneIndex].GeneID, len(counts), columns)
			}
			for column, count := range counts {
				if !protocol.GlobalMAFKeepPlain(count, dataParams.CohortN, globalMaxMAF) {
					genotypes.SetCol(column, make([]float64, rows))
				}
			}
		}
		return public, private, nil
	}

	done = secureSession.metrics.start(
		"packed_statistics",
		chromosome.Chromosome,
		networks,
	)
	gpQ, gpL, gvQ, gvL, geneV, geneInvS1, geneS2, geneS3, err :=
		protocol.ComputePackedStatistics(
			secureSession.mpcObjects,
			secureSession.heContext,
			dataParams,
			chromosome.CryptoParams,
			loadGenotypes,
			input.X,
			input.Y,
			model.beta,
			model.omega,
			weight,
			signedWeight,
			seed,
			laneObservers,
			observeBatch,
		)
	done()
	if err != nil {
		return chromosomeResult{}, err
	}

	done = secureSession.metrics.start(
		"finalize",
		chromosome.Chromosome,
		networks,
	)
	b, z := protocol.Finalize(
		secureSession.mpcObjects,
		dataParams,
		gpQ,
		gpL,
		gvQ,
		gvL,
		model.rssOverN,
		geneV,
		geneInvS1,
		geneS2,
		geneS3,
		observe,
	)
	done()

	done = secureSession.metrics.start(
		"release",
		chromosome.Chromosome,
		networks,
	)
	burdenP, skatWHP := protocol.Release(
		mpcObject,
		secureSession.heContext,
		dataParams,
		b,
		z,
		observe,
	)
	done()

	return chromosomeResult{
		Chromosome: chromosome.Chromosome,
		Genes:      dataParams.Genes,
		BurdenP:    burdenP,
		SKATWHP:    skatWHP,
	}, nil
}

func runParty(config *Config, partyID int) error {
	var networks []*mpc.Network
	defer func() {
		for _, network := range networks {
			network.CloseAll()
		}
	}()

	for _, ancestry := range config.Ancestries {
		if err := runAncestry(
			config,
			partyID,
			filepath.Join(config.SharedKeysPath, ancestry),
			ancestry,
			&networks,
		); err != nil {
			return fmt.Errorf("run ancestry %s: %w", ancestry, err)
		}
	}
	return nil
}

func runAncestry(
	config *Config,
	partyID int,
	sharedKeysPath string,
	ancestry string,
	sharedNetworks *[]*mpc.Network,
) error {
	printLogs := partyID == cohortAPartyID
	if config.PrintLogs != nil {
		printLogs = *config.PrintLogs
	}
	metrics := newMetricRecorder(
		fmt.Sprintf("party%d", partyID),
		ancestry,
		printLogs,
	)

	done := metrics.start("input_loading", 0, nil)
	input, err := loadPartyInput(config, partyID, ancestry)
	done()
	if err != nil {
		return fmt.Errorf(
			"load party %d input: %w",
			partyID,
			err,
		)
	}

	if *sharedNetworks == nil {
		done = metrics.start("network_init", 0, nil)
		*sharedNetworks = mpc.InitCommunication(
			config.BindingIP,
			config.Servers,
			partyID,
			partyCount,
			config.MpcNumThreads,
			sharedKeysPath,
		)
		done()
	}

	secureSession, err := openSession(
		config,
		partyID,
		input,
		sharedKeysPath,
		*sharedNetworks,
		metrics,
	)
	if err != nil {
		return fmt.Errorf(
			"open party %d session: %w",
			partyID,
			err,
		)
	}

	traceKey, err := os.ReadFile(filepath.Join(sharedKeysPath, "shared_key_global.bin"))
	if err != nil {
		return err
	}
	if len(traceKey) != 32 {
		return fmt.Errorf("trace key must contain 32 bytes")
	}
	networks := mpc.ParallelNetworks(secureSession.networks)
	done = metrics.start("null_model", 0, networks)
	beta, omega, rssOverN := protocol.SetupNull(
		secureSession.mpcObjects[0],
		secureSession.dataParams[0],
		input.X,
		input.Y,
		metrics.observe(0, networks),
	)
	done()
	model := nullModel{
		beta:     beta,
		omega:    omega,
		rssOverN: rssOverN,
	}

	// One session and null model per ancestry; each global MAF threshold reruns the
	// chromosomes with its own weights. Metrics accumulate over thresholds.
	thresholds := config.GlobalMaxMAF
	summaries := make(map[string][]chromosomeResult, len(thresholds))
	var resultErr error
	for _, threshold := range thresholds {
		label := analysisLabel(ancestry, threshold)
		results := make(
			[]chromosomeResult,
			0,
			len(input.Chromosomes),
		)
		for chromosomeIndex := range input.Chromosomes {
			chromosome := input.Chromosomes[chromosomeIndex].Chromosome
			done := metrics.start(
				"chromosome_total",
				chromosome,
				networks,
			)
			traceSeed := deriveTraceSeed(traceKey, ancestry, chromosome)
			result, err := runChromosome(
				secureSession,
				model,
				input,
				chromosomeIndex,
				traceSeed,
				threshold,
			)
			done()
			if err != nil {
				return fmt.Errorf("%s: run chromosome %d: %w", label, chromosome, err)
			}
			if partyID == cohortAPartyID || partyID == cohortBPartyID {
				results = append(results, result)
			}
			runtime.GC()
		}

		if partyID == cohortAPartyID || partyID == cohortBPartyID {
			done := metrics.start("write_results", 0, nil)
			if err := writeSecureResults(
				config.RunDir,
				label,
				results,
				config.PhenotypeColumns,
			); err != nil && resultErr == nil {
				resultErr = err
			}
			done()
		}
		summaries[label] = results
	}

	metricsPath := filepath.Join(
		config.RunDir,
		"metrics",
		ancestry,
		fmt.Sprintf("metrics_party%d.csv", partyID),
	)
	metricsErr := metrics.writeCSV(metricsPath)
	if resultErr != nil {
		return resultErr
	}
	if metricsErr != nil {
		return fmt.Errorf("write party %d metrics: %w", partyID, metricsErr)
	}
	metrics.printTimeTree()
	if printLogs && partyID != auxiliaryPartyID {
		for _, threshold := range thresholds {
			label := analysisLabel(ancestry, threshold)
			printResultSummary(summaries[label], label, partyID)
		}
	}
	return nil
}

// analysisLabel is the output subdirectory: the ancestry, plus maf<threshold> when filtered.
func analysisLabel(ancestry string, globalMaxMAF float64) string {
	return filepath.Join(ancestry, "maf"+strconv.FormatFloat(globalMaxMAF, 'g', -1, 64))
}

func printResultSummary(results []chromosomeResult, ancestry string, partyID int) {
	for index, name := range []string{"Burden", "SKAT WH"} {
		var finite []float64
		count, nonfinite, outside, zeros, ones := 0, 0, 0, 0, 0
		for _, result := range results {
			values := result.BurdenP
			if index == 1 {
				values = result.SKATWHP
			}
			for _, value := range values {
				count++
				if math.IsNaN(value) || math.IsInf(value, 0) {
					nonfinite++
					continue
				}
				finite = append(finite, value)
				if value < 0 || value > 1 {
					outside++
				}
				if value == 0 {
					zeros++
				}
				if value == 1 {
					ones++
				}
			}
		}
		rangeText, q1, median, q3 := "NA", "NA", "NA", "NA"
		if len(finite) > 0 {
			sort.Float64s(finite)
			quantile := func(p float64) string {
				position := p * float64(len(finite)-1)
				lower, upper := int(position), int(math.Ceil(position))
				weight := position - float64(lower)
				return fmt.Sprintf("%.12g", finite[lower]*(1-weight)+finite[upper]*weight)
			}
			rangeText = fmt.Sprintf("[%.12g, %.12g]", finite[0], finite[len(finite)-1])
			q1, median, q3 = quantile(0.25), quantile(0.5), quantile(0.75)
		}
		fmt.Printf("[party%d %s] %s p-values: range=%s, Q1=%s, median=%s, Q3=%s, count=%d, nonfinite=%d, outside[0,1]=%d, zeros=%d, ones=%d\n",
			partyID, ancestry, name, rangeText, q1, median, q3, count, nonfinite, outside, zeros, ones)
	}
}

func deriveTraceSeed(key []byte, ancestry string, chromosome int) int64 {
	mac := hmac.New(sha256.New, key)
	fmt.Fprintf(mac, "secure-rvas/trace/v1\x00%s\x00%d", ancestry, chromosome)
	return int64(binary.BigEndian.Uint64(mac.Sum(nil)[:8]) & ((1 << 63) - 1))
}
