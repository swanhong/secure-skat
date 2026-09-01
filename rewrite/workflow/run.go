package workflow

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"

	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rewrite/protocol"
	"gonum.org/v1/gonum/mat"
)

type chromosomeResult struct {
	Chromosome int
	Genes      []protocol.Gene
	BurdenP    []float64
	SKATWHP    []float64
}

func ancestryComplete(config *Config, ancestry string) (bool, error) {
	paths := []string{
		filepath.Join(
			config.RunDir,
			"secure",
			ancestry,
			"all_secure_results.tsv",
		),
	}
	for _, chromosome := range config.Chromosomes {
		paths = append(paths, filepath.Join(
			config.RunDir,
			"secure",
			ancestry,
			fmt.Sprintf("chr%d.tsv", chromosome),
		))
	}
	for partyID := 0; partyID < partyCount; partyID++ {
		paths = append(paths, filepath.Join(
			config.RunDir,
			"metrics",
			ancestry,
			fmt.Sprintf("metrics_party%d.csv", partyID),
		))
	}

	for _, path := range paths {
		info, err := os.Stat(path)
		if os.IsNotExist(err) {
			return false, nil
		}
		if err != nil {
			return false, fmt.Errorf("inspect ancestry output %s: %w", path, err)
		}
		if !info.Mode().IsRegular() || info.Size() == 0 {
			return false, nil
		}
	}
	return true, nil
}

func runChromosome(
	secureSession *session,
	input *PartyInput,
	chromosomeIndex int,
	seed int64,
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
	)
	done()

	loadGenotypes := func(
		batch protocol.GeneBatch,
	) ([]*mat.Dense, []*mat.Dense, error) {
		return readGeneBatch(
			chromosome.GenotypeDirectory,
			chromosome.Genes,
			batch,
			input.SampleCount,
			mpcObject.GetPid() == cohortBPartyID,
		)
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
			secureSession.beta,
			secureSession.xtxInv,
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
		secureSession.rss,
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

func runParty(
	config *Config,
	partyID int,
	sharedKeysPath string,
) error {
	if err := ValidateConfig(config); err != nil {
		return fmt.Errorf("validate config: %w", err)
	}
	for _, ancestry := range config.Ancestries {
		complete, err := ancestryComplete(config, ancestry)
		if err != nil {
			return err
		}
		if complete {
			if partyID == cohortAPartyID {
				fmt.Printf("Skip completed ancestry %s\n", ancestry)
			}
			continue
		}
		if err := runAncestry(
			config,
			partyID,
			filepath.Join(sharedKeysPath, ancestry),
			ancestry,
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
) error {
	metrics := newMetricRecorder(
		fmt.Sprintf("party%d", partyID),
		ancestry,
		partyID == cohortAPartyID,
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

	secureSession, err := openSession(
		config,
		partyID,
		input,
		sharedKeysPath,
		metrics,
	)
	if err != nil {
		return fmt.Errorf(
			"open party %d session: %w",
			partyID,
			err,
		)
	}
	defer secureSession.close()

	results := make(
		[]chromosomeResult,
		0,
		len(input.Chromosomes),
	)
	networks := mpc.ParallelNetworks(secureSession.networks)
	for chromosomeIndex := range input.Chromosomes {
		chromosome := input.Chromosomes[chromosomeIndex].Chromosome
		done := metrics.start(
			"chromosome_total",
			chromosome,
			networks,
		)
		result, err := runChromosome(
			secureSession,
			input,
			chromosomeIndex,
			config.Seed,
		)
		done()
		if err != nil {
			return fmt.Errorf("run chromosome %d: %w", chromosome, err)
		}
		if partyID == cohortAPartyID {
			results = append(results, result)
		}
		runtime.GC()
	}

	var resultErr error
	if partyID == cohortAPartyID {
		done := metrics.start("write_results", 0, nil)
		resultErr = writeSecureResults(
			config.RunDir,
			ancestry,
			results,
			config.PhenotypeColumns,
		)
		done()
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
	return nil
}
