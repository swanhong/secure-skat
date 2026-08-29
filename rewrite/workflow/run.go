package workflow

import (
	"fmt"
	"path/filepath"

	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rewrite/protocol"
)

type chromosomeResult struct {
	Chromosome int
	Genes      []protocol.Gene
	BurdenP    []float64
	SKATWHP    []float64
}

func runChromosome(
	secureSession *session,
	input *PartyInput,
	chromosomeIndex int,
	seed int64,
) chromosomeResult {
	chromosome := input.Chromosomes[chromosomeIndex]
	dataParams := secureSession.dataParams[chromosomeIndex]
	networks := mpc.ParallelNetworks(secureSession.networks)
	observe := secureSession.metrics.observe(
		chromosome.Chromosome,
		networks,
	)
	observeBatch := secureSession.metrics.observePackedWidth(
		chromosome.Chromosome,
		chromosome.CryptoParams.R,
	)

	done := secureSession.metrics.start(
		"compute_weights",
		chromosome.Chromosome,
		networks,
	)
	weight, signedWeight := protocol.ComputeWeights(
		secureSession.mpcObject,
		secureSession.heContext,
		dataParams,
		chromosome.PublicGenotypes,
	)
	done()

	done = secureSession.metrics.start(
		"packed_statistics",
		chromosome.Chromosome,
		networks,
	)
	gpQ, gpL, gvQ, gvL, geneV, geneS1, geneS2, geneS3 :=
		protocol.ComputePackedStatistics(
			secureSession.mpcObject,
			secureSession.heContext,
			dataParams,
			chromosome.CryptoParams,
			chromosome.PublicGenotypes,
			chromosome.PrivateGenotypes,
			input.X,
			input.Y,
			secureSession.beta,
			secureSession.xtxInv,
			weight,
			signedWeight,
			seed,
			observe,
			observeBatch,
		)
	done()

	done = secureSession.metrics.start(
		"finalize",
		chromosome.Chromosome,
		networks,
	)
	b, z := protocol.Finalize(
		secureSession.mpcObject,
		dataParams,
		gpQ,
		gpL,
		gvQ,
		gvL,
		secureSession.rss,
		geneV,
		geneS1,
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
		secureSession.mpcObject,
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
	}
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
		result := runChromosome(
			secureSession,
			input,
			chromosomeIndex,
			config.Seed,
		)
		done()
		if partyID == cohortAPartyID {
			results = append(results, result)
		}
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
