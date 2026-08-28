package workflow

import (
	"fmt"

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

	weight, signedWeight := protocol.ComputeWeights(
		secureSession.mpcObject,
		secureSession.heContext,
		dataParams,
		chromosome.PublicGenotypes,
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
	)
	burdenP, skatWHP := protocol.Release(
		secureSession.mpcObject,
		secureSession.heContext,
		dataParams,
		b,
		z,
	)

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

	input, err := loadPartyInput(config, partyID)
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
	for chromosomeIndex := range input.Chromosomes {
		result := runChromosome(
			secureSession,
			input,
			chromosomeIndex,
			config.Seed,
		)
		if partyID == cohortAPartyID {
			results = append(results, result)
		}
	}

	if partyID == cohortAPartyID {
		return writeSecureResults(
			config.RunDir,
			results,
			config.PhenotypeColumns,
		)
	}
	return nil
}
