package workflow

import (
	"crypto/hmac"
	"crypto/sha256"
	"encoding/binary"
	"fmt"
	"os"
	"path/filepath"
	"runtime"

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
	beta   mpc_core.RMat
	xtxInv mpc_core.RMat
	rss    mpc_core.RVec
}

func runChromosome(
	secureSession *session,
	model nullModel,
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
			chromosome.PrivateGenotypeDirectory,
			chromosome.Genes,
			batch,
			input.SampleCount,
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
			model.beta,
			model.xtxInv,
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
		model.rss,
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
	beta, xtxInv, rss := protocol.SetupNull(
		secureSession.mpcObjects[0],
		secureSession.dataParams[0],
		input.X,
		input.Y,
		metrics.observe(0, networks),
	)
	done()
	model := nullModel{
		beta:   beta,
		xtxInv: xtxInv,
		rss:    rss,
	}

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
		)
		done()
		if err != nil {
			return fmt.Errorf("run chromosome %d: %w", chromosome, err)
		}
		if partyID == cohortAPartyID || partyID == cohortBPartyID {
			results = append(results, result)
		}
		runtime.GC()
	}

	var resultErr error
	if partyID == cohortAPartyID || partyID == cohortBPartyID {
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

func deriveTraceSeed(key []byte, ancestry string, chromosome int) int64 {
	mac := hmac.New(sha256.New, key)
	fmt.Fprintf(mac, "secure-rvas/trace/v1\x00%s\x00%d", ancestry, chromosome)
	return int64(binary.BigEndian.Uint64(mac.Sum(nil)[:8]) & ((1 << 63) - 1))
}
