package workflow

import (
	"fmt"
	"sort"

	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rewrite/protocol"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	partyCount       = 3
	mpcFieldBits     = 256
	divSqrtMaxLength = 64
)

type session struct {
	networks   []*mpc.Network
	heContext  *securecrypto.CryptoParams
	mpcObjects []*mpc.MPC
	metrics    *metricRecorder
	dataParams []protocol.DataParams
	beta       mpc_core.RMat
	xtxInv     mpc_core.RMat
	rss        mpc_core.RVec
}

func exchangePublicParameters(
	network *mpc.Network,
	partyID int,
	localSampleCount int,
	chromosomes []ChromosomeInput,
) (nA, nB int) {
	switch partyID {
	case auxiliaryPartyID:
		nA = network.ReceiveInt(cohortAPartyID)
		nB = network.ReceiveInt(cohortBPartyID)
		network.SendInt(nB, cohortAPartyID)
		network.SendInt(nA, cohortBPartyID)

		for index := range chromosomes {
			geneCount := network.ReceiveInt(cohortAPartyID)
			variantCounts := network.ReceiveIntVector(
				geneCount,
				cohortAPartyID,
			)
			genes := make([]protocol.Gene, geneCount)
			for geneIndex, variantCount := range variantCounts {
				genes[geneIndex].VariantCount = int(variantCount)
			}
			chromosomes[index].Genes = genes
		}
	case cohortAPartyID:
		nA = localSampleCount
		network.SendInt(nA, auxiliaryPartyID)
		nB = network.ReceiveInt(auxiliaryPartyID)

		for _, chromosome := range chromosomes {
			variantCounts := make([]uint64, len(chromosome.Genes))
			for index, gene := range chromosome.Genes {
				variantCounts[index] = uint64(gene.VariantCount)
			}
			network.SendInt(len(variantCounts), auxiliaryPartyID)
			network.SendIntVector(variantCounts, auxiliaryPartyID)
		}
	case cohortBPartyID:
		nB = localSampleCount
		network.SendInt(nB, auxiliaryPartyID)
		nA = network.ReceiveInt(auxiliaryPartyID)
	}

	return nA, nB
}

func openSession(
	config *Config,
	partyID int,
	input *PartyInput,
	sharedKeysPath string,
	metrics *metricRecorder,
) (*session, error) {
	if sharedKeysPath == "" {
		return nil, fmt.Errorf("shared keys path is required")
	}

	literal, err := securecrypto.ResolveCKKSParametersLiteral(config.CKKS)
	if err != nil {
		return nil, err
	}
	heParameters, err := ckks.NewParametersFromLiteral(literal)
	if err != nil {
		return nil, err
	}
	done := metrics.start("network_init", 0, nil)
	networks := mpc.InitCommunication(
		config.BindingIP,
		config.Servers,
		partyID,
		partyCount,
		config.MpcNumThreads,
		sharedKeysPath,
	)
	done()

	parallelNetworks := mpc.ParallelNetworks(networks)
	for _, network := range networks {
		network.EnableLogging()
	}

	done = metrics.start(
		"sample_count_exchange",
		0,
		parallelNetworks,
	)
	nA, nB := exchangePublicParameters(
		networks[0],
		partyID,
		input.SampleCount,
		input.Chromosomes,
	)
	done()

	for index, chromosome := range input.Chromosomes {
		if chromosome.CryptoParams.Slots != 0 {
			continue
		}
		dataParams := protocol.DataParams{
			Genes:          chromosome.Genes,
			C:              config.NumCov + 1,
			PhenotypeCount: len(config.PhenotypeColumns),
		}
		cryptoParams, err := protocol.BuildCryptoParams(
			dataParams,
			config.Probes,
			heParameters.MaxSlots(),
		)
		if err != nil {
			for _, network := range networks {
				network.CloseAll()
			}
			return nil, err
		}
		input.Chromosomes[index].CryptoParams = cryptoParams
	}
	galoisElements := requiredGaloisElements(
		heParameters,
		input.Chromosomes,
	)

	dataParams := make(
		[]protocol.DataParams,
		len(input.Chromosomes),
	)
	for index, chromosome := range input.Chromosomes {
		dataParams[index] = protocol.DataParams{
			Genes:          chromosome.Genes,
			NA:             nA,
			NB:             nB,
			N:              nA + nB,
			C:              config.NumCov + 1,
			PhenotypeCount: len(config.PhenotypeColumns),
		}
	}

	done = metrics.start("collective_setup", 0, parallelNetworks)
	heContext := parallelNetworks.CollectiveInit(
		&heParameters,
		mpcFieldBits,
		true,
		galoisElements,
	)
	metrics.addDuration("pubkey_gen", 0, mpc.SetupTiming.PubKey)
	metrics.addDuration("relin_key_gen", 0, mpc.SetupTiming.RelinKey)
	metrics.addDuration("rotkey_gen", 0, mpc.SetupTiming.RotKey)
	mpcObjects := mpc.InitParallelMPCEnv(
		networks,
		mpc_core.LElem256Zero,
		config.DataBits,
		config.FractionalBits,
	)

	booleanShares := mpcObjects[0].GetBooleanShareFlag()
	for _, mpcObject := range mpcObjects {
		mpcObject.SetHubPid(cohortAPartyID)
		mpcObject.SetBooleanShareFlag(booleanShares)
		mpcObject.SetDivSqrtMaxLen(divSqrtMaxLength)
	}
	done()

	done = metrics.start("null_model", 0, parallelNetworks)
	beta, xtxInv, rss := protocol.SetupNull(
		mpcObjects[0],
		dataParams[0],
		input.X,
		input.Y,
		metrics.observe(0, parallelNetworks),
	)
	done()

	return &session{
		networks:   networks,
		heContext:  heContext,
		mpcObjects: mpcObjects,
		metrics:    metrics,
		dataParams: dataParams,
		beta:       beta,
		xtxInv:     xtxInv,
		rss:        rss,
	}, nil
}

func (session *session) close() {
	for _, network := range session.networks {
		network.CloseAll()
	}
}

func requiredGaloisElements(
	heParameters ckks.Parameters,
	chromosomes []ChromosomeInput,
) []uint64 {
	unique := make(map[uint64]struct{})

	for _, chromosome := range chromosomes {
		cryptoParams := chromosome.CryptoParams

		for _, batch := range cryptoParams.Batches {
			nu := cryptoParams.Slots / batch.W
			for _, element := range heParameters.GaloisElementsForInnerSum(
				nu,
				batch.W,
			) {
				unique[element] = struct{}{}
			}
		}

		for _, element := range protocol.GtGTransformGaloisElements(
			heParameters,
			cryptoParams,
		) {
			unique[element] = struct{}{}
		}
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
