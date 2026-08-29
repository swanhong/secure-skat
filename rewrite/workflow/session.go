package workflow

import (
	"fmt"
	"sort"
	"strconv"

	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rewrite/protocol"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const (
	partyCount       = 3
	numThreads       = 1
	mpcFieldBits     = 256
	divSqrtMaxLength = 64
)

type session struct {
	networks   []*mpc.Network
	heContext  *securecrypto.CryptoParams
	mpcObject  *mpc.MPC
	metrics    *metricRecorder
	dataParams []protocol.DataParams
	beta       mpc_core.RMat
	xtxInv     mpc_core.RMat
	rss        mpc_core.RVec
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
	galoisElements := requiredGaloisElements(
		heParameters,
		input.Chromosomes,
	)

	done := metrics.start("network_init", 0, nil)
	networks := mpc.InitCommunication(
		"127.0.0.1",
		localhostServers(config.PortBase),
		partyID,
		partyCount,
		numThreads,
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
	nA, nB := exchangeSampleCounts(
		networks[0],
		partyID,
		input.SampleCount,
	)
	done()

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
			C:              len(config.CovariateColumns) + 1,
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
	mpcObject := mpc.InitParallelMPCEnv(
		networks,
		mpc_core.LElem256Zero,
		config.DataBits,
		config.FractionalBits,
	)[0]
	mpcObject.SetHubPid(cohortAPartyID)
	mpcObject.SetDivSqrtMaxLen(divSqrtMaxLength)
	done()

	done = metrics.start("null_model", 0, parallelNetworks)
	beta, xtxInv, rss := protocol.SetupNull(
		mpcObject,
		dataParams[0],
		input.X,
		input.Y,
		metrics.observe(0, parallelNetworks),
	)
	done()

	return &session{
		networks:   networks,
		heContext:  heContext,
		mpcObject:  mpcObject,
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

func exchangeSampleCounts(
	network *mpc.Network,
	partyID, localCount int,
) (nA, nB int) {
	switch partyID {
	case auxiliaryPartyID:
		nA = network.ReceiveInt(cohortAPartyID)
		nB = network.ReceiveInt(cohortBPartyID)

		network.SendInt(nB, cohortAPartyID)
		network.SendInt(nA, cohortBPartyID)

	case cohortAPartyID:
		nA = localCount
		network.SendInt(nA, auxiliaryPartyID)
		nB = network.ReceiveInt(auxiliaryPartyID)

	case cohortBPartyID:
		nB = localCount
		network.SendInt(nB, auxiliaryPartyID)
		nA = network.ReceiveInt(auxiliaryPartyID)
	}

	return nA, nB
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
