package protocol

import (
	"math"

	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
)

func Release(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	b, z mpc_core.RVec,
) (
	burdenP, skatP []float64,
) {
	/*
		Release only the two final gene-phenotype vectors.

		1. Convert b and z from shares to ciphertexts.
		2. Mask inactive slots in the final ciphertext.
		3. Collectively decrypt the masked vectors.
		4. Compute p-values locally at cohorts A and B.

		burdenP = erfc(b)
		skatP   = 0.5 * erfc(z / sqrt(2))

		Outputs use gene-major order g*q+t.
	*/

	logicalLength :=
		len(dataParams.Genes) *
			dataParams.PhenotypeCount

	// 1. Convert the two released vectors from shares to ciphertexts.
	encryptedB := mpcObj.SSToCVec(heParams, b)
	encryptedZ := mpcObj.SSToCVec(heParams, z)

	// 2. Mask inactive slots in the final ciphertext.
	tail := logicalLength % heParams.GetSlots()
	if mpcObj.GetPid() != auxiliaryPartyID && tail != 0 {
		last := len(encryptedB) - 1

		encryptedB[last] = securecrypto.MaskTrunc(
			heParams,
			encryptedB[last],
			tail,
		)
		encryptedZ[last] = securecrypto.MaskTrunc(
			heParams,
			encryptedZ[last],
			tail,
		)
	}

	// 3. Collectively decrypt only b and z.
	plainB := mpcObj.Network.CollectiveDecryptVec(
		heParams,
		encryptedB,
		0,
	)
	plainZ := mpcObj.Network.CollectiveDecryptVec(
		heParams,
		encryptedZ,
		0,
	)

	if mpcObj.GetPid() == auxiliaryPartyID {
		return nil, nil
	}

	releasedB := securecrypto.DecodeFloatVector(
		heParams,
		plainB,
	)[:logicalLength]
	releasedZ := securecrypto.DecodeFloatVector(
		heParams,
		plainZ,
	)[:logicalLength]

	// 4. Compute the two p-values locally.
	burdenP = make([]float64, logicalLength)
	skatP = make([]float64, logicalLength)

	for index := 0; index < logicalLength; index++ {
		burdenP[index] = math.Erfc(releasedB[index])
		skatP[index] =
			0.5 * math.Erfc(releasedZ[index]/math.Sqrt2)
	}

	return burdenP, skatP
}
