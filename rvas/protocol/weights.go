package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

func mux(
	mpcObj *mpc.MPC,
	condition, whenTrue, whenFalse mpc_core.RVec,
) mpc_core.RVec {
	fixedCondition := condition.Copy()
	fixedCondition.MulScalar(mpcObj.GetRType().FromFloat64(1, mpcObj.GetFracBits()))

	difference := whenTrue.Copy()
	difference.Sub(whenFalse)

	selected := mpcObj.SSMultElemVec(fixedCondition, difference)
	selected = mpcObj.TruncVec(selected, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	result := whenFalse.Copy()
	result.Add(selected)
	return result
}

func ComputeWeights(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	localDosage *mat.Dense,
) (weight, signedWeight []mpc_core.RVec) {
	total := 0
	for _, gene := range dataParams.Genes {
		total += gene.VariantCount
	}
	weight = make([]mpc_core.RVec, len(dataParams.Genes))
	signedWeight = make([]mpc_core.RVec, len(dataParams.Genes))
	if total == 0 {
		for geneIndex := range dataParams.Genes {
			weight[geneIndex] = mpc_core.RVec{}
			signedWeight[geneIndex] = mpc_core.RVec{}
		}
		return weight, signedWeight
	}

	count := ShareSum(mpcObj, localDosage, 1, total)[0]
	rtype := mpcObj.GetRType()
	fracBits := mpcObj.GetFracBits()

	altMajor := mpcObj.NotLessThanPublic(
		count,
		rtype.FromFloat64(float64(dataParams.N), fracBits),
		mpcObj.GetBooleanShareFlag(),
	)

	// complement = 2N - count
	totalAlleles := rtype.FromFloat64(float64(2*dataParams.N), fracBits)
	complement := mpc_core.InitRVec(rtype.Zero(), total)
	complement.Sub(count)
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		complement.AddScalar(totalAlleles)
	}
	base := mux(mpcObj, altMajor, count, complement)
	base.MulScalar(rtype.FromFloat64(1/float64(2*dataParams.N), fracBits))
	base = mpcObj.TruncVec(base, mpcObj.GetDataBits(), fracBits)

	// weight <- 25 * base^24
	baseCipher := mpcObj.SSToCVec(heParams, base)
	var base8 securecrypto.CipherVector
	if mpcObj.GetPid() != auxiliaryPartyID {
		base2 := securecrypto.CMult(heParams, baseCipher, baseCipher)
		base4 := securecrypto.CMult(heParams, base2, base2)
		base8 = securecrypto.CMult(heParams, base4, base4)
	}

	// assume baselevel is not enough. Can be improved if we choose a larger CKKS modulus
	base8 = mpcObj.Network.CollectiveBootstrapVec(heParams, base8, auxiliaryPartyID)

	var base24 securecrypto.CipherVector
	if mpcObj.GetPid() != auxiliaryPartyID {
		base16 := securecrypto.CMult(heParams, base8, base8)
		base24 = securecrypto.CMult(heParams, base16, base8)
	}

	numCiphertexts := (total + heParams.GetSlots() - 1) / heParams.GetSlots()
	flatWeight := mpcObj.CVecToSS(
		heParams,
		rtype,
		base24,
		auxiliaryPartyID,
		numCiphertexts,
		total,
	)
	flatWeight.MulScalar(rtype.FromInt(25))

	negativeWeight := mpc_core.InitRVec(rtype.Zero(), len(flatWeight))
	negativeWeight.Sub(flatWeight)

	flatSignedWeight := mux(mpcObj, altMajor, negativeWeight, flatWeight)

	offset := 0
	for geneIndex, gene := range dataParams.Genes {
		end := offset + gene.VariantCount
		weight[geneIndex] = flatWeight[offset:end].Copy()
		signedWeight[geneIndex] = flatSignedWeight[offset:end].Copy()
		offset = end
	}

	return weight, signedWeight
}
