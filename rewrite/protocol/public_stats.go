package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/mat"
)

func PackBeta(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	beta mpc_core.RVec,
) securecrypto.CipherVector {
	coefficientCount := len(beta)
	slots := heParams.GetSlots()
	packed := mpc_core.InitRMat(mpcObj.GetRType().Zero(), coefficientCount, slots)
	for coefficient := range beta {
		for slot := 0; slot < slots; slot++ {
			packed[coefficient][slot] = beta[coefficient].Copy()
		}
	}

	packedMatrix := mpcObj.SSToCMat(heParams, packed)
	packedBeta := make(securecrypto.CipherVector, coefficientCount)
	if mpcObj.GetPid() == auxiliaryPartyID {
		return packedBeta
	}
	for coefficient := range packedBeta {
		packedBeta[coefficient] = packedMatrix[coefficient][0]
	}
	return packedBeta
}

func ComputeGtX(
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	gp []*mat.Dense,
	x *mat.Dense,
) ([]*mat.Dense, securecrypto.PlainVector) {
	gtx := make([]*mat.Dense, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}
		gtx[position] = new(mat.Dense)
		gtx[position].Mul(gp[geneIndex].T(), x)
	}

	_, covariateCount := x.Dims()
	gtxEncoded := make(securecrypto.PlainVector, covariateCount)
	for covariate := 0; covariate < covariateCount; covariate++ {
		columnByGene := make([][]float64, len(batch.GeneIndices))
		for position := range batch.GeneIndices {
			if gtx[position] == nil {
				continue
			}
			columnByGene[position] = mat.Col(nil, covariate, gtx[position])
		}
		packed := PackGeneBatch(dataParams, cryptoParams, batch, columnByGene)
		encoded, _ := securecrypto.EncodeFloatVector(heParams, packed)
		gtxEncoded[covariate] = encoded[0]
	}
	return gtx, gtxEncoded
}

func ComputeScore(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	gp []*mat.Dense,
	y0 mat.Vector,
	gtxEncoded securecrypto.PlainVector,
	packedBeta securecrypto.CipherVector,
) securecrypto.CipherVector {
	/*
		A: Enc(gty[A,t]) - Σa Enc(beta[a,t]) * Plain(gtx[A,a])
		B: Enc(gty[B,t]) - Σa Enc(beta[a,t]) * Plain(gtx[B,a])

		score = Aggregate(A, B)
	*/
	if mpcObj.GetPid() == auxiliaryPartyID {
		return make(securecrypto.CipherVector, 1)
	}

	gty := make([][]float64, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		gty[position] = make([]float64, variantCount)
		for variant := 0; variant < variantCount; variant++ {
			gty[position][variant] = mat.Dot(gp[geneIndex].ColView(variant), y0)
		}
	}

	packedGty := PackGeneBatch(dataParams, cryptoParams, batch, gty)
	localScore, _ := securecrypto.EncryptFloatVector(heParams, packedGty)
	for covariate := range packedBeta {
		term := securecrypto.CPMult(
			heParams,
			securecrypto.CipherVector{packedBeta[covariate]},
			securecrypto.PlainVector{gtxEncoded[covariate]},
		)
		localScore = securecrypto.CSub(heParams, localScore, term)
	}
	return mpcObj.Network.AggregateCVec(heParams, localScore)
}

func PackWeight(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	signedWeight []mpc_core.RVec,
) securecrypto.CipherVector {
	weightByGene := make([][]mpc_core.RElem, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		weightByGene[position] = signedWeight[geneIndex]
	}
	packed := packSharedGeneBatch(
		mpcObj.GetRType(),
		dataParams,
		cryptoParams,
		batch,
		weightByGene,
	)
	return mpcObj.SSToCVec(heParams, packed)
}

func ReduceQL(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	score securecrypto.CipherVector,
	packedWeight securecrypto.CipherVector,
	maskEncoded securecrypto.PlainVector,
) (geneBatchQ, geneBatchL mpc_core.RVec) {
	/*
		For each gene g and current phenotype t,
		with variants indexed by j:

		u[g,j,t] = signedWeight[g,j] * score[g,j,t]
		L[g,t]   = Σj u[g,j,t]
		Q[g,t]   = Σj u[g,j,t]^2
	*/
	// u = weighted score
	var encPackedU securecrypto.CipherVector
	if mpcObj.GetPid() != auxiliaryPartyID {
		encPackedU = securecrypto.CMult(heParams, score, packedWeight)
	}

	packedUShares := mpcObj.CVecToSS(
		heParams,
		mpcObj.GetRType(),
		encPackedU,
		auxiliaryPartyID,
		1,
		cryptoParams.Slots,
	)
	uSharesByGene := UnpackGeneBatch(
		dataParams, cryptoParams, batch, packedUShares,
	)

	// L = sum of u shares by variants
	geneBatchL = mpc_core.InitRVec(
		mpcObj.GetRType().Zero(), len(batch.GeneIndices),
	)
	for position := range uSharesByGene {
		for _, value := range uSharesByGene[position] {
			geneBatchL[position] = geneBatchL[position].Add(value)
		}
	}

	repackedUshares := packSharedGeneBatch(
		mpcObj.GetRType(),
		dataParams,
		cryptoParams,
		batch,
		uSharesByGene,
	)
	freshEncPackedU := mpcObj.SSToCVec(
		heParams, repackedUshares,
	)

	var encQRows securecrypto.CipherVector
	if mpcObj.GetPid() != auxiliaryPartyID {
		encUSquare := securecrypto.CMult(
			heParams, freshEncPackedU, freshEncPackedU,
		)
		encUSquare = securecrypto.CPMult(heParams, encUSquare, maskEncoded)

		encQRows = make(securecrypto.CipherVector, 1)
		encQRows[0] = encUSquare[0].CopyNew()
		nu := cryptoParams.Slots / batch.W
		if err := heParams.WithEvaluator(func(evaluator *ckks.Evaluator) error {
			return evaluator.InnerSum(
				encUSquare[0], nu, batch.W, encQRows[0],
			)
		}); err != nil {
			panic(err)
		}
	}

	nu := cryptoParams.Slots / batch.W
	qRowShares := mpcObj.CVecToSS(
		heParams,
		mpcObj.GetRType(),
		encQRows,
		auxiliaryPartyID,
		1,
		nu,
	)
	h := nu / len(batch.GeneIndices)
	geneBatchQ = mpc_core.InitRVec(
		mpcObj.GetRType().Zero(), len(batch.GeneIndices),
	)
	for position := range geneBatchQ {
		geneBatchQ[position] = qRowShares[position*h].Copy()
	}
	return geneBatchQ, geneBatchL
}

func ComputeGpQL(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	gp []*mat.Dense,
	y0 mat.Vector,
	gtxEncoded securecrypto.PlainVector,
	packedBeta securecrypto.CipherVector,
	packedWeight securecrypto.CipherVector,
	maskEncoded securecrypto.PlainVector,
) (geneBatchQ, geneBatchL mpc_core.RVec) {
	score := ComputeScore(
		mpcObj,
		heParams,
		dataParams,
		cryptoParams,
		batch,
		gp,
		y0,
		gtxEncoded,
		packedBeta,
	)
	return ReduceQL(
		mpcObj,
		heParams,
		dataParams,
		cryptoParams,
		batch,
		score,
		packedWeight,
		maskEncoded,
	)
}
