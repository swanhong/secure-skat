package gwas

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"gonum.org/v1/gonum/mat"
)

// --- CKKS level alignment (mixed-provenance ciphertexts in the low-rank path) ---

func minCipherVectorLevel(v crypto.CipherVector) int {
	if len(v) == 0 || v[0] == nil {
		return 0
	}
	min := v[0].Level()
	for _, ct := range v {
		if ct != nil && ct.Level() < min {
			min = ct.Level()
		}
	}
	return min
}

// alignCipherVectorLevels drops whichever of left/right is higher to the common
// minimum level, so CMult/CSub can combine them.
func alignCipherVectorLevels(cps *crypto.CryptoParams, left, right crypto.CipherVector) (crypto.CipherVector, crypto.CipherVector) {
	leftLevel, rightLevel := minCipherVectorLevel(left), minCipherVectorLevel(right)
	target := min(leftLevel, rightLevel)
	if leftLevel != target {
		left = crypto.DropLevel(cps, crypto.CipherMatrix{left}, target)[0]
	}
	if rightLevel != target {
		right = crypto.DropLevel(cps, crypto.CipherMatrix{right}, target)[0]
	}
	return left, right
}

// --- low-rank key-free secure score + blind locally-oriented weight ---

// scoreHE returns this party's encrypted score s = Enc(Gᵀy₀) − Σ_ℓ (GᵀX)[:,ℓ]·Enc(β̂_ℓ)
// from its plaintext contraction × the shared β̂. Each term is plaintext×cipher (key-free). pid-0 → nil.
func (ast *AssocTest) scoreHE(GtX *mat.Dense, Gty0 []float64, null skatNull) crypto.CipherVector {
	cps := ast.general.cps
	m := len(Gty0)
	if m == 0 {
		return nil
	}

	sEnc, _ := crypto.EncryptFloatVector(cps, Gty0)

	for j := 0; j < null.c; j++ {
		col := make([]float64, m)
		for k := 0; k < m; k++ {
			col[k] = GtX.At(k, j)
		}
		colPlain, _ := crypto.EncodeFloatVector(cps, col)
		term := crypto.CPMult(cps, crypto.CipherVector{null.betaRep[j]}, colPlain)

		sEnc, term = alignCipherVectorLevels(cps, sEnc, term)
		sEnc = crypto.CSub(cps, sEnc, term)
	}
	return sEnc
}

// scoreSS computes the global secure score s = Gᵀy₀ − (GᵀX)·β̂ in secret shares.
func (ast *AssocTest) scoreSS(GtX *mat.Dense, Gty0 []float64, betaSS mpc_core.RVec, m int) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	fracBits := mpcObj.GetFracBits()
	c := len(betaSS)

	gty0SS := mpc_core.InitRVec(rtype.Zero(), m)
	gtxSS := mpc_core.InitRMat(rtype.Zero(), m, c)
	if GtX != nil && len(Gty0) == m {
		for j := 0; j < m; j++ {
			gty0SS[j] = rtype.FromFloat64(Gty0[j], fracBits)
			for k := 0; k < c; k++ {
				gtxSS[j][k] = rtype.FromFloat64(GtX.At(j, k), fracBits)
			}
		}
	}

	prod := mpcObj.TruncVec(ast.ssMatVec(gtxSS, betaSS), mpcObj.GetDataBits(), fracBits) // GᵀX·β̂

	gty0SS.Sub(prod) // exact SS subtraction — cancellation is lossless here
	return gty0SS
}

// blindWeightCKKS computes w_j = 25·(1-MAF_j)^24
func (ast *AssocTest) blindWeightCKKS(dosageSum []float64, nsnps int) (crypto.CipherVector, mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()
	pid := mpcObj.GetPid()
	if nsnps == 0 {
		return nil, mpc_core.RVec{}
	}

	var mafEnc crypto.CipherVector
	if pid > 0 {
		localMAF := make([]float64, nsnps)
		inv2N := 1.0 / float64(2*ast.skatTotalNumInds())
		for j := range localMAF {
			localMAF[j] = dosageSum[j] * inv2N
		}
		mafEnc, _ = crypto.EncryptFloatVector(cps, localMAF)
		mafEnc = mpcObj.Network.AggregateCVec(cps, mafEnc)
	}

	var base24 crypto.CipherVector
	if pid > 0 {
		base := crypto.CAddConst(cps, crypto.CNeg(cps, mafEnc, false), 1.0)
		w2 := crypto.CMult(cps, base, base)
		w4 := crypto.CMult(cps, w2, w2)
		w8 := crypto.CMult(cps, w4, w4)
		w16 := crypto.CMult(cps, w8, w8)
		base24 = crypto.CMult(cps, w16, w8)
	}

	base24SS := mpcObj.CVecToSS(cps, rtype, base24, -1, len(base24), nsnps)
	weightSS := base24SS.Copy()
	weightSS.MulScalar(rtype.FromFloat64(25.0, 0))

	// Re-encryption restores the full level chain before the two multiplications in scoreCalculation.
	weightEnc := mpcObj.SSToCVec(cps, base24SS)
	if pid > 0 {
		weightEnc = crypto.CMultConst(cps, weightEnc, 25, false)
	}
	return weightEnc, weightSS
}

// scalarCiphertextToShares converts an encrypted scalar statistic to a 1-elem RVec (zero on pid 0).
func (ast *AssocTest) scalarCiphertextToShares(stat crypto.CipherVector) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	if len(stat) == 0 || stat[0] == nil {
		return mpc_core.InitRVec(rtype.Zero(), 1)
	}
	return mpcObj.CiphertextToSS(ast.general.cps, rtype, stat[0], -1, 1)
}

// scoreCalculation returns Q=Σ(ws)² and Burden=Σws as 1-elem CipherVectors.
// Reusing ws cuts the equivalent w²s² path from four ciphertext multiplications to two.
func (ast *AssocTest) scoreCalculation(score, weight crypto.CipherVector) (q, burden crypto.CipherVector) {
	cps := ast.general.cps
	weightedScore := crypto.CMult(cps, weight, score)
	return crypto.CipherVector{crypto.SqSum(cps, weightedScore)},
		crypto.CipherVector{crypto.InnerSumAll(cps, weightedScore)}
}
