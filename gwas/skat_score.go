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

func dropCipherVectorToLevel(cps *crypto.CryptoParams, v crypto.CipherVector, level int) crypto.CipherVector {
	return crypto.DropLevel(cps, crypto.CipherMatrix{v}, level)[0]
}

// alignCipherVectorLevels drops whichever of left/right is higher to the common
// minimum level, so CMult/CSub can combine them.
func alignCipherVectorLevels(cps *crypto.CryptoParams, left, right crypto.CipherVector) (crypto.CipherVector, crypto.CipherVector) {
	target := minCipherVectorLevel(left)
	if r := minCipherVectorLevel(right); r < target {
		target = r
	}
	if minCipherVectorLevel(left) != target {
		left = dropCipherVectorToLevel(cps, left, target)
	}
	if minCipherVectorLevel(right) != target {
		right = dropCipherVectorToLevel(cps, right, target)
	}
	return left, right
}

// --- low-rank key-free secure score + oriented weight ---

// partyScore returns this party's encrypted score s = Enc(Gᵀy₀) − Σ_ℓ (GᵀX)[:,ℓ]·Enc(β̂_ℓ)
// from its plaintext contraction × the shared β̂. Each term is plaintext×cipher (key-free). pid-0 → nil.
func (ast *AssocTest) partyScore(GtX *mat.Dense, Gty0 []float64, null skatNull) crypto.CipherVector {
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

	betaCol := make(mpc_core.RMat, c)
	for k := 0; k < c; k++ {
		betaCol[k] = mpc_core.RVec{betaSS[k]}
	}
	gtxBeta := mpcObj.SSMultMat(gtxSS, betaCol) // (m×c)·(c×1) = global GᵀX·β̂ (frac 2×)
	prod := make(mpc_core.RVec, m)
	for j := 0; j < m; j++ {
		prod[j] = gtxBeta[j][0]
	}
	prod = mpcObj.TruncVec(prod, mpcObj.GetDataBits(), fracBits)

	s := gty0SS.Copy()
	s.Sub(prod) // exact SS subtraction — cancellation is lossless here
	return s
}

// signedWeight returns the minor-allele-oriented weight ŵ_j = t_j·w_j (t_j=−1 iff p̄_j>½),
// folding the orientation into the weight so ScoreCalculation gives both Q (sign²=1) and the
// R::SKAT-oriented Burden from one vector. Returned in SS.
func (ast *AssocTest) signedWeight(dosageSum []float64, nsnps int) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	fracBits := mpcObj.GetFracBits()

	if nsnps == 0 { // all-private gene: no public list → empty signed weight (NotLessThan would deref a[0])
		return mpc_core.RVec{}
	}

	pVec, pBarVec, w24 := ast.weightsCalculation(dosageSum, nsnps)

	// noFlip = 1 iff 1−p̄ ≥ p̄ (p̄ ≤ ½); NotLessThan(≥) matches the oracle's strict p̄>½ flip.
	noFlip := mpcObj.NotLessThan(pVec, pBarVec, mpcObj.GetBooleanShareFlag())
	noFlip.MulScalar(rtype.FromFloat64(1.0, fracBits))

	// sign = 2·noFlip − 1 ∈ {+1,−1}; public −1 subtracted on the HUB ONLY (SS convention).
	sign := noFlip.Copy()
	sign.MulScalar(rtype.FromFloat64(2.0, 0))
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		one := mpc_core.InitRVec(rtype.FromFloat64(1.0, fracBits), nsnps)
		sign.Sub(one)
	}

	signed := mpcObj.SSMultElemVec(sign, w24)
	return mpcObj.TruncVec(signed, mpcObj.GetDataBits(), fracBits)
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

// weightsCalculation returns the SKAT beta-density weights in secret shares:
// pVec=1−p̄, pBarVec=p̄, w24=25(1−MAF)^24.
func (ast *AssocTest) weightsCalculation(dosageSum []float64, nsnps_block int) (mpc_core.RVec, mpc_core.RVec, mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()

	if nsnps_block == 0 { // all-private gene: no public list; return typed empty (LessThan/Type() would deref a[0])
		return mpc_core.RVec{}, mpc_core.RVec{}, mpc_core.RVec{}
	}

	xSumRVec := mpc_core.InitRVec(rtype.Zero(), nsnps_block)

	if pid > 0 {
		for j := 0; j < nsnps_block; j++ {
			xSumRVec[j] = rtype.FromFloat64(dosageSum[j], 0) // exact integer sums
		}
	}

	// 2N is public (N public across parties), so p_j = xSum/(2N) is a secret-share ×
	// public-scalar multiply — no party interaction.
	totalIndivs := ast.skatTotalNumInds()
	inv2N := rtype.FromFloat64(1.0/float64(2*totalIndivs), mpcObj.GetFracBits())

	p_j := xSumRVec.Copy()
	p_j.MulScalar(inv2N)
	// p_j carries GetFracBits() precision (0 bits × frac bits = frac bits).

	one := rtype.FromFloat64(1.0, mpcObj.GetFracBits())
	var onesRVec mpc_core.RVec
	if pid == mpcObj.GetHubPid() {
		onesRVec = mpc_core.InitRVec(one, len(p_j))
	} else {
		onesRVec = mpc_core.InitRVec(rtype.Zero(), len(p_j))
	}
	onesRVec.Sub(p_j)
	pVec := onesRVec // p = 1 - p_j

	// SKAT weight base is 1-MAF = max(p_j, 1-p_j), invariant to allele orientation.
	// betaBase = p_j + [p_j < p]*(p - p_j).
	useBoolean := mpcObj.GetBooleanShareFlag()
	majorSelect := mpcObj.LessThan(p_j, pVec, useBoolean)
	majorSelect.MulScalar(one)

	betaBase := p_j.Copy()
	majorDelta := pVec.Copy()
	majorDelta.Sub(p_j)
	majorDelta = mpcObj.SSMultElemVec(majorDelta, majorSelect)
	majorDelta = mpcObj.TruncVec(majorDelta, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	betaBase.Add(majorDelta)

	// Compute betaBase^24 = (1 - MAF)^24 via squaring
	w2 := mpcObj.SSMultElemVec(betaBase, betaBase)
	w2 = mpcObj.TruncVec(w2, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w4 := mpcObj.SSMultElemVec(w2, w2)
	w4 = mpcObj.TruncVec(w4, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w8 := mpcObj.SSMultElemVec(w4, w4)
	w8 = mpcObj.TruncVec(w8, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w16 := mpcObj.SSMultElemVec(w8, w8)
	w16 = mpcObj.TruncVec(w16, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w24 := mpcObj.SSMultElemVec(w16, w8)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	// In standard SKAT, the weight is dbeta(MAF, 1, 25)
	// The beta density is f(x) = x^(a-1)*(1-x)^(b-1) / B(a,b)
	// B(1, 25) = 1/25. So f(x) = 25 * (1-x)^24
	betaConst := rtype.FromFloat64(25.0, mpcObj.GetFracBits())
	w24.MulScalar(betaConst)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	return pVec, p_j, w24
}

// ScoreCalculation returns Q=Σw²s² and Burden=Σw·s (1-elem CipherVectors); the
// remaining returns (S2, w2, w2S2, wS) are intermediates for debugging.
func (ast *AssocTest) ScoreCalculation(S_vec crypto.CipherVector, w_enc crypto.CipherVector) (
	crypto.CipherVector, crypto.CipherVector, crypto.CipherVector, crypto.CipherVector, crypto.CipherVector, crypto.CipherVector) {
	cryptoParams := ast.general.cps

	S2 := crypto.CMult(cryptoParams, S_vec, S_vec)
	w2 := crypto.CMult(cryptoParams, w_enc, w_enc)
	w2S2 := crypto.CMult(cryptoParams, w2, S2)
	qSkatBlock := crypto.InnerSumAll(cryptoParams, w2S2)

	wS := crypto.CMult(cryptoParams, w_enc, S_vec)
	qBurdenBlock := crypto.InnerSumAll(cryptoParams, wS)

	return crypto.CipherVector{qSkatBlock}, crypto.CipherVector{qBurdenBlock}, S2, w2, w2S2, wS
}
