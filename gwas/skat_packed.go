package gwas

import (
	"math"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"gonum.org/v1/gonum/mat"
)

// localPublicDosage reads the public blocks once for the chromosome-wide flat weight circuit.
func (ast *AssocTest) localPublicDosage(publicSizes []int) []float64 {
	genes := make([][]float64, len(publicSizes))
	if ast.general.mpcObj[0].GetPid() == 0 {
		return packFlatVariants(genes)
	}
	for gene, size := range publicSizes {
		genes[gene] = make([]float64, size)
		if size == 0 {
			continue
		}
		G := orientGenotypeLocal(ast.readGenoBlockLocal(gene))
		n, _ := G.Dims()
		for j := 0; j < size; j++ {
			for i := 0; i < n; i++ {
				genes[gene][j] += G.At(i, j)
			}
		}
	}
	return packFlatVariants(genes)
}

// packedPublicWeights evaluates all public weights in flat ciphertexts and returns flat shares.
func (ast *AssocTest) packedPublicWeights(dosage []float64, n int) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	if n == 0 {
		return mpc_core.RVec{}
	}

	slots := cps.GetSlots()
	nCtx := (n + slots - 1) / slots
	var base24 crypto.CipherVector
	if mpcObj.GetPid() > 0 {
		p := make([]float64, n)
		inv2N := 1.0 / float64(2*ast.skatTotalNumInds())
		for i := range p {
			p[i] = dosage[i] * inv2N
		}
		pEnc, _ := crypto.EncryptFloatVector(cps, p)
		pEnc = mpcObj.Network.AggregateCVec(cps, pEnc)

		base := crypto.CAddConst(cps, crypto.CNeg(cps, pEnc, false), 1)
		w2 := crypto.CMult(cps, base, base)
		w4 := crypto.CMult(cps, w2, w2)
		w8 := crypto.CMult(cps, w4, w4)
		w8 = mpcObj.Network.CollectiveBootstrapVec(cps, w8, -1)
		w16 := crypto.CMult(cps, w8, w8)
		base24 = crypto.CMult(cps, w16, w8)
	}

	weight := mpcObj.CVecToSS(cps, mpcObj.GetRType(), base24, -1, nCtx, n, 1)
	weight.MulScalar(mpcObj.GetRType().FromFloat64(25, 0))
	return weight
}

func scatterFlatVariantShares(rtype mpc_core.RElem, packed mpc_core.RVec, sizes []int) []mpc_core.RVec {
	genes := make([]mpc_core.RVec, len(sizes))
	offset := 0
	for gene, size := range sizes {
		genes[gene] = mpc_core.InitRVec(rtype.Zero(), size)
		copy(genes[gene], packed[offset:offset+size])
		offset += size
	}
	if offset != len(packed) {
		panic("flat variant share size mismatch")
	}
	return genes
}

// packedWindowScore returns one aggregated score ciphertext for a public window.
func (ast *AssocTest) packedWindowScore(bucket GeneBatchBucket, window GeneBatchWindow, local []windowLocalContraction, null skatNull) crypto.CipherVector {
	mpcObj := ast.general.mpcObj[0]
	if mpcObj.GetPid() == 0 {
		return nil
	}

	slots := ast.general.cps.GetSlots()
	a := make([]float64, slots)
	H := mat.NewDense(slots, null.c, nil)
	for tileIndex, tile := range window.Tiles {
		for variant := 0; variant < tile.Variants; variant++ {
			slot := variant*bucket.L + tile.LaneBase
			a[slot] = local[tileIndex].Gty0[variant]
			for k := 0; k < null.c; k++ {
				H.Set(slot, k, local[tileIndex].GtX.At(variant, k))
			}
		}
	}

	score := ast.scoreHE(H, a, null)
	zero := crypto.CZeros(ast.general.cps, len(score))
	zero = crypto.DropLevel(ast.general.cps, crypto.CipherMatrix{zero}, score[0].Level())[0]
	score = crypto.CAdd(ast.general.cps, score, zero)
	return mpcObj.Network.AggregateCVec(ast.general.cps, score)
}

// packedWindowQL computes normalized Q=sum(x^2) and L=sum(x), x=w*s/sqrt(N).
func (ast *AssocTest) packedWindowQL(bucket GeneBatchBucket, window GeneBatchWindow, score crypto.CipherVector, weights []mpc_core.RVec) (q, l mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	tileWeights := make([]mpc_core.RMat, len(weights))
	for tile := range weights {
		tileWeights[tile] = mpc_core.InitRMat(rtype.Zero(), len(weights[tile]), 1)
		for variant := range weights[tile] {
			tileWeights[tile][variant][0] = weights[tile][variant]
		}
	}

	weightEnc := ast.windowSharesToCiphertexts(bucket, window, tileWeights, 1)
	xEnc := make(crypto.CipherVector, 1)
	if mpcObj.GetPid() > 0 {
		normalizedMask := windowActiveMask(bucket, window, 1)
		invSqrtN := 1.0 / math.Sqrt(float64(ast.skatTotalNumInds()))
		for i := range normalizedMask {
			normalizedMask[i] *= invSqrtN
		}
		weightEnc = ast.applyPackedMask(weightEnc, normalizedMask)
		score, weightEnc = alignCipherVectorLevels(ast.general.cps, score, weightEnc)
		xEnc = crypto.CMult(ast.general.cps, score, weightEnc)
	}

	xShares := ast.windowCiphertextsToShares(bucket, window, xEnc, 1)
	l = mpc_core.InitRVec(rtype.Zero(), len(window.Tiles))
	for tile := range xShares {
		for variant := range xShares[tile] {
			l[tile] = l[tile].Add(xShares[tile][variant][0])
		}
	}

	xEnc = ast.windowSharesToCiphertexts(bucket, window, xShares, 1)
	qEnc := make(crypto.CipherVector, 1)
	if mpcObj.GetPid() > 0 {
		qEnc = crypto.CMult(ast.general.cps, xEnc, xEnc)
		qEnc = ast.sumWindowRows(bucket, window, qEnc, 1)
	}
	qRows := ast.windowSumsToShares(bucket, window, qEnc, 1)
	q = mpc_core.InitRVec(rtype.Zero(), len(window.Tiles))
	for tile := range qRows {
		q[tile] = qRows[tile][0]
	}
	return q, l
}

func (ast *AssocTest) packedWindowU(window GeneBatchWindow, local []windowLocalContraction, c int) []mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	u := make([]mpc_core.RMat, len(window.Tiles))
	invN := 1.0 / float64(ast.skatTotalNumInds())
	if mpcObj.GetPid() > 0 {
		for tile, entry := range window.Tiles {
			u[tile] = mpc_core.InitRMat(rtype.Zero(), entry.Variants, c)
			for j := 0; j < entry.Variants; j++ {
				for k := 0; k < c; k++ {
					u[tile][j][k] = rtype.FromFloat64(local[tile].GtX.At(j, k)*invN, mpcObj.GetFracBits())
				}
			}
		}
		return u
	}
	for tile, entry := range window.Tiles {
		u[tile] = mpc_core.InitRMat(rtype.Zero(), entry.Variants, c)
	}
	return u
}

// packedWindowTheta computes every Theta=U Omega' in one window multiplication.
func (ast *AssocTest) packedWindowTheta(u []mpc_core.RMat, null skatNull) []mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	theta := make([]mpc_core.RMat, len(u))
	rows := 0
	for tile := range u {
		rows += len(u[tile])
	}
	if rows == 0 {
		return theta
	}
	stacked := mpc_core.InitRMat(mpcObj.GetRType().Zero(), rows, null.c)
	offset := 0
	for tile := range u {
		copy(stacked[offset:], u[tile])
		offset += len(u[tile])
	}
	flat := mpcObj.TruncMat(mpcObj.SSMultMat(stacked, null.omp), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	offset = 0
	for tile := range u {
		theta[tile] = flat[offset : offset+len(u[tile])]
		offset += len(u[tile])
	}
	return theta
}

type privateQLInput struct {
	a, ell mpc_core.RVec
	vh, M  mpc_core.RMat
}

func newPrivateQLInput(rtype mpc_core.RElem, nGenes, c int) privateQLInput {
	return privateQLInput{
		a:   mpc_core.InitRVec(rtype.Zero(), nGenes),
		ell: mpc_core.InitRVec(rtype.Zero(), nGenes),
		vh:  mpc_core.InitRMat(rtype.Zero(), 2*nGenes, c),
		M:   mpc_core.InitRMat(rtype.Zero(), nGenes*c, c),
	}
}

func (ast *AssocTest) addPrivateQL(input *privateQLInput, gene int, local *privateGeneLocal) {
	if local == nil || len(local.Weight) == 0 {
		return
	}
	mpcObj := ast.general.mpcObj[0]
	N := float64(ast.skatTotalNumInds())
	invN, invSqrtN := 1.0/N, 1.0/math.Sqrt(N)
	c := len(input.vh[gene])
	a, ell := 0.0, 0.0
	v, h := make([]float64, c), make([]float64, c)
	M := make([]float64, c*c)
	for j, w := range local.Weight {
		wa := w * local.Gty0[j]
		a += wa * wa * invN
		ell += wa * invSqrtN
		for k := 0; k < c; k++ {
			wh := w * local.GtX.At(j, k)
			v[k] += wh * wa * invN
			h[k] += wh * invSqrtN
			for k2 := 0; k2 < c; k2++ {
				M[k*c+k2] += wh * w * local.GtX.At(j, k2) * invN
			}
		}
	}
	fb := mpcObj.GetFracBits()
	rtype := mpcObj.GetRType()
	input.a[gene] = rtype.FromFloat64(a, fb)
	input.ell[gene] = rtype.FromFloat64(ell, fb)
	for k := 0; k < c; k++ {
		input.vh[gene][k] = rtype.FromFloat64(v[k], fb)
		input.vh[len(input.a)+gene][k] = rtype.FromFloat64(h[k], fb)
		for k2 := 0; k2 < c; k2++ {
			input.M[gene*c+k][k2] = rtype.FromFloat64(M[k*c+k2], fb)
		}
	}
}

// packedPrivateQL evaluates every private gene from fixed nGenes-by-c shares.
func (ast *AssocTest) packedPrivateQL(input privateQLInput, null skatNull) (q, l mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	nGenes, c := len(input.a), null.c

	beta := asCol(null.betaSS)
	vhBeta := mpcObj.TruncVec(col0(mpcObj.SSMultMat(input.vh, beta)), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	mBeta := mpcObj.TruncVec(col0(mpcObj.SSMultMat(input.M, beta)), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	betaRep := make(mpc_core.RVec, nGenes*c)
	for gene := 0; gene < nGenes; gene++ {
		copy(betaRep[gene*c:(gene+1)*c], null.betaSS)
	}
	betaMbetaTerms := ast.ssMul(betaRep, mBeta)

	q = mpc_core.InitRVec(rtype.Zero(), nGenes)
	l = mpc_core.InitRVec(rtype.Zero(), nGenes)
	for gene := 0; gene < nGenes; gene++ {
		betaMbeta := rtype.Zero()
		for k := 0; k < c; k++ {
			betaMbeta = betaMbeta.Add(betaMbetaTerms[gene*c+k])
		}
		q[gene] = input.a[gene].Sub(vhBeta[gene]).Sub(vhBeta[gene]).Add(betaMbeta)
		l[gene] = input.ell[gene].Sub(vhBeta[nGenes+gene])
	}
	return q, l
}

// computePackedFederated is the linear packed SKAT path.
func (ast *AssocTest) computePackedFederated(privateOnly []*mat.Dense, privatePid int) (skatStat, burdenStat, skatZStat crypto.CipherVector) {
	started := time.Now()
	fedTimings.nullTotal, fedTimings.blocks, fedTimings.total = 0, 0, 0
	fedTimings.blockSecs = nil

	manifest := ast.general.skatManifest
	publicSizes := ast.general.skatGeneSizes
	if len(publicSizes) != ast.general.config.GenoNumBlocks || manifest.Probes != ast.general.config.SkatPValueProbes {
		panic("packed gene manifest unavailable")
	}
	nullStarted := time.Now()
	null, X, y0 := ast.nullSetup()
	fedTimings.nullTotal = time.Since(nullStarted)
	nullClassified := ast.fedMetrics.parentLeafDuration("null_model", "null_other")
	ast.fedMetrics.addDurationCount("null_other", nonNegativeDuration(fedTimings.nullTotal-nullClassified), 1)
	if ast.general.mpcObj[0].GetPid() == privatePid && len(privateOnly) != len(publicSizes) {
		panic("private block count mismatch")
	}

	weightMark := ast.metricMark()
	nVariants := 0
	for _, size := range publicSizes {
		nVariants += size
	}
	flatWeights := ast.packedPublicWeights(ast.localPublicDosage(publicSizes), nVariants)
	weights := scatterFlatVariantShares(ast.general.mpcObj[0].GetRType(), flatWeights, publicSizes)
	ast.metricEnd("packed_weights", weightMark)

	rtype := ast.general.mpcObj[0].GetRType()
	q := mpc_core.InitRVec(rtype.Zero(), len(publicSizes))
	l := mpc_core.InitRVec(rtype.Zero(), len(publicSizes))
	zpz := mpc_core.InitRVec(rtype.Zero(), len(publicSizes))
	private := newPrivateQLInput(rtype, len(publicSizes), null.c)
	probes := ast.general.config.SkatPValueProbes
	var moments []packedMomentGene
	if probes > 0 {
		moments = make([]packedMomentGene, len(publicSizes))
	}
	hasHutchinson := false
	blocksStarted := time.Now()
	blockSecs := make([]float64, 0)
	for _, bucket := range manifest.Buckets {
		for _, window := range bucket.Windows {
			windowMark := ast.metricMark()
			local := ast.computeWindowLocal(window, X, y0, privateOnly, privatePid, probes > 0)
			score := ast.packedWindowScore(bucket, window, local, null)
			windowWeights := make([]mpc_core.RVec, len(window.Tiles))
			for tile, entry := range window.Tiles {
				windowWeights[tile] = weights[entry.Gene]
			}
			windowQ, windowL := ast.packedWindowQL(bucket, window, score, windowWeights)
			states := ast.preparePackedMomentGenes(window, local, windowWeights, null, probes)

			var gammaW []mpc_core.RVec
			if bucket.Mode == geneBatchHutchinson {
				hasHutchinson = true
				gammaW = ast.packedHutchinsonWave1(bucket, window, local, states, null.c)
			} else {
				if bucket.Mode == geneBatchRaw {
					gammaW = ast.rawGammaWeights(window, local, windowWeights)
				} else {
					gamma := ast.localGammaShares(window, local)
					weightColumns := make([]mpc_core.RMat, len(states))
					for tile := range states {
						weightColumns[tile] = asCol(states[tile].weight)
					}
					gammaWeight := ast.batchMatMul(gamma, weightColumns)
					gammaW = make([]mpc_core.RVec, len(states))
					for tile := range states {
						gammaW[tile] = col0(gammaWeight[tile])
					}
					ast.packedExactMoments(states, gamma)
				}
			}
			windowZpz := ast.packedBurdenVariance(states, gammaW, local, null)
			if probes > 0 {
				ast.packedMomentCorrections(states, null)
			}
			for tile, entry := range window.Tiles {
				q[entry.Gene] = windowQ[tile]
				l[entry.Gene] = windowL[tile]
				zpz[entry.Gene] = windowZpz[tile]
				if probes > 0 {
					moments[entry.Gene] = states[tile]
				}
				ast.addPrivateQL(&private, entry.Gene, local[tile].Private)
			}
			stage := "packed_raw_first_pass"
			if bucket.Mode == geneBatchExact {
				stage = "packed_exact_first_pass"
			} else if bucket.Mode == geneBatchHutchinson {
				stage = "packed_hutch_first_pass"
			}
			blockSecs = append(blockSecs, ast.metricEnd(stage, windowMark).Seconds())
		}
	}

	if hasHutchinson {
		ast.general.mpcObj[0].AssertSync()
		for _, bucket := range manifest.Buckets {
			if bucket.Mode != geneBatchHutchinson {
				continue
			}
			for _, window := range bucket.Windows {
				wave2Mark := ast.metricMark()
				local := ast.computeWindowLocal(window, X, y0, nil, -1, false)
				u := ast.packedWindowU(window, local, null.c)
				states := make([]packedMomentGene, len(window.Tiles))
				for tile, entry := range window.Tiles {
					states[tile] = moments[entry.Gene]
				}
				ast.packedHutchinsonWave2(bucket, window, local, states, u)
				for tile, entry := range window.Tiles {
					moments[entry.Gene] = states[tile]
				}
				ast.metricEnd("packed_hutch_wave2", wave2Mark)
			}
		}
	}
	fedTimings.blocks = time.Since(blocksStarted)
	fedTimings.blockSecs = blockSecs

	finalMark := ast.metricMark()
	privateQ, privateL := ast.packedPrivateQL(private, null)
	q.Add(privateQ)
	l.Add(privateL)
	skatStat, burdenStat, skatZStat = ast.finalizePackedFederated(q, l, zpz, moments, null)
	ast.metricEnd("packed_finalize", finalMark)
	fedTimings.total = time.Since(started)
	return
}

func (ast *AssocTest) finalizePackedFederated(q, l, zpz mpc_core.RVec, moments []packedMomentGene, null skatNull) (skatStat, burdenStat, skatZStat crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	scale, ok := ast.general.rareVariantScaleShares(null.rssSS)
	if !ok {
		panic("packed SKAT requires more samples than covariates")
	}

	burden := ast.general.scaleRareVariantShareStat(ast.ssSquare(l), scale)
	sqrtBurden, _ := mpcObj.SqrtAndSqrtInverse(burden, false)
	_, invSqrtZpz := mpcObj.SqrtAndSqrtInverse(zpz, false)
	burdenStat = ast.maskPackedOutputTail(mpcObj.SSToCVec(ast.general.cps, ast.ssMul(sqrtBurden, invSqrtZpz)), len(q))

	if ast.general.config.SkatPValueProbes == 0 {
		q.MulScalar(rtype.FromInt(ast.skatTotalNumInds()))
		q = ast.general.scaleRareVariantShareStat(q, scale)
		skatStat = ast.maskPackedOutputTail(mpcObj.SSToCVec(ast.general.cps, q), len(q))
		return
	}

	q = ast.general.scaleRareVariantShareStat(q, scale)
	s1 := mpc_core.InitRVec(rtype.Zero(), len(moments))
	s2 := mpc_core.InitRVec(rtype.Zero(), len(moments))
	s3 := mpc_core.InitRVec(rtype.Zero(), len(moments))
	for gene := range moments {
		s1[gene] = moments[gene].tau1.Add(moments[gene].delta1)
		s2[gene] = moments[gene].tau2.Add(moments[gene].delta2)
		s3[gene] = moments[gene].tau3.Add(moments[gene].delta3)
	}
	skatZStat = ast.maskPackedOutputTail(mpcObj.SSToCVec(ast.general.cps, ast.skatZSSVec(q, s1, s2, s3)), len(q))
	return
}

func (ast *AssocTest) maskPackedOutputTail(values crypto.CipherVector, n int) crypto.CipherVector {
	if ast.general.mpcObj[0].GetPid() == 0 || len(values) == 0 {
		return values
	}
	remainder := n % ast.general.cps.GetSlots()
	if remainder != 0 {
		values[len(values)-1] = crypto.MaskTrunc(ast.general.cps, values[len(values)-1], remainder)
	}
	return values
}
