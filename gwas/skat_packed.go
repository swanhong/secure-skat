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
	if ast.general.mpcObj[0].GetPid() == 0 {
		return nil
	}
	total := 0
	for _, size := range publicSizes {
		total += size
	}
	packed := make([]float64, total)
	offset := 0
	for gene, size := range publicSizes {
		if size == 0 {
			continue
		}
		gfs := ast.openBlockGenoStream(gene)
		if gfs == nil {
			panic("missing public genotype block")
		}
		copy(packed[offset:offset+size], orientedDosageFromStream(gfs, size))
		offset += size
	}
	return packed
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
	return mpcObj.Network.AggregateCVec(ast.general.cps, score)
}

func (ast *AssocTest) packedWindowWeight(bucket GeneBatchBucket, window GeneBatchWindow, weights []mpc_core.RVec) crypto.CipherVector {
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
	if mpcObj.GetPid() > 0 {
		normalizedMask := windowActiveMask(bucket, window, 1)
		invSqrtN := 1.0 / math.Sqrt(float64(ast.skatTotalNumInds()))
		for i := range normalizedMask {
			normalizedMask[i] *= invSqrtN
		}
		weightEnc = ast.applyPackedMask(weightEnc, normalizedMask)
	}
	return weightEnc
}

// packedWindowQL computes normalized Q=sum(x^2) and L=sum(x), x=w*s/sqrt(N).
func (ast *AssocTest) packedWindowQL(bucket GeneBatchBucket, window GeneBatchWindow, score, weightEnc crypto.CipherVector) (q, l mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	xEnc := make(crypto.CipherVector, 1)
	if mpcObj.GetPid() > 0 {
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

type privateQLBase struct {
	h, M mpc_core.RMat // phenotype-independent weighted GᵀX contractions
}

type privateQLInput struct {
	a, ell mpc_core.RVec // phenotype-specific weighted Gᵀy contractions
	v      mpc_core.RMat
}

func newPrivateQLBase(rtype mpc_core.RElem, nGenes, c int) privateQLBase {
	return privateQLBase{
		h: mpc_core.InitRMat(rtype.Zero(), nGenes, c),
		M: mpc_core.InitRMat(rtype.Zero(), nGenes*c, c),
	}
}

func newPrivateQLInput(rtype mpc_core.RElem, nGenes, c int) privateQLInput {
	return privateQLInput{
		a:   mpc_core.InitRVec(rtype.Zero(), nGenes),
		ell: mpc_core.InitRVec(rtype.Zero(), nGenes),
		v:   mpc_core.InitRMat(rtype.Zero(), nGenes, c),
	}
}

func (ast *AssocTest) addPrivateQLBase(input *privateQLBase, gene int, local *privateGeneLocal) {
	if local == nil || len(local.Weight) == 0 {
		return
	}
	mpcObj := ast.general.mpcObj[0]
	N := float64(ast.skatTotalNumInds())
	invN, invSqrtN := 1.0/N, 1.0/math.Sqrt(N)
	c := len(input.h[gene])
	h := make([]float64, c)
	M := make([]float64, c*c)
	for j, w := range local.Weight {
		for k := 0; k < c; k++ {
			wh := w * local.GtX.At(j, k)
			h[k] += wh * invSqrtN
			for k2 := 0; k2 < c; k2++ {
				M[k*c+k2] += wh * w * local.GtX.At(j, k2) * invN
			}
		}
	}
	fb := mpcObj.GetFracBits()
	rtype := mpcObj.GetRType()
	for k := 0; k < c; k++ {
		input.h[gene][k] = rtype.FromFloat64(h[k], fb)
		for k2 := 0; k2 < c; k2++ {
			input.M[gene*c+k][k2] = rtype.FromFloat64(M[k*c+k2], fb)
		}
	}
}

func (ast *AssocTest) addPrivateQL(input *privateQLInput, gene int, local *privateGeneLocal) {
	if local == nil || len(local.Weight) == 0 {
		return
	}
	mpcObj := ast.general.mpcObj[0]
	N := float64(ast.skatTotalNumInds())
	invN, invSqrtN := 1.0/N, 1.0/math.Sqrt(N)
	c := len(input.v[gene])
	v := make([]float64, c)
	a, ell := 0.0, 0.0
	for j, w := range local.Weight {
		wa := w * local.Gty0[j]
		a += wa * wa * invN
		ell += wa * invSqrtN
		for k := 0; k < c; k++ {
			v[k] += w * local.GtX.At(j, k) * wa * invN
		}
	}
	fb := mpcObj.GetFracBits()
	rtype := mpcObj.GetRType()
	input.a[gene] = rtype.FromFloat64(a, fb)
	input.ell[gene] = rtype.FromFloat64(ell, fb)
	for k := 0; k < c; k++ {
		input.v[gene][k] = rtype.FromFloat64(v[k], fb)
	}
}

// packedPrivateQL evaluates one public window of private genes from fixed-shape shares.
func (ast *AssocTest) packedPrivateQL(base privateQLBase, input privateQLInput, null skatNull) (q, l mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	nGenes, c := len(input.a), null.c

	beta := asCol(null.betaSS)
	vh := make(mpc_core.RMat, 0, 2*nGenes)
	vh = append(vh, input.v...)
	vh = append(vh, base.h...)
	vhBeta := mpcObj.TruncVec(col0(mpcObj.SSMultMat(vh, beta)), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	mBeta := mpcObj.TruncVec(col0(mpcObj.SSMultMat(base.M, beta)), mpcObj.GetDataBits(), mpcObj.GetFracBits())
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
	fedTimings.distributionName = "packed_first_pass_window_distribution"
	qCount := ast.general.config.NumPhenos
	if qCount == 0 {
		qCount = 1
	}
	resetPhenotypeTimings(qCount)

	manifest := ast.general.skatManifest
	publicSizes := ast.general.skatGeneSizes
	if len(publicSizes) != ast.general.config.GenoNumBlocks || manifest.Probes != ast.general.config.SkatPValueProbes {
		panic("packed gene manifest unavailable")
	}
	nullStarted := time.Now()
	nulls, X, y0 := ast.nullSetupMulti()
	qCount = len(nulls)
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
	q := make([]mpc_core.RVec, qCount)
	l := make([]mpc_core.RVec, qCount)
	for phenotype := 0; phenotype < qCount; phenotype++ {
		q[phenotype] = mpc_core.InitRVec(rtype.Zero(), len(publicSizes))
		l[phenotype] = mpc_core.InitRVec(rtype.Zero(), len(publicSizes))
	}
	zpz := mpc_core.InitRVec(rtype.Zero(), len(publicSizes))
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
			windowWeights := make([]mpc_core.RVec, len(window.Tiles))
			for tile, entry := range window.Tiles {
				windowWeights[tile] = weights[entry.Gene]
			}
			windowWeight := ast.packedWindowWeight(bucket, window, windowWeights)
			states := ast.preparePackedMomentGenes(window, local, windowWeights, nulls[0], probes)

			var gammaW []mpc_core.RVec
			if bucket.Mode == geneBatchHutchinson {
				hasHutchinson = true
				gammaW = ast.packedHutchinsonWave1(bucket, window, local, states, nulls[0].c)
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
			windowZpz := ast.packedBurdenVariance(states, gammaW, local, nulls[0])
			if probes > 0 {
				ast.packedMomentCorrections(states, nulls[0])
			}
			privateBase := newPrivateQLBase(rtype, len(window.Tiles), nulls[0].c)
			for tile := range window.Tiles {
				ast.addPrivateQLBase(&privateBase, tile, local[tile].Private)
			}
			for phenotype, null := range nulls {
				phenotypeStarted := time.Now()
				phenotypeLocal := phenotypeWindowLocal(local, phenotype)
				score := ast.packedWindowScore(bucket, window, phenotypeLocal, null)
				windowQ, windowL := ast.packedWindowQL(bucket, window, score, windowWeight)
				private := newPrivateQLInput(rtype, len(window.Tiles), null.c)
				for tile := range window.Tiles {
					ast.addPrivateQL(&private, tile, phenotypeLocal[tile].Private)
				}
				privateQ, privateL := ast.packedPrivateQL(privateBase, private, null)
				windowQ.Add(privateQ)
				windowL.Add(privateL)
				for tile, entry := range window.Tiles {
					q[phenotype][entry.Gene] = windowQ[tile]
					l[phenotype][entry.Gene] = windowL[tile]
				}
				fedTimings.phenotypeScoreQL[phenotype] += time.Since(phenotypeStarted)
			}
			for tile, entry := range window.Tiles {
				zpz[entry.Gene] = windowZpz[tile]
				if probes > 0 {
					moments[entry.Gene] = states[tile]
				}
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
		// Complete every Wave 1 window before any Wave 2 window.
		ast.general.mpcObj[0].AssertSync()
		for _, bucket := range manifest.Buckets {
			if bucket.Mode != geneBatchHutchinson {
				continue
			}
			for _, window := range bucket.Windows {
				wave2Mark := ast.metricMark()
				// Recompute local Gram instead of retaining every window across the barrier.
				local := ast.computeWindowLocal(window, X, nil, nil, -1, false)
				u := ast.packedWindowU(window, local, nulls[0].c)
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
	skatStat, burdenStat, skatZStat = ast.finalizePackedFederated(q, l, zpz, moments, nulls)
	ast.metricEnd("packed_finalize", finalMark)
	fedTimings.total = time.Since(started)
	return
}

func flattenGenePhenotypes(rtype mpc_core.RElem, values []mpc_core.RVec) mpc_core.RVec {
	if len(values) == 0 {
		return nil
	}
	genes := len(values[0])
	out := mpc_core.InitRVec(rtype.Zero(), genes*len(values))
	for phenotype := range values {
		if len(values[phenotype]) != genes {
			panic("phenotype statistic gene count mismatch")
		}
		for gene := 0; gene < genes; gene++ {
			out[gene*len(values)+phenotype] = values[phenotype][gene]
		}
	}
	return out
}

func repeatGenes(rtype mpc_core.RElem, values mpc_core.RVec, q int) mpc_core.RVec {
	out := mpc_core.InitRVec(rtype.Zero(), len(values)*q)
	for gene := range values {
		for phenotype := 0; phenotype < q; phenotype++ {
			out[gene*q+phenotype] = values[gene]
		}
	}
	return out
}

func (ast *AssocTest) finalizePackedFederated(qByPhenotype, lByPhenotype []mpc_core.RVec, zpz mpc_core.RVec, moments []packedMomentGene, nulls []skatNull) (skatStat, burdenStat, skatZStat crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	qCount := len(nulls)
	rss := make(mpc_core.RVec, qCount)
	for phenotype := range nulls {
		rss[phenotype] = nulls[phenotype].rssSS
	}
	scale, ok := ast.general.rareVariantScaleVectorShares(rss)
	if !ok {
		panic("packed SKAT requires more samples than covariates")
	}
	q := flattenGenePhenotypes(rtype, qByPhenotype)
	l := flattenGenePhenotypes(rtype, lByPhenotype)
	zpzFlat := repeatGenes(rtype, zpz, qCount)

	burden := ast.general.scaleRareVariantShareStats(ast.ssSquare(l), scale)
	sqrtBurden, _ := mpcObj.SqrtAndSqrtInverse(burden, false)
	_, invSqrtZpz := mpcObj.SqrtAndSqrtInverse(zpzFlat, false)
	burdenStat = ast.maskPackedOutputTail(mpcObj.SSToCVec(ast.general.cps, ast.ssMul(sqrtBurden, invSqrtZpz)), len(q))

	if ast.general.config.SkatPValueProbes == 0 {
		q.MulScalar(rtype.FromInt(ast.skatTotalNumInds()))
		q = ast.general.scaleRareVariantShareStats(q, scale)
		skatStat = ast.maskPackedOutputTail(mpcObj.SSToCVec(ast.general.cps, q), len(q))
		return
	}

	q = ast.general.scaleRareVariantShareStats(q, scale)
	s1Gene := mpc_core.InitRVec(rtype.Zero(), len(moments))
	s2Gene := mpc_core.InitRVec(rtype.Zero(), len(moments))
	s3Gene := mpc_core.InitRVec(rtype.Zero(), len(moments))
	for gene := range moments {
		s1Gene[gene] = moments[gene].tau1.Add(moments[gene].delta1)
		s2Gene[gene] = moments[gene].tau2.Add(moments[gene].delta2)
		s3Gene[gene] = moments[gene].tau3.Add(moments[gene].delta3)
	}
	s1 := repeatGenes(rtype, s1Gene, qCount)
	s2 := repeatGenes(rtype, s2Gene, qCount)
	s3 := repeatGenes(rtype, s3Gene, qCount)
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
