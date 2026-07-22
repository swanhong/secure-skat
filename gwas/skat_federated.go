package gwas

import (
	"fmt"
	"math"
	"sort"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

// --- federated SKAT with party-private variants ---
//
// Two parties hold partially-overlapping variants; one party's list is PUBLIC. Per gene a variant
// is shared (both), public_only (public-list party only — other contributes 0) or private (private
// party only, never shared). MVP+AoU is the motivating instance (MVP=public list, AoU=private).
//
//	per-gene Q = PART A (secure over the public list) + PART B (private, computed locally)
//
// Both parts are raw Σw²s² over locally minor-oriented G; the common 1/(2σ̂²) is applied once
// at the end.

// skatBetaWeight is the plaintext SKAT weight w_j = 25·(1−MAF)^24, MAF = min(p̄,1−p̄),
// p̄ = dosageSum_j/(2N); unsigned (Q only, sign²=1).
func skatBetaWeight(dosageSum []float64, totalInds int) []float64 {
	w := make([]float64, len(dosageSum))
	twoN := float64(2 * totalInds)
	for j, d := range dosageSum {
		p := d / twoN
		maf := math.Min(p, 1-p)
		w[j] = 25.0 * math.Pow(1-maf, 24)
	}
	return w
}

// privateGeneLocal memoizes the private party's per-gene plaintext work shared by Q, Burden, and
// trace moments. It never crosses the party boundary; in particular len(Weight) still reveals no
// private-variant count through an MPC message shape.
type privateGeneLocal struct {
	LocalContraction
	Weight      []float64
	BurdenZZ    float64     // w_vᵀG_vᵀG_vw_v/N
	BurdenCross []float64   // G_pᵀG_vw_v/N (m_p), empty when m_p=0
	BurdenXtz   []float64   // XᵀG_vw_v/√N (c)
	Svv         [][]float64 // G_vᵀG_v/N (m_v×m_v), built only for trace moments
	A1          [][]float64 // G_pᵀG_v/N (m_p×m_v), built only for trace moments
}

func (ast *AssocTest) computePrivateGeneLocal(G, X *mat.Dense, y0 []float64, gl *geneLocal, needMoments bool) *privateGeneLocal {
	pl := &privateGeneLocal{}
	if G == nil {
		return pl
	}
	_, m := G.Dims()
	if m == 0 {
		return pl
	}
	pl.LocalContraction = localGenotypeContract(G, X, y0)
	pl.Weight = skatBetaWeight(pl.DosageSum, ast.skatTotalNumInds())
	N := float64(ast.skatTotalNumInds())
	invN := 1.0 / N
	invSqrtN := 1.0 / math.Sqrt(N)
	_, c := X.Dims()
	pl.BurdenXtz = make([]float64, c)
	for l := 0; l < c; l++ {
		for k, w := range pl.Weight {
			pl.BurdenXtz[l] += pl.GtX.At(k, l) * w * invSqrtN
		}
	}

	// Burden-only mode needs just G_vw_v and its contractions, not the O(m_v²) moment tables.
	// This keeps the default nProbes=0 path linear in the number of private variants.
	if !needMoments {
		n, _ := G.Dims()
		z := make([]float64, n)
		for i := 0; i < n; i++ {
			for k, w := range pl.Weight {
				z[i] += G.At(i, k) * w
			}
			pl.BurdenZZ += z[i] * z[i] * invN
		}
		if gl != nil && gl.Gloc != nil {
			ng, mp := gl.Gloc.Dims()
			if ng != n {
				panic("computePrivateGeneLocal: public/private row mismatch")
			}
			pl.BurdenCross = make([]float64, mp)
			for j := 0; j < mp; j++ {
				for i := 0; i < n; i++ {
					pl.BurdenCross[j] += gl.Gloc.At(i, j) * z[i] * invN
				}
			}
		}
		return pl
	}

	toRows := func(d *mat.Dense) [][]float64 {
		r, c := d.Dims()
		out := make([][]float64, r)
		for i := 0; i < r; i++ {
			out[i] = make([]float64, c)
			for j := 0; j < c; j++ {
				out[i][j] = d.At(i, j) * invN
			}
		}
		return out
	}
	var svv mat.Dense
	svv.Mul(G.T(), G)
	pl.Svv = toRows(&svv)
	if gl != nil && gl.Gloc != nil {
		var a1 mat.Dense
		a1.Mul(gl.Gloc.T(), G)
		pl.A1 = toRows(&a1)
	}
	for k, w := range pl.Weight {
		for k2, w2 := range pl.Weight {
			pl.BurdenZZ += w * pl.Svv[k][k2] * w2
		}
	}
	if len(pl.A1) > 0 {
		pl.BurdenCross = make([]float64, len(pl.A1))
		for j := range pl.A1 {
			for k, w := range pl.Weight {
				pl.BurdenCross[j] += pl.A1[j][k] * w
			}
		}
	}
	return pl
}

// privateRawStats returns the private party's local raw SKAT = Σw²s² and Burden linear term Σw·s.
// G is already locally minor-oriented, so no signed weight or comparison is needed. The private
// variant count remains hidden: only the two scalar ciphertexts leave the owner.
func (ast *AssocTest) privateRawStats(pl *privateGeneLocal, null skatNull) (skat, burdenLin crypto.CipherVector) {
	cps := ast.general.cps

	m := 0
	if pl != nil {
		m = len(pl.Weight)
	}
	if m == 0 {
		z, _ := crypto.EncryptFloatVector(cps, []float64{0})
		zb, _ := crypto.EncryptFloatVector(cps, []float64{0})
		return z, zb
	}

	s := ast.scoreHE(pl.GtX, pl.Gty0, null)
	wEnc, _ := crypto.EncryptFloatVector(cps, pl.Weight)

	s, wEnc = alignCipherVectorLevels(cps, s, wEnc)
	skat, burdenLin = ast.scoreCalculation(s, wEnc)
	return skat, burdenLin
}

// privateBlockStat returns one gene's private-variant raw SKAT and Burden linear term as 1-elem SS
// RVecs. The private party computes them locally (privateRawStats); CiphertextToSS reshares each from
// privatePid (the other party gets only the ciphertext, pid 0 sits out). Must run for EVERY gene —
// the Enc(0) on empty genes keeps which genes have private variants indistinguishable.
func (ast *AssocTest) privateBlockStat(pl *privateGeneLocal, null skatNull, privatePid int) (skatSS, burdenSS mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()

	var skatCt, burdenCt *rlwe.Ciphertext
	if mpcObj.GetPid() == privatePid {
		skat, burdenLin := ast.privateRawStats(pl, null)
		if len(skat) > 0 {
			skatCt = skat[0]
		}
		if len(burdenLin) > 0 {
			burdenCt = burdenLin[0]
		}
	}
	skatSS = mpcObj.CiphertextToSS(ast.general.cps, rtype, skatCt, privatePid, 1)
	burdenSS = mpcObj.CiphertextToSS(ast.general.cps, rtype, burdenCt, privatePid, 1)
	return skatSS, burdenSS
}

// ComputeSKATFederatedPrivate returns per-gene Q (slot b = gene/block b). genoBlocks hold the
// public-list genotypes (the private party's aligned, public_only cols = 0); privateOnly[b] is the
// private party's private-variant genotypes for gene b (only read on privatePid, must be nil or
// length GenoNumBlocks). Caller collectively decrypts the result.
// fedTimings keeps only inclusive spans and the block distribution for the end-of-run report.
var fedTimings struct {
	initTotal, nullTotal, blocks, total time.Duration
	blockSecs                           []float64
}

// blockSecStats formats min/Q1/avg/Q3/max of per-block seconds (nearest-rank quantiles).
func blockSecStats(secs []float64) string {
	n := len(secs)
	if n == 0 {
		return "(none)"
	}
	sorted := append([]float64(nil), secs...)
	sort.Float64s(sorted)
	sum := 0.0
	for _, v := range secs {
		sum += v
	}
	q := func(p float64) float64 { return sorted[int(p*float64(n-1)+0.5)] }
	return fmt.Sprintf("min %.1fs Q1 %.1fs avg %.1fs Q3 %.1fs max %.1fs",
		sorted[0], q(0.25), sum/float64(n), q(0.75), sorted[n-1])
}

func (ast *AssocTest) ComputeSKATFederatedPrivate(privateOnly []*mat.Dense, privatePid int) (skatStat, burdenSqrtT2Stat, skatZStat crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()

	tStart := time.Now()
	null, X, y0 := ast.nullSetup()
	fedTimings.nullTotal = time.Since(tStart)
	preBlockMark := ast.metricMark()
	nullClassified := ast.fedMetrics.parentLeafDuration("null_model", "null_other")
	ast.fedMetrics.addDuration("null_other", nonNegativeDuration(fedTimings.nullTotal-nullClassified))
	log.LLvl1(fmt.Sprintf("[skat_fed] null model: %v", fedTimings.nullTotal.Round(time.Millisecond)))

	nB := ast.general.config.GenoNumBlocks
	if mpcObj.GetPid() == privatePid && privateOnly != nil && len(privateOnly) != nB {
		panic(fmt.Sprintf("ComputeSKATFederatedPrivate: privateOnly has %d blocks, want %d", len(privateOnly), nB))
	}
	skatBlockSS := mpc_core.InitRVec(rtype.Zero(), nB) // Σw²s² (SKAT raw, per gene)
	bLinBlockSS := mpc_core.InitRVec(rtype.Zero(), nB) // Σw·s  (Burden linear term, per gene)
	zpzBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)  // zᵀPz     (Burden variance, per gene; unscaled)
	nProbes := ast.general.config.SkatPValueProbes     // trace-column budget (exact basis when m<=budget); 0 = disabled
	s1B := mpc_core.InitRVec(rtype.Zero(), nB)         // SKAT kernel moments per gene (N-normalized), if enabled
	s2B := mpc_core.InitRVec(rtype.Zero(), nB)
	s3B := mpc_core.InitRVec(rtype.Zero(), nB)
	blockSecs := make([]float64, 0, nB)
	// ETA weight: the moment cost is O(m_pub²·probes), so weight each gene by nsnps² of its public block.
	var totalWork, doneWork float64
	for b := 0; b < nB; b++ {
		m := float64(ast.general.genoBlockSizes[b])
		totalWork += m * m
	}
	ast.metricEnd("pre_block_setup", preBlockMark)
	tBlocks := time.Now()
	for b := 0; b < nB; b++ {
		tb := time.Now()
		accSkat := mpc_core.InitRVec(rtype.Zero(), 1)
		accBurden := mpc_core.InitRVec(rtype.Zero(), 1)
		shapeMark := ast.metricMark()
		nsnps := ast.skatBlockNumSnps(b) // collective; public-list size for gene b
		tShape := ast.metricEnd("gene_shape_sync", shapeMark)
		// Not hub-gated: run_fed.sh tees party 2 (not the hub) to the terminal, so every party logs this.
		if nProbes > 0 { // moments: exact basis or Hutchinson + block-contracted private corrections
			log.LLvl1(fmt.Sprintf("[skat_fed] gene %d/%d start (m_pub=%d, probes=%d)", b+1, nB, nsnps, nProbes))
		} else {
			log.LLvl1(fmt.Sprintf("[skat_fed] gene %d/%d start (m_pub=%d)", b+1, nB, nsnps))
		}

		// Read + contract + Gram the public block ONCE per gene; blockStat/burdenVarSS/skatMomentsSS reuse it.
		localMark := ast.metricMark()
		gl := ast.computeGeneLocal(b, nsnps, X, y0)
		tLocal := ast.metricEnd("gene_local_public_gtg_gtx_gty", localMark)

		// PART A: secure SKAT over the public list (existing per-block path).
		var wA mpc_core.RVec // unsigned weight from PART A, reused by burden/moment paths below
		var tPublic, tPrivateLocal, tPrivateShare, tBurden, tMoments time.Duration
		if nsnps > 0 {
			t := time.Now()
			skatA, burdenA, wPub := ast.blockStat(b, nsnps, null, X, y0, gl)
			tPublic = time.Since(t)
			accSkat.Add(skatA)
			accBurden.Add(burdenA)
			wA = wPub
		}

		// PART B: private variants for this gene (uniform across all genes).
		privateLocalMark := ast.metricMark()
		var G *mat.Dense
		if mpcObj.GetPid() == privatePid && b < len(privateOnly) {
			G = orientedGenotypeLocalCopy(privateOnly[b])
		}
		pl := ast.computePrivateGeneLocal(G, X, y0, gl, nProbes > 0)
		tPrivateLocal = ast.metricEnd("gene_local_private_gtg_gtx_gty", privateLocalMark)
		privateShareMark := ast.metricMark()
		skatB, burdenB := ast.privateBlockStat(pl, null, privatePid)
		tPrivateShare = ast.metricEnd("gene_private_stat_share", privateShareMark)
		accSkat.Add(skatB)
		accBurden.Add(burdenB)

		skatBlockSS[b] = accSkat[0]
		bLinBlockSS[b] = accBurden[0]
		tBurStart := time.Now()
		zpzBlockSS[b] = ast.burdenVarSS(b, nsnps, null, X, y0, pl, privatePid, wA, gl) // Burden variance zᵀPz
		tBurden = time.Since(tBurStart)
		if nProbes > 0 { // SKAT p-value kernel moments (N-normalized)
			t := time.Now()
			s1B[b], s2B[b], s3B[b] = ast.skatMomentsSS(b, nsnps, nProbes, null, X, y0, pl, gl, wA)
			tMoments = time.Since(t)
		}
		dt := time.Since(tb)
		tPrivate := tPrivateLocal + tPrivateShare
		tOther := nonNegativeDuration(dt - tShape - tLocal - tPublic - tPrivate - tBurden - tMoments)
		blockSecs = append(blockSecs, dt.Seconds())
		// ETA weights remaining genes by m² (the O(m_pub²·R) moment cost dominates).
		doneWork += float64(nsnps) * float64(nsnps)
		elapsed := time.Since(tBlocks).Seconds()
		eta := 0.0
		pct := 100 * float64(b+1) / float64(nB)
		if doneWork > 0 {
			eta = elapsed * (totalWork - doneWork) / doneWork
		}
		if totalWork > 0 {
			pct = 100 * doneWork / totalWork
		} else { // hub has no local genotype sizes; fall back to gene-count progress
			eta = elapsed * float64(nB-b-1) / float64(b+1)
		}
		log.LLvl1(fmt.Sprintf("[skat_fed] gene %d/%d done %.3fs [shape %.3f | local %.3f | public %.3f | private-local %.3f | private-share %.3f | burden %.3f | moments %.3f | loop-overhead %.3f] | elapsed %.0fs  ETA ~%.0fm (%.0f%%)",
			b+1, nB, dt.Seconds(), tShape.Seconds(), tLocal.Seconds(), tPublic.Seconds(),
			tPrivateLocal.Seconds(), tPrivateShare.Seconds(), tBurden.Seconds(), tMoments.Seconds(),
			tOther.Seconds(), elapsed, eta/60, pct))
	}
	fedTimings.blocks = time.Since(tBlocks)
	fedTimings.blockSecs = blockSecs
	geneClassified := ast.fedMetrics.parentLeafDuration("genes", "gene_other")
	ast.fedMetrics.addDurationCount("gene_other", nonNegativeDuration(fedTimings.blocks-geneClassified), nB)
	log.LLvl1(fmt.Sprintf("[skat_fed] %d blocks: %v total | per-block %s",
		nB, fedTimings.blocks.Round(time.Millisecond), blockSecStats(blockSecs)))

	// Burden = (Σw·s)² — square the accumulated linear term (A+B) per gene before scaling.
	// Q_b/N to match zᵀPz/N: scale the score by 1/√N BEFORE squaring (b² hits the fixed-point wall at the square).
	scaleMark := ast.metricMark()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	invSqrtN := rtype.FromFloat64(1.0/math.Sqrt(float64(ast.skatTotalNumInds())), fb)
	bLinNorm := bLinBlockSS.Copy()
	bLinNorm.MulScalar(invSqrtN)
	bLinNorm = mpcObj.TruncVec(bLinNorm, db, fb)
	burdenBlockSS := mpcObj.TruncVec(mpcObj.SSSquareElemVec(bLinNorm), db, fb)

	// Common 1/(2σ̂²) applied once to both stats (linear, distributes over A+B).
	scaleSS, ok := ast.general.rareVariantScaleShares(null.rssSS)
	if !ok {
		panic("ComputeSKATFederatedPrivate: div undefined (dof = n-c <= 0); need more samples than covariates")
	}
	skatBlockSS = ast.general.scaleRareVariantShareStat(skatBlockSS, scaleSS)
	burdenBlockSS = ast.general.scaleRareVariantShareStat(burdenBlockSS, scaleSS)
	ast.metricEnd("finalize_scale", scaleMark)

	// Burden p-value: reveal ONLY √(T/2) = √(Burden/zᵀPz) = √Burden·(1/√zᵀPz). Both the Burden
	// statistic and zᵀPz stay secret-shared (never decrypted), so neither leaves and zᵀPz=2·Burden/T
	// is not derivable. √Burden and 1/√zᵀPz come from one SqrtAndSqrtInverse each; driver applies erfc.
	// SqrtAndSqrtInverse assumes strictly-positive inputs; Burden=B²/(2σ̂²)≥0 and zᵀPz≥0 hold with a
	// wide margin for any polymorphic gene (MAF-filtered upstream). A degenerate (monomorphic) gene
	// would give an unreliable — but bounded, non-crashing — p; the driver's √(T/2)>0 guard bounds output.
	burdenPMark := ast.metricMark()
	sqrtBurden, _ := mpcObj.SqrtAndSqrtInverse(burdenBlockSS, false)
	_, invSqrtZpz := mpcObj.SqrtAndSqrtInverse(zpzBlockSS, false)
	sqrtT2 := mpcObj.TruncVec(mpcObj.SSMultElemVec(sqrtBurden, invSqrtZpz), db, fb)
	ast.metricEnd("finalize_burden_pvalue", burdenPMark)

	// SKAT p-value pivot z per gene: Q_norm = (scaled SKAT stat)/N, then WH z = skatZSS(Q,S1,S2,S3).
	// Reveal only z (Q_skat, moments stay secret); driver applies p = ½erfc(z/√2).
	if nProbes > 0 {
		skatPMark := ast.metricMark()
		invN := rtype.FromFloat64(1.0/float64(ast.skatTotalNumInds()), fb)
		qNorm := skatBlockSS.Copy()
		qNorm.MulScalar(invN)
		qNorm = mpcObj.TruncVec(qNorm, db, fb)
		zB := ast.skatZSSVec(qNorm, s1B, s2B, s3B)
		skatZStat = mpcObj.SSToCVec(cps, zB)
		ast.metricEnd("finalize_skat_pvalue", skatPMark)
	}
	// SKAT-p mode reveals ONLY z (the equivalent p) — the raw SKAT statistic Q_skat must not leave.
	// Suppress it (nil ⇒ driver skips decrypt/save); the Burden path (nProbes==0) still returns Q_skat.
	packMark := ast.metricMark()
	var skatStatOut crypto.CipherVector
	if nProbes == 0 {
		skatStatOut = mpcObj.SSToCVec(cps, skatBlockSS)
	}
	burdenOut := mpcObj.SSToCVec(cps, sqrtT2)
	ast.metricEnd("finalize_output_pack", packMark)
	fedTimings.total = time.Since(tStart)
	postBlock := nonNegativeDuration(fedTimings.total - fedTimings.nullTotal -
		ast.fedMetrics.stageDuration("pre_block_setup") - fedTimings.blocks)
	finalizeClassified := ast.fedMetrics.parentLeafDuration("post_block_finalize", "finalize_other")
	ast.fedMetrics.addDuration("finalize_other", nonNegativeDuration(postBlock-finalizeClassified))
	log.LLvl1(fmt.Sprintf("[skat_fed] total compute: %v", fedTimings.total.Round(time.Millisecond)))
	return skatStatOut, burdenOut, skatZStat
}
