package gwas

import (
	"math"

	mpc_core "github.com/hhcho/mpc-core"
	"gonum.org/v1/gonum/mat"
)

// gtgChunkRows bounds how many GᵀG rows enter one secret SSMultMat in the Burden-variance path, so
// a large-m gene's m×m gram is revealed in O(gtgChunkRows·m) pieces instead of one O(m²) blob (OOM).
const gtgChunkRows = 256

// burdenVarSS returns one gene's Burden-variance factor zᵀPz = z_fullᵀP z_full (P = I − X(XᵀX)⁻¹Xᵀ,
// z_full = Σⱼ wⱼ Gⱼ over the public list ∪ private variants) as an SS scalar. It never forms z (which
// is n-dim); instead it works in the small m/c space:
//
//	zᵀPz = w_pubᵀ(GᵀG)_pub w_pub  +  2·w_pubᵀ d  +  z_privᵀz_priv  −  (Xᵀz)ᵀ(XtX)⁻¹(Xᵀz)
//
// where the public GᵀG/GᵀX are federated by the "local contraction = SS share → global sum" trick
// (same as scoreSS), and the private party contributes d = G_pubᵀz_priv (m_pub), z_privᵀz_priv, and
// Xᵀz_priv (all n-contracted locally, so the private variant count m_priv stays hidden).
func (ast *AssocTest) burdenVarSS(b, nsnps int, null skatNull, X *mat.Dense, y0 []float64, priv *privateGeneLocal, privatePid int, wPubIn mpc_core.RVec, gl *geneLocal) mpc_core.RElem {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	pid := mpcObj.GetPid()
	c := null.c

	// zᵀPz/N (else wᵀGᵀGw ~2^31 at AoU overflows the fixed-point wall): quadratic pieces /N, Xᵀz pieces
	// /√N so the XtX correction also lands at /N (XtX raw). Caller normalizes Q_b by the same N → √(T/2) unchanged.
	N := float64(ast.skatTotalNumInds())
	sqrtN := math.Sqrt(N)

	// vdot = Σ aᵢbᵢ in SS (one mult+truncate, then local ring sum).
	vdot := func(a, bb mpc_core.RVec) mpc_core.RElem {
		p := mpcObj.TruncVec(mpcObj.SSMultElemVec(a, bb), db, fb)
		acc := rtype.Zero()
		for i := range p {
			acc = acc.Add(p[i])
		}
		return acc
	}

	// --- public list: GᵀX (m×c, small) as SS shares; the m×m GᵀG is kept as a local plaintext gram
	// (gg) and streamed in row-chunks below, so a full m×m *secret* matrix never materializes. ---
	var dosage []float64
	var gg *mat.Dense
	gtxT := mpc_core.InitRMat(rtype.Zero(), c, nsnps)
	if pid > 0 && nsnps > 0 {
		g := ast.localFor(b, nsnps, X, y0, gl)
		gg, dosage = g.gg, g.DosageSum
		for j := 0; j < nsnps; j++ {
			for l := 0; l < c; l++ {
				gtxT[l][j] = rtype.FromFloat64(g.GtX.At(j, l)/sqrtN, fb) // XᵀG/√N
			}
		}
	}
	wPub := wPubIn // unsigned weight reused from PART A's blockStat (same oriented dosage → same value)
	if wPub == nil {
		_, wPub = ast.blindWeightCKKS(dosage, nsnps) // collective fallback for standalone callers
	}

	// pubZZ = w_pubᵀ(GᵀG)w_pub ; pubXtz = (XᵀG_pub)w_pub (c-vector)
	pubZZ := rtype.Zero()
	pubXtz := mpc_core.InitRVec(rtype.Zero(), c)
	gtgMark := ast.metricMark()
	var wCol mpc_core.RMat
	if nsnps > 0 {
		wCol = make(mpc_core.RMat, nsnps)
		for j := 0; j < nsnps; j++ {
			wCol[j] = mpc_core.RVec{wPub[j]}
		}
		// (GᵀG)·w in row-chunks: forming the whole m×m secret gram (and its single Beaver reveal)
		// costs O(m²) memory and OOMs for m~thousands; chunking caps peak memory at O(gtgChunkRows·m)
		// while the Beaver comm total stays the same.
		gw := make(mpc_core.RVec, nsnps)
		for start := 0; start < nsnps; start += gtgChunkRows {
			end := start + gtgChunkRows
			if end > nsnps {
				end = nsnps
			}
			gtgChunk := mpc_core.InitRMat(rtype.Zero(), end-start, nsnps)
			if gg != nil {
				for j := start; j < end; j++ {
					for k := 0; k < nsnps; k++ {
						gtgChunk[j-start][k] = rtype.FromFloat64(gg.At(j, k)/N, fb) // GᵀG/N
					}
				}
			}
			res := mpcObj.SSMultMat(gtgChunk, wCol)
			for j := start; j < end; j++ {
				gw[j] = res[j-start][0]
			}
		}
		gw = mpcObj.TruncVec(gw, db, fb)
		pubZZ = vdot(wPub, gw)
	}
	ast.metricEnd("gene_burden_public_gtg", gtgMark)

	gtxMark := ast.metricMark()
	if nsnps > 0 {
		// pubXtz = (GᵀX)ᵀ·w_pub (c-vector) as one SSMultMat over the c×nsnps transpose.
		pxM := mpcObj.SSMultMat(gtxT, wCol) // c×1 (untruncated 2·fb)
		for l := 0; l < c; l++ {
			pubXtz[l] = pxM[l][0]
		}
		pubXtz = mpcObj.TruncVec(pubXtz, db, fb)
	}
	ast.metricEnd("gene_burden_public_gtx", gtxMark)

	// --- private (privatePid only): z_priv = G_priv·w_priv, contracted locally (count hidden) ---
	privateCrossMark := ast.metricMark()
	privZZ := rtype.Zero()
	crossD := mpc_core.InitRVec(rtype.Zero(), nsnps) // d = G_pubᵀ z_priv (m_pub)
	privXtz := mpc_core.InitRVec(rtype.Zero(), c)    // Xᵀ z_priv (c)
	if pid == privatePid && priv != nil {
		mp := len(priv.Weight)
		if mp > 0 {
			privZZ = rtype.FromFloat64(priv.BurdenZZ, fb)
			if len(priv.BurdenXtz) != c {
				panic("burdenVarSS: private cache missing covariate contraction")
			}
			for l, s := range priv.BurdenXtz {
				privXtz[l] = rtype.FromFloat64(s, fb)
			}
			if nsnps > 0 {
				if len(priv.BurdenCross) != nsnps {
					panic("burdenVarSS: private cache missing public-private contraction")
				}
				for j, s := range priv.BurdenCross {
					crossD[j] = rtype.FromFloat64(s, fb) // (G_pubᵀG_v/N)w_v
				}
			}
		}
	}

	crossZZ := rtype.Zero()
	if nsnps > 0 {
		crossZZ = vdot(wPub, crossD)
	}
	zz := pubZZ.Add(crossZZ).Add(crossZZ).Add(privZZ) // crossZZ added twice = the 2· term
	xtz := make(mpc_core.RVec, c)
	for l := 0; l < c; l++ {
		xtz[l] = pubXtz[l].Add(privXtz[l])
	}
	ast.metricEnd("gene_burden_private_cross", privateCrossMark)

	// zᵀPz/N = zz − (1/N)·xtzᵀΩ'xtz, Ω'=N(XtX)⁻¹ cached once by the null model.
	projectionMark := ast.metricMark()
	xtzCol := make(mpc_core.RMat, c)
	for l := 0; l < c; l++ {
		xtzCol[l] = mpc_core.RVec{xtz[l]}
	}
	aMat := mpcObj.TruncMat(mpcObj.SSMultMat(null.omp, xtzCol), db, fb)
	a := make(mpc_core.RVec, c)
	for l := 0; l < c; l++ {
		a[l] = aMat[l][0]
	}
	corr := vdot(xtz, a)
	corr = mpcObj.TruncVec(mpc_core.RVec{corr.Mul(rtype.FromFloat64(1.0/N, fb))}, db, fb)[0]
	ast.metricEnd("gene_burden_projection", projectionMark)
	return zz.Sub(corr)
}
