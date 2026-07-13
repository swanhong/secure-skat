package gwas

import (
	"fmt"
	"math"
	"math/rand"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

// secureClamp returns min(max(x, loPub), hiPub) for SECRET x via two secure compares + branch-free
// selects (x ← cond ? bound : x). Keeps the cube-root arg in secureCbrt's window: out of it the
// inverse-Newton ring-wraps to a wrong z (a spurious hit); clamped genes are extreme-tail (p≈1 or capped).
func (ast *AssocTest) secureClampVec(x mpc_core.RVec, loPub, hiPub float64) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	hub := mpcObj.GetPid() == mpcObj.GetHubPid()
	bv := mpcObj.GetBooleanShareFlag()
	sel := func(cond mpc_core.RVec, boundPub float64, x mpc_core.RVec) mpc_core.RVec {
		cond.MulScalar(rtype.FromFloat64(1.0, fb)) // lift 0/1 bits to fixed-point {0, 1.0}
		diff := make(mpc_core.RVec, len(x))
		for i := range x {
			diff[i] = x[i].Neg()
			if hub {
				diff[i] = diff[i].Add(rtype.FromFloat64(boundPub, fb)) // bound − x (public bound on hub only)
			}
		}
		d := mpcObj.TruncVec(mpcObj.SSMultElemVec(cond, diff), db, fb)
		out := make(mpc_core.RVec, len(x))
		for i := range x {
			out[i] = x[i].Add(d[i]) // cond ? bound : x
		}
		return out
	}
	x = sel(mpcObj.LessThanPublic(x, rtype.FromFloat64(loPub, fb), bv), loPub, x)    // x<lo → lo
	x = sel(mpcObj.NotLessThanPublic(x, rtype.FromFloat64(hiPub, fb), bv), hiPub, x) // x≥hi → hi
	return x
}

func (ast *AssocTest) secureClamp(x mpc_core.RElem, loPub, hiPub float64) mpc_core.RElem {
	return ast.secureClampVec(mpc_core.RVec{x}, loPub, hiPub)[0]
}

// secureCbrt returns x^(1/3) for a SECRET x. The fixed-seed inverse-Newton converges only on
// ~[0.05, 10.1] (outside it diverges and RING-WRAPS to a finite wrong value — no NaN in the ring),
// so the caller (skatZSS) MUST clamp the arg into that window first (secureClamp to [0.1, 9]).
// Branch-free inverse-cube-root Newton — NO data-dependent range reduction (the bounded arg converges):
// y → x^(-1/3) via y ← y(4 − x·y³)/3 from seed 0.7, then x^(1/3) = x·y². 8 iters (P0: rel err ~3e-9).
func (ast *AssocTest) secureCbrtVec(x mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	hub := mpcObj.GetPid() == mpcObj.GetHubPid()
	n := len(x)

	mul := func(a, b mpc_core.RVec) mpc_core.RVec { // secret × secret → truncate to fb
		return mpcObj.TruncVec(mpcObj.SSMultElemVec(a, b), db, fb)
	}
	pmul := func(a mpc_core.RVec, cf float64) mpc_core.RVec { // secret × public constant → fb
		cfE := rtype.FromFloat64(cf, fb)
		out := make(mpc_core.RVec, len(a))
		for i := range a {
			out[i] = a[i].Mul(cfE)
		}
		return mpcObj.TruncVec(out, db, fb)
	}

	y := mpc_core.InitRVec(rtype.Zero(), n) // seed 0.7 as an additive share (hub holds it, others 0)
	if hub {
		for i := range y {
			y[i] = rtype.FromFloat64(0.7, fb)
		}
	}
	four := rtype.FromFloat64(4.0, fb)
	for it := 0; it < 8; it++ {
		y3 := mul(mul(y, y), y)
		t := mul(x, y3)
		for i := range t {
			t[i] = t[i].Neg() // −x·y³
			if hub {
				t[i] = t[i].Add(four) // 4 − x·y³
			}
		}
		y = pmul(mul(y, t), 1.0/3.0) // y·(4 − x·y³)/3
	}
	return mul(x, mul(y, y)) // x·y² = x^(1/3)
}

func (ast *AssocTest) secureCbrt(x mpc_core.RElem) mpc_core.RElem {
	return ast.secureCbrtVec(mpc_core.RVec{x})[0]
}

// Small plaintext float64 matrix ops: matMul (a·b), matMulT (a·bᵀ), matMulD (a·diag(d), column scale),
// matTrace. subMat subtracts two secret-shared matrices.
func matMul(a, b [][]float64) [][]float64 {
	p := len(a)
	q := 0
	if p > 0 {
		q = len(a[0])
	}
	r := 0
	if len(b) > 0 {
		r = len(b[0])
	}
	out := make([][]float64, p)
	for i := 0; i < p; i++ {
		out[i] = make([]float64, r)
		for k := 0; k < q; k++ {
			aik := a[i][k]
			if aik == 0 {
				continue
			}
			for j := 0; j < r; j++ {
				out[i][j] += aik * b[k][j]
			}
		}
	}
	return out
}

func matMulT(a, b [][]float64) [][]float64 { // a · bᵀ  (a: p×q, b: r×q → p×r)
	p, r := len(a), len(b)
	q := 0
	if p > 0 {
		q = len(a[0])
	}
	out := make([][]float64, p)
	for i := 0; i < p; i++ {
		out[i] = make([]float64, r)
		for j := 0; j < r; j++ {
			s := 0.0
			for k := 0; k < q; k++ {
				s += a[i][k] * b[j][k]
			}
			out[i][j] = s
		}
	}
	return out
}

func matMulD(a [][]float64, d []float64) [][]float64 { // a · diag(d): scale column k by d[k]
	out := make([][]float64, len(a))
	for i := range a {
		out[i] = make([]float64, len(a[i]))
		for k := range a[i] {
			out[i][k] = a[i][k] * d[k]
		}
	}
	return out
}

func matTrace(a [][]float64) float64 {
	s := 0.0
	for i := range a {
		if i < len(a[i]) {
			s += a[i][i]
		}
	}
	return s
}

func subMat(a, b mpc_core.RMat) mpc_core.RMat {
	out := make(mpc_core.RMat, len(a))
	for i := range a {
		out[i] = make(mpc_core.RVec, len(a[i]))
		for j := range a[i] {
			out[i][j] = a[i][j].Sub(b[i][j])
		}
	}
	return out
}

// skatMomentsSS returns the SKAT-kernel power sums S₁,S₂,S₃ = tr(Kᵏ) for one gene, K = ½D(GᵀPG/N)D,
// over the public list (nsnps, both parties) plus party B's private variants (privG). The public
// moments τₖ=tr(K_ppᵏ) come from a Hutchinson estimator over nsnps-length probes; the private and
// cross terms are added from per-gene tables (Ψ,Ξ,Π,s) party B builds in plaintext and shares. Moments
// are N-normalized (Sₖ/Nᵏ); w is the unsigned SKAT weight; S₄ is unneeded since K is PSD.
func (ast *AssocTest) skatMomentsSS(b, nsnps, nProbes int, null skatNull, X *mat.Dense, y0 []float64, privG *mat.Dense, gl *geneLocal) (S1, S2, S3 mpc_core.RElem) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	pid := mpcObj.GetPid()
	hub := pid == mpcObj.GetHubPid()
	c := null.c
	m := nsnps // public block only; the private variants enter via the contracted corrections below

	// Kernel normalized by N (public) so SS intermediates stay O(1) and the WH ratios are unchanged.
	// Hutchinson uses gg/N and GᵀX/√N; the combine uses Uₚ=GᵀX/N and Ω'=(XᵀX/N)⁻¹.
	N := float64(ast.skatTotalNumInds())
	sqrtN := math.Sqrt(N)

	// ---- SS helpers ----
	pmul := func(a mpc_core.RElem, cf float64) mpc_core.RElem { // secret × public const → truncated
		return mpcObj.TruncVec(mpc_core.RVec{a.Mul(rtype.FromFloat64(cf, fb))}, db, fb)[0]
	}
	ssMulM := func(a, bb mpc_core.RMat) mpc_core.RMat { return mpcObj.TruncMat(mpcObj.SSMultMat(a, bb), db, fb) }
	elemMulM := func(a, bb mpc_core.RMat) mpc_core.RMat { return mpcObj.TruncMat(mpcObj.SSMultElemMat(a, bb), db, fb) }
	transpose := func(a mpc_core.RMat) mpc_core.RMat {
		if len(a) == 0 {
			return a
		}
		out := mpc_core.InitRMat(rtype.Zero(), len(a[0]), len(a))
		for i := range a {
			for j := range a[i] {
				out[j][i] = a[i][j]
			}
		}
		return out
	}
	sumAllMat := func(a mpc_core.RMat) mpc_core.RElem {
		acc := rtype.Zero()
		for j := range a {
			for p := range a[j] {
				acc = acc.Add(a[j][p])
			}
		}
		return acc
	}
	trProd := func(a, bb mpc_core.RMat) mpc_core.RElem { // tr(A·B) = Σ_ij A[i][j]B[j][i] = sumAll(A⊙Bᵀ)
		if len(a) == 0 || len(bb) == 0 {
			return rtype.Zero()
		}
		return sumAllMat(elemMulM(a, transpose(bb)))
	}
	vdot := func(a, bb mpc_core.RVec) mpc_core.RElem {
		p := mpcObj.TruncVec(mpcObj.SSMultElemVec(a, bb), db, fb)
		acc := rtype.Zero()
		for i := range p {
			acc = acc.Add(p[i])
		}
		return acc
	}
	scaleRows := func(d mpc_core.RVec, M mpc_core.RMat) mpc_core.RMat { // row j × d[j]
		if len(M) == 0 {
			return M
		}
		dm := mpc_core.InitRMat(rtype.Zero(), len(M), len(M[0]))
		for j := range M {
			for k := range M[j] {
				dm[j][k] = d[j].Copy()
			}
		}
		return elemMulM(dm, M)
	}
	// shareMat wraps a party-local plaintext matrix as an additive SS share (FromFloat64); a nil vals
	// (parties without the data) contributes the zero share, so the reconstruction is the owner's value.
	shareMat := func(vals [][]float64, r, cols int) mpc_core.RMat {
		out := mpc_core.InitRMat(rtype.Zero(), r, cols)
		if vals != nil {
			for i := 0; i < r; i++ {
				for j := 0; j < cols; j++ {
					out[i][j] = rtype.FromFloat64(vals[i][j], fb)
				}
			}
		}
		return out
	}
	shareScalar := func(v float64, have bool) mpc_core.RElem {
		if have {
			return rtype.FromFloat64(v, fb)
		}
		return rtype.Zero()
	}

	// ---- public-block gram (local plaintext): gg = Gᵀpub·Gpub, gtxSS = GᵀX/√N (SS), and lc.GtX for
	// the /N combine. Gloc = this party's aligned public genotype (reused for the private A₁ term). ----
	var gg *mat.Dense
	var dosage []float64
	var Gloc *mat.Dense
	var lc LocalContraction
	haveGram := false
	gtxSS := mpc_core.InitRMat(rtype.Zero(), m, c)
	if pid > 0 && nsnps > 0 {
		g := ast.localFor(b, nsnps, X, y0, gl)
		Gloc, gg, lc, dosage = g.Gloc, g.gg, g.LocalContraction, g.DosageSum
		haveGram = true
		for j := 0; j < m; j++ {
			for l := 0; l < c; l++ {
				gtxSS[j][l] = rtype.FromFloat64(lc.GtX.At(j, l)/sqrtN, fb)
			}
		}
	}
	_, _, w := ast.weightsCalculation(dosage, nsnps) // unsigned w24 (SS), public list

	// ---- τₖ = tr(K_ppᵏ): Hutchinson over nsnps probes on the PUBLIC block (M_pp·A via chunked GᵀG). ----
	gtxT := mpc_core.InitRMat(rtype.Zero(), c, m)
	for l := 0; l < c; l++ {
		for j := 0; j < m; j++ {
			gtxT[l][j] = gtxSS[j][l]
		}
	}
	// tqdm-style progress over the gene's secure matvec work: total column-chunks across all mActionMat
	// calls = ⌈m/chunk⌉·(τ 2·probes + Ψ₁ probes + quads 2·c); log at each 10% (hub only).
	nChunks := (m + gtgChunkRows - 1) / gtgChunkRows
	momTotal := float64(nChunks * (3*nProbes + 2*c))
	momDone := 0.0
	momStart := time.Now()
	momPct := 0
	mActionMat := func(A mpc_core.RMat) mpc_core.RMat { // M_pp·A = (GᵀG/N)A − (GᵀX/√N)(XtX)⁻¹(XᵀG/√N)A
		R := len(A[0])
		ga := mpc_core.InitRMat(rtype.Zero(), m, R)
		for start := 0; start < m; start += gtgChunkRows {
			end := start + gtgChunkRows
			if end > m {
				end = m
			}
			gtgChunk := mpc_core.InitRMat(rtype.Zero(), end-start, m)
			if gg != nil {
				for j := start; j < end; j++ {
					for k := 0; k < m; k++ {
						gtgChunk[j-start][k] = rtype.FromFloat64(gg.At(j, k)/N, fb)
					}
				}
			}
			res := mpcObj.SSMultMat(gtgChunk, A)
			for j := start; j < end; j++ {
				for p := 0; p < R; p++ {
					ga[j][p] = res[j-start][p]
				}
			}
			momDone += float64(R)
			if hub && momTotal > 0 {
				if p := int(10 * momDone / momTotal); p > momPct {
					momPct = p
					log.LLvl1(fmt.Sprintf("[skat_fed]   moments %3d%% (elapsed %.0fs)", p*10, time.Since(momStart).Seconds()))
				}
			}
		}
		ga = mpcObj.TruncMat(ga, db, fb)
		xtGA := mpcObj.TruncMat(mpcObj.SSMultMat(gtxT, A), db, fb)
		solMat := ast.choleskySolveMat(null.xtxL, null.xtxDinv, xtGA)
		gxsol := mpcObj.TruncMat(mpcObj.SSMultMat(gtxSS, solMat), db, fb)
		Mv := mpc_core.InitRMat(rtype.Zero(), m, R)
		for j := 0; j < m; j++ {
			for p := 0; p < R; p++ {
				Mv[j][p] = ga[j][p].Sub(gxsol[j][p])
			}
		}
		return Mv
	}
	vmulMat := func(a, bb mpc_core.RMat) mpc_core.RMat { return mpcObj.TruncMat(mpcObj.SSMultElemMat(a, bb), db, fb) }
	vpmulMat := func(a mpc_core.RMat, cf float64) mpc_core.RMat {
		cfE := rtype.FromFloat64(cf, fb)
		out := mpc_core.InitRMat(rtype.Zero(), len(a), len(a[0]))
		for j := range a {
			for p := range a[j] {
				out[j][p] = a[j][p].Mul(cfE)
			}
		}
		return mpcObj.TruncMat(out, db, fb)
	}

	tTau := time.Now()
	tau1, tau2, tau3 := rtype.Zero(), rtype.Zero(), rtype.Zero()
	if m > 0 {
		DvMat := mpc_core.InitRMat(rtype.Zero(), m, nProbes)
		prng := rand.New(rand.NewSource(int64(b)*1000003 + 1))
		for p := 0; p < nProbes; p++ {
			for j := 0; j < m; j++ {
				if prng.Intn(2) == 0 {
					DvMat[j][p] = w[j].Copy()
				} else {
					DvMat[j][p] = w[j].Copy().Neg()
				}
			}
		}
		wMat := mpc_core.InitRMat(rtype.Zero(), m, nProbes)
		for j := 0; j < m; j++ {
			for p := 0; p < nProbes; p++ {
				wMat[j][p] = w[j].Copy()
			}
		}
		MvMat := mActionMat(DvMat)
		u1 := vpmulMat(vmulMat(wMat, MvMat), 0.5) // Kv
		Du1 := vmulMat(wMat, u1)
		Mu1Mat := mActionMat(Du1)
		u2 := vpmulMat(vmulMat(wMat, Mu1Mat), 0.5) // K²v
		inv := 1.0 / float64(nProbes)
		tau1 = pmul(sumAllMat(vmulMat(DvMat, MvMat)), 0.5*inv)
		tau2 = pmul(sumAllMat(vmulMat(u1, u1)), inv)
		tau3 = pmul(sumAllMat(vmulMat(u1, u2)), inv)
	}
	tauSecs := time.Since(tTau).Seconds() // τ = pub Hutchinson (2 matvec passes); the dominant secure cost
	tRest := time.Now()

	// ---- private/cross corrections Δₖ, added to τₖ. Party B builds the per-gene tables in plaintext and
	// shares them: diagonals of Ψ₁,Ξ₁ (length m), Ψ₂,Ξ₂ (m×c), Π₀,Π₁,Π₂ (c×c), scalars sₖ. The one m×m
	// term tr(D_p²M_ppD_p²Ψ₁) uses Hutchinson: B forms Ψ₁·u = A₁D_v²(A₁ᵀu) locally for public probes u. ----
	// Federated Up = GᵀpubX/N (m×c); Ω' = N(XtX)⁻¹; Θ = Up·Ω'.
	upVals := [][]float64(nil)
	if haveGram {
		upVals = make([][]float64, m)
		for j := 0; j < m; j++ {
			upVals[j] = make([]float64, c)
			for l := 0; l < c; l++ {
				upVals[j][l] = lc.GtX.At(j, l) / N
			}
		}
	}
	Up := shareMat(upVals, m, c)
	NIc := mpc_core.InitRMat(rtype.Zero(), c, c)
	if hub {
		nE := rtype.FromFloat64(N, fb)
		for l := 0; l < c; l++ {
			NIc[l][l] = nE
		}
	}
	Omp := ast.choleskySolveMat(null.xtxL, null.xtxDinv, NIc) // c×c = N(XtX)⁻¹
	// m==0 (all-private gene): the public-block objects are empty; skip their zero-row SS ops (SSMultMat
	// and SSMultElemVec deref a[0]) and let only the pure-private Π/s corrections below contribute.
	Theta := mpc_core.InitRMat(rtype.Zero(), m, c) // m×c (empty when m==0)
	Dp2 := mpc_core.InitRVec(rtype.Zero(), m)      // m (empty when m==0)
	if m > 0 {
		Theta = ssMulM(Up, Omp)
		Dp2 = mpcObj.TruncVec(mpcObj.SSMultElemVec(w, w), db, fb)
	}

	// Public raw Rademacher probe signs (m×R) for the Ψ₁ Hutchinson term; byte-identical across parties.
	uSign := make([][]float64, m)
	for j := 0; j < m; j++ {
		uSign[j] = make([]float64, nProbes)
	}
	{
		prng := rand.New(rand.NewSource(int64(b)*1000003 + 7))
		for p := 0; p < nProbes; p++ {
			for j := 0; j < m; j++ {
				if prng.Intn(2) == 0 {
					uSign[j][p] = 1.0
				} else {
					uSign[j][p] = -1.0
				}
			}
		}
	}

	// Party B's reduced tables (plaintext; nil on other parties → zero share).
	var psi1diagv, xi1diagv []float64
	var psi2v, xi2v, pi0v, pi1v, pi2v, psi1Uv [][]float64
	var s1v, s2v, s3v float64
	haveB := false
	if pid > 0 && privG != nil {
		np, mp := privG.Dims()
		if mp > 0 {
			haveB = true
			dpriv := make([]float64, mp)
			for k := 0; k < mp; k++ {
				for i := 0; i < np; i++ {
					dpriv[k] += privG.At(i, k)
				}
			}
			wv := skatBetaWeight(dpriv, ast.skatTotalNumInds())
			dv2 := make([]float64, mp)
			for k := range wv {
				dv2[k] = wv[k] * wv[k]
			}
			// raw aggregates (/N): Svv (mp×mp), Uv (c×mp), A₁ (m×mp)
			svv := make([][]float64, mp)
			for k := 0; k < mp; k++ {
				svv[k] = make([]float64, mp)
				for k2 := 0; k2 < mp; k2++ {
					s := 0.0
					for i := 0; i < np; i++ {
						s += privG.At(i, k) * privG.At(i, k2)
					}
					svv[k][k2] = s / N
				}
			}
			uv := make([][]float64, c)
			for l := 0; l < c; l++ {
				uv[l] = make([]float64, mp)
				for k := 0; k < mp; k++ {
					s := 0.0
					for i := 0; i < np; i++ {
						s += X.At(i, l) * privG.At(i, k)
					}
					uv[l][k] = s / N
				}
			}
			a1 := make([][]float64, m)
			for j := 0; j < m; j++ {
				a1[j] = make([]float64, mp)
				if Gloc != nil {
					for k := 0; k < mp; k++ {
						s := 0.0
						for i := 0; i < np; i++ {
							s += Gloc.At(i, j) * privG.At(i, k)
						}
						a1[j][k] = s / N
					}
				}
			}
			// s_k = tr((D_v² S_vv)^k)
			dvSvv := make([][]float64, mp)
			for k := 0; k < mp; k++ {
				dvSvv[k] = make([]float64, mp)
				for k2 := 0; k2 < mp; k2++ {
					dvSvv[k][k2] = dv2[k] * svv[k][k2]
				}
			}
			dvSvv2 := matMul(dvSvv, dvSvv)
			s1v, s2v, s3v = matTrace(dvSvv), matTrace(dvSvv2), matTrace(matMul(dvSvv2, dvSvv))
			// Π₀,Π₁,Π₂ (c×c)
			uvD2 := matMulD(uv, dv2)                           // Uv·D²
			pi0v = matMulT(uvD2, uv)                           // Uv D² Uvᵀ
			uvD2Svv := matMul(uvD2, svv)                       // Uv D² Svv
			pi1v = matMulT(matMulD(uvD2Svv, dv2), uv)          // Uv D² Svv D² Uvᵀ
			uvD2SvvD2Svv := matMul(matMulD(uvD2Svv, dv2), svv) // Uv D² Svv D² Svv
			pi2v = matMulT(matMulD(uvD2SvvD2Svv, dv2), uv)     // ·D² Uvᵀ
			if nsnps > 0 {
				a1D2 := matMulD(a1, dv2) // A₁·D²
				// diag(Ψ₁)[j] = Σ_k A₁[j][k]² D²[k]
				psi1diagv = make([]float64, m)
				for j := 0; j < m; j++ {
					s := 0.0
					for k := 0; k < mp; k++ {
						s += a1D2[j][k] * a1[j][k]
					}
					psi1diagv[j] = s
				}
				// diag(Ξ₁)[j] = A₁[j]·(D² Svv D²)·A₁[j]ᵀ
				bmat := make([][]float64, mp) // D² Svv D²
				for k := 0; k < mp; k++ {
					bmat[k] = make([]float64, mp)
					for k2 := 0; k2 < mp; k2++ {
						bmat[k][k2] = dv2[k] * svv[k][k2] * dv2[k2]
					}
				}
				a1B := matMul(a1, bmat) // m×mp
				xi1diagv = make([]float64, m)
				for j := 0; j < m; j++ {
					s := 0.0
					for k := 0; k < mp; k++ {
						s += a1[j][k] * a1B[j][k]
					}
					xi1diagv[j] = s
				}
				psi2v = matMulT(a1D2, uv) // Ψ₂ = A₁ D² Uvᵀ (m×c)
				a1D2Svv := matMul(a1D2, svv)
				xi2v = matMulT(matMulD(a1D2Svv, dv2), uv) // Ξ₂ = A₁ D² Svv D² Uvᵀ (m×c)
				// Ψ₁·U = A₁ D²(A₁ᵀU) (m×R), computed in B's plaintext (A₁ᵀU is mp×R).
				a1T := make([][]float64, mp)
				for k := 0; k < mp; k++ {
					a1T[k] = make([]float64, m)
					for j := 0; j < m; j++ {
						a1T[k][j] = a1[j][k]
					}
				}
				a1tUd := matMul(a1T, uSign) // A₁ᵀU (mp×R)
				for k := 0; k < mp; k++ {
					for p := 0; p < nProbes; p++ {
						a1tUd[k][p] *= dv2[k]
					}
				}
				psi1Uv = matMul(a1, a1tUd) // m×R
			}
		}
	}
	shareVec := func(vals []float64, n int) mpc_core.RVec {
		out := mpc_core.InitRVec(rtype.Zero(), n)
		if vals != nil {
			for i := 0; i < n; i++ {
				out[i] = rtype.FromFloat64(vals[i], fb)
			}
		}
		return out
	}
	Psi1diag := shareVec(psi1diagv, m)
	Xi1diag := shareVec(xi1diagv, m)
	Psi2 := shareMat(psi2v, m, c)
	Xi2 := shareMat(xi2v, m, c)
	Pi0 := shareMat(pi0v, c, c)
	Pi1 := shareMat(pi1v, c, c)
	Pi2 := shareMat(pi2v, c, c)
	Psi1U := shareMat(psi1Uv, m, nProbes)
	s1 := shareScalar(s1v, haveB)
	s2 := shareScalar(s2v, haveB)
	s3 := shareScalar(s3v, haveB)

	// rowdot(A,B)[j] = Σ_l A[j][l]B[j][l]  (m×c inputs → m vector)
	rowdot := func(A, B mpc_core.RMat) mpc_core.RVec {
		prod := elemMulM(A, B)
		out := make(mpc_core.RVec, len(prod))
		for j := range prod {
			acc := rtype.Zero()
			for l := range prod[j] {
				acc = acc.Add(prod[j][l])
			}
			out[j] = acc
		}
		return out
	}
	OmPi0 := ssMulM(Omp, Pi0) // c×c = Ω'Π₀

	// tr(D_p²C) from diag(C): Ψ₁diag − 2·rowdot(Ψ₂,Θ) + rowdot(ΘΠ₀,Θ)
	trDp2C := rtype.Zero()
	trDp2C2 := rtype.Zero()
	if m > 0 {
		rd1 := rowdot(Psi2, Theta)
		rd2 := rowdot(ssMulM(Theta, Pi0), Theta)
		diagC := make(mpc_core.RVec, m)
		for j := 0; j < m; j++ {
			diagC[j] = Psi1diag[j].Sub(rd1[j]).Sub(rd1[j]).Add(rd2[j])
		}
		trDp2C = vdot(Dp2, diagC)
		// diag(C₂): (Ξ₁diag − rowdot(Ψ₂Ω',Ψ₂)) − 2·rowdot(termB,Θ) + rowdot(ΘtermD,Θ)
		termB := subMat(Xi2, ssMulM(Psi2, OmPi0)) // m×c
		termD := subMat(Pi1, ssMulM(Pi0, OmPi0))  // c×c
		rdA := rowdot(ssMulM(Psi2, Omp), Psi2)
		rdB := rowdot(termB, Theta)
		rdD := rowdot(ssMulM(Theta, termD), Theta)
		diagC2 := make(mpc_core.RVec, m)
		for j := 0; j < m; j++ {
			diagC2[j] = Xi1diag[j].Sub(rdA[j]).Sub(rdB[j]).Sub(rdB[j]).Add(rdD[j])
		}
		trDp2C2 = vdot(Dp2, diagC2)
	}

	// tr(D_p² M_pp D_p² C) = tr(D_p²M_ppD_p²Ψ₁) − 2·quad(Ψ₂) + tr(quad(Θ)·Π₀), M_pp via mActionMat.
	trDp2MppDp2C := rtype.Zero()
	if m > 0 {
		// Ψ₁ Hutchinson: (1/R) Σ_jp uSign·(D_p²·M_pp·(D_p²·Ψ₁U))
		DB1 := scaleRows(Dp2, Psi1U)
		DMDB1 := scaleRows(Dp2, mActionMat(DB1))
		psiAcc := rtype.Zero()
		for j := 0; j < m; j++ {
			for p := 0; p < nProbes; p++ {
				if uSign[j][p] > 0 {
					psiAcc = psiAcc.Add(DMDB1[j][p])
				} else {
					psiAcc = psiAcc.Sub(DMDB1[j][p])
				}
			}
		}
		trPsi1 := pmul(psiAcc, 1.0/float64(nProbes))
		// quad(Ψ₂) = tr(Θᵀ·D_p²M_ppD_p²·Ψ₂) = sumAll(Θ ⊙ D_p²M_pp(D_p²Ψ₂))
		DMV2 := scaleRows(Dp2, mActionMat(scaleRows(Dp2, Psi2)))
		quadPsi2 := sumAllMat(elemMulM(Theta, DMV2))
		// quad(Θ) = Θᵀ·D_p²M_ppD_p²·Θ (c×c); tr(quad(Θ)·Π₀)
		DMVt := scaleRows(Dp2, mActionMat(scaleRows(Dp2, Theta)))
		quadTh := ssMulM(transpose(Theta), DMVt) // c×c
		trDp2MppDp2C = trPsi1.Sub(pmul(quadPsi2, 2.0)).Add(trProd(quadTh, Pi0))
	}

	// pure-private traces (c×c / scalar): tr(Ω'Πₖ), tr((Ω'Π₀)²), tr((Ω'Π₀)³), tr(Π₁Ω'Π₀Ω')
	trOmPi0 := trProd(Omp, Pi0)
	trOmPi1 := trProd(Omp, Pi1)
	trOmPi2 := trProd(Omp, Pi2)
	trM0sq := trProd(OmPi0, OmPi0)
	trM0cu := trProd(ssMulM(OmPi0, OmPi0), OmPi0)
	trP1OP0O := trProd(ssMulM(Pi1, Omp), ssMulM(Pi0, Omp))

	// Δ₁ = ½(s₁ − tr(Ω'Π₀))
	d1 := pmul(s1.Sub(trOmPi0), 0.5)
	// Δ₂ = ½ tr(D_p²C) + ¼(s₂ − 2tr(Ω'Π₁) + tr((Ω'Π₀)²))
	d2 := pmul(trDp2C, 0.5).Add(pmul(s2.Sub(trOmPi1).Sub(trOmPi1).Add(trM0sq), 0.25))
	// Δ₃ = 3/8 tr(D_p²M_ppD_p²C) + 3/8 tr(D_p²C₂) + ⅛(s₃ − 3tr(Ω'Π₂) + 3tr(Π₁Ω'Π₀Ω') − tr((Ω'Π₀)³))
	pureP := s3.Sub(pmul(trOmPi2, 3.0)).Add(pmul(trP1OP0O, 3.0)).Sub(trM0cu)
	d3 := pmul(trDp2MppDp2C, 3.0/8.0).Add(pmul(trDp2C2, 3.0/8.0)).Add(pmul(pureP, 1.0/8.0))

	S1 = tau1.Add(d1)
	S2 = tau2.Add(d2)
	S3 = tau3.Add(d3)
	if hub {
		log.LLvl1(fmt.Sprintf("[skat_fed]   moments split: τ-hutch(2 pass) %.0fs | private+combine(Ψ₁ pass) %.0fs",
			tauSecs, time.Since(tRest).Seconds()))
	}
	return
}

// skatZSS assembles the Wilson-Hilferty pivot z from (Q, S1, S2, S3) in secret shares (δ=0, S4 unused).
// Uses the ORIGINAL Liu form (s1=S3/S2^1.5, l=1/s1², t=(Q−S1)√(S2)⁻¹/s1+l, arg=t/l=1+u, h=2/(9l)):
//
//	z = (∛(arg) − 1 + h)/√h.
//
// 1/x via (1/√x)²; √S2,√h via SqrtAndSqrtInverse; ∛ via secureCbrt. Reveal z ⇒ p = ½erfc(z/√2).
func (ast *AssocTest) skatZSSVec(Q, S1, S2, S3 mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	hub := mpcObj.GetPid() == mpcObj.GetHubPid()
	mul := func(a, b mpc_core.RVec) mpc_core.RVec {
		return mpcObj.TruncVec(mpcObj.SSMultElemVec(a, b), db, fb)
	}
	pmul := func(a mpc_core.RVec, cf float64) mpc_core.RVec {
		cfE := rtype.FromFloat64(cf, fb)
		out := make(mpc_core.RVec, len(a))
		for i := range a {
			out[i] = a[i].Mul(cfE)
		}
		return mpcObj.TruncVec(out, db, fb)
	}
	inv := func(x mpc_core.RVec) mpc_core.RVec { // 1/x = (1/√x)²
		_, si := mpcObj.SqrtAndSqrtInverse(x, false)
		return mul(si, si)
	}

	// Floor S2 off zero: an all-common window underflows S2/N² to 0 → SqrtAndSqrtInverse(0) garbage.
	// Real rare genes sit far above 1e-8, so this only catches the degenerate/underflow case.
	S2 = ast.secureClampVec(S2, 1e-8, 1e12)

	// Original Liu form (NOT the collapsed 1+u): avoids the huge S2², S2³ and their tiny inverses.
	// Interleaving the huge S3 with 1/√S2 one factor at a time keeps every intermediate O(1) — no g³
	// underflow, no S2³ overflow past mpc_data_bits. δ=0 ⇒ l=1/s1², k=l, arg = t/k = 1+u.
	_, si := mpcObj.SqrtAndSqrtInverse(S2, false) // 1/√S2
	s1 := mul(mul(mul(S3, si), si), si)           // S3/S2^1.5 (interleaved: huge→O(1))
	invS1 := inv(s1)
	l := mul(invS1, invS1) // 1/s1²
	qmS1 := make(mpc_core.RVec, len(Q))
	for i := range Q {
		qmS1[i] = Q[i].Sub(S1[i])
	}
	t := mul(mul(qmS1, si), invS1) // (Q−S1)·(1/√S2)·(1/s1)
	for i := range t {
		t[i] = t[i].Add(l[i]) // + l
	}
	invL := inv(l)
	arg := mul(t, invL)                     // t/k = t/l = 1+u
	arg = ast.secureClampVec(arg, 0.1, 9.0) // into secureCbrt's convergence window (no Newton divergence)
	h := pmul(invL, 2.0/9.0)                // 2/(9l)
	cr := ast.secureCbrtVec(arg)
	_, invSqrtH := mpcObj.SqrtAndSqrtInverse(h, false)
	num := make(mpc_core.RVec, len(cr))
	for i := range cr {
		num[i] = cr[i].Add(h[i]) // ∛(t/k) + h
		if hub {
			num[i] = num[i].Sub(rtype.FromFloat64(1.0, fb)) // − 1
		}
	}
	return mul(num, invSqrtH) // (∛(t/k) − 1 + h)/√h
}

func (ast *AssocTest) skatZSS(Q, S1, S2, S3 mpc_core.RElem) mpc_core.RElem {
	return ast.skatZSSVec(mpc_core.RVec{Q}, mpc_core.RVec{S1}, mpc_core.RVec{S2}, mpc_core.RVec{S3})[0]
}

// skatPValueSS returns the SKAT WH pivot z (SS) for gene b's public list: Q (scaled SKAT statistic) +
// Hutchinson moments + skatZSS. Caller reveals z, then p = ½erfc(z/√2). Public-list only (no private yet).
func (ast *AssocTest) skatPValueSS(b, nsnps, nProbes int, null skatNull, nullRSS crypto.CipherVector, X *mat.Dense, y0 []float64, privG *mat.Dense, privatePid int) mpc_core.RElem {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	// full-gene SKAT statistic Q = Σŵ²s² over public list (PART A) + private variants (PART B)
	qPub, _, _ := ast.blockStat(b, nsnps, null, X, y0, nil)
	qPriv, _ := ast.privateBlockStat(privG, null, X, y0, privatePid)
	qRaw := mpc_core.RVec{qPub[0].Add(qPriv[0])}
	scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS)
	if !ok {
		panic("skatPValueSS: scale undefined (dof ≤ 0)")
	}
	Q := mpcObj.TruncVec(mpcObj.SSMultElemVec(qRaw, mpc_core.RVec{scaleSS[0]}), db, fb)[0] // Q/(2σ̂²)
	// Normalize Q by N to match the N-normalized moments (kernel /N); WH ratios are scale-invariant.
	invN := rtype.FromFloat64(1.0/float64(ast.skatTotalNumInds()), fb)
	Q = mpcObj.TruncVec(mpc_core.RVec{Q.Mul(invN)}, db, fb)[0]
	S1, S2, S3 := ast.skatMomentsSS(b, nsnps, nProbes, null, X, y0, privG, nil)
	return ast.skatZSS(Q, S1, S2, S3)
}
