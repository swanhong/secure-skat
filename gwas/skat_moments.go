package gwas

import (
	"fmt"
	"math"
	"math/rand"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

// secureClampVec returns min(max(x, loPub), hiPub) for SECRET x via two secure compares and
// branch-free selects. It is shared by the moment-domain guards and the cube-root Newton core.
func (ast *AssocTest) secureClampVec(x mpc_core.RVec, loPub, hiPub float64) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	fb := mpcObj.GetFracBits()
	bv := mpcObj.GetBooleanShareFlag()
	n := len(x)
	x = ast.secureSelectVec(mpcObj.LessThanPublic(x, rtype.FromFloat64(loPub, fb), bv), ast.hubVec(loPub, n), x)    // x<lo → lo
	x = ast.secureSelectVec(mpcObj.NotLessThanPublic(x, rtype.FromFloat64(hiPub, fb), bv), ast.hubVec(hiPub, n), x) // x≥hi → hi
	return x
}

// secureSelectVec returns cond ? whenTrue : whenFalse without revealing the scale-0 secret bit cond.
func (ast *AssocTest) secureSelectVec(cond, whenTrue, whenFalse mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	fb := mpcObj.GetFracBits()
	cond = cond.Copy()
	cond.MulScalar(rtype.FromFloat64(1.0, fb))
	diff := whenTrue.Copy()
	diff.Sub(whenFalse)
	delta := ast.ssMul(cond, diff)
	out := whenFalse.Copy()
	out.Add(delta)
	return out
}

// secureSelectOrPublicVec returns valid ? x : fallbackPub without revealing valid.
// valid is an arithmetic sharing of bits at scale 0; x/fallback use the configured fixed-point scale.
func (ast *AssocTest) secureSelectOrPublicVec(valid, x mpc_core.RVec, fallbackPub float64) mpc_core.RVec {
	return ast.secureSelectVec(valid, x, ast.hubVec(fallbackPub, len(x)))
}

// secureCbrtVec returns the real signed cube root of SECRET x. It classifies |x| against a fixed,
// public ladder of powers of eight, constructs the corresponding secret range/root multipliers,
// and applies them in a fixed number of public groups. Thus every representable nonzero input is
// moved into the Newton core's [0.1,9] interval without revealing its magnitude; zero is exact.
// Branch-free inverse-cube-root Newton on the positive reduced input:
// y → x^(-1/3) via y ← (4y − x·y⁴)/3 from seed 0.7, then x^(1/3) = x·y². The rearranged
// update is identical to y(4−xy³)/3 but saves one secret multiplication per iteration.
func (ast *AssocTest) secureCbrtVec(x mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	n := len(x)

	mul, square, pmul := ast.ssMul, ast.ssSquare, ast.ssPMul
	addScaledBits := func(dst, bits mpc_core.RVec, cf float64) {
		cfE := rtype.FromFloat64(cf, fb)
		for i := range dst {
			dst[i] = dst[i].Add(bits[i].Mul(cfE)) // scale-0 bit × public fixed-point constant
		}
	}

	bv := mpcObj.GetBooleanShareFlag()
	neg := mpcObj.LessThanPublic(x, rtype.Zero(), bv)
	negX := x.Copy()
	for i := range negX {
		negX[i] = negX[i].Neg()
	}
	absX := ast.secureSelectVec(neg, negX, x)

	// If n nested ladder predicates hold, these linear identities give the secret factors:
	//   8^-n = 1 - Σ 7/8^(k+1),  2^n = 1 + Σ 2^k,
	//   8^n  = 1 + Σ 7·8^k,      2^-n = 1 - Σ 2^-(k+1).
	// A group is capped so its smallest fractional and largest integer coefficient are representable.
	// Repeating public groups supports asymmetric db/fb configurations.
	integerBits := db - fb - 1
	if integerBits < 4 || fb < 3 {
		panic("secureCbrtVec: fixed-point format needs at least 4 integer and 3 fractional bits")
	}
	ceilPositive := func(v float64) int {
		if v <= 0 {
			return 0
		}
		return int(math.Ceil(v))
	}
	highSteps := ceilPositive((float64(integerBits) - math.Log2(9.0)) / 3.0)
	lowSteps := ceilPositive((float64(fb) + math.Log2(0.1)) / 3.0)
	rootScales := make([]mpc_core.RVec, 0, 2)
	reduce := func(steps, maxGroup int, high bool) {
		for steps > 0 {
			group := steps
			if group > maxGroup {
				group = maxGroup
			}
			rangeScale, rootScale := ast.hubVec(1.0, n), ast.hubVec(1.0, n)
			highThreshold, highRootCoeff := 8.0, 1.0
			lowThreshold, lowRangeCoeff, lowRootCoeff := 0.1, 7.0, 0.5
			for k := 0; k < group; k++ {
				if high {
					bit := mpcObj.NotLessThanPublic(absX, rtype.FromFloat64(highThreshold, fb), bv)
					addScaledBits(rangeScale, bit, -7.0/highThreshold)
					addScaledBits(rootScale, bit, highRootCoeff)
					highThreshold *= 8.0
					highRootCoeff *= 2.0
				} else {
					bit := mpcObj.LessThanPublic(absX, rtype.FromFloat64(lowThreshold, fb), bv)
					addScaledBits(rangeScale, bit, lowRangeCoeff)
					addScaledBits(rootScale, bit, -lowRootCoeff)
					lowThreshold /= 8.0
					lowRangeCoeff *= 8.0
					lowRootCoeff *= 0.5
				}
			}
			absX = mul(absX, rangeScale)
			rootScales = append(rootScales, rootScale)
			steps -= group
		}
	}
	reduce(highSteps, fb/3, true)              // 8^-group must remain representable
	reduce(lowSteps, (integerBits-1)/3, false) // 8^group must fit the signed integer range
	// After the full public ladder, every representable nonzero magnitude is >=0.1; only zero is tiny.
	tiny := mpcObj.LessThanPublic(absX, rtype.FromFloat64(0.1, fb), bv)
	x = ast.secureSelectVec(tiny, ast.hubVec(0.1, n), absX)

	y := ast.hubVec(0.7, n) // seed 0.7 as an additive share (hub holds it, others 0)
	four := rtype.FromFloat64(4.0, 0)
	for it := 0; it < 8; it++ {
		y2 := square(y)
		y4 := square(y2)
		xy4 := mul(x, y4)
		next := make(mpc_core.RVec, n)
		for i := range next {
			next[i] = y[i].Mul(four).Sub(xy4[i]) // 4y − x·y⁴
		}
		y = pmul(next, 1.0/3.0)
	}
	root := mul(x, square(y)) // x·y² = x^(1/3) on the reduced interval
	for _, rootScale := range rootScales {
		root = mul(root, rootScale)
	}
	// signNonzero is +1 (positive), -1 (negative), or 0 (zero); neg and tiny are disjoint.
	signNonzero := ast.hubVec(1.0, n)
	addScaledBits(signNonzero, neg, -2.0)
	addScaledBits(signNonzero, tiny, -1.0)
	return mul(root, signNonzero)
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

// matTraceProduct returns tr(a*b) without materializing the full product.
// It is especially useful for tr(A^3)=tr((A^2)A), where A^2 is already available.
func matTraceProduct(a, b [][]float64) float64 {
	s := 0.0
	for i := range a {
		for j := range a[i] {
			s += a[i][j] * b[j][i]
		}
	}
	return s
}

// skatTraceProbes uses the standard basis for an exact trace when requested >= m and deterministic
// Rademacher probes otherwise. The branch depends only on public dimensions.
func skatTraceProbes(m, requested int, seed int64) (values [][]float64, multiplier float64, exact bool) {
	if requested <= 0 {
		panic("skatTraceProbes: requested must be positive")
	}
	if m == 0 {
		return nil, 1.0, true
	}
	n := requested
	if n >= m {
		n = m
		exact = true
	}
	values = make([][]float64, m)
	for j := range values {
		values[j] = make([]float64, n)
	}
	if exact {
		for j := 0; j < m; j++ {
			values[j][j] = 1.0
		}
		return values, 1.0, true
	}
	prng := rand.New(rand.NewSource(seed))
	for p := 0; p < n; p++ {
		for j := 0; j < m; j++ {
			if prng.Intn(2) == 0 {
				values[j][p] = 1.0
			} else {
				values[j][p] = -1.0
			}
		}
	}
	return values, 1.0 / float64(n), false
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
// moment τ₁=tr(K_pp) is exact; τ₂,τ₃ use either an exact public basis (m<=requested probes) or
// Hutchinson probes. Private/cross terms come from per-gene tables (Ψ,Ξ,Π,s) party B builds in
// plaintext and shares. Moments
// are N-normalized (Sₖ/Nᵏ); w is the unsigned SKAT weight; S₄ is unneeded since K is PSD.
func (ast *AssocTest) skatMomentsSS(b, nsnps, nProbes int, null skatNull, priv *privateGeneLocal, gl *geneLocal, w mpc_core.RVec) (S1, S2, S3 mpc_core.RElem) {
	if nProbes <= 0 {
		panic("skatMomentsSS: nProbes must be positive")
	}
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	pid := mpcObj.GetPid()
	hub := pid == mpcObj.GetHubPid()
	c := null.c
	m := nsnps // public block only; the private variants enter via the contracted corrections below
	setupMark := ast.metricMark()

	// Kernel normalized by N (public) so SS intermediates stay O(1) and the WH ratios are unchanged.
	// Both Hutchinson and the private combine reuse Uₚ=GᵀX/N, Ω'=(XᵀX/N)⁻¹, and Θ=UₚΩ'.
	N := float64(ast.skatTotalNumInds())

	// ---- SS helpers ----
	pmul := func(a mpc_core.RElem, cf float64) mpc_core.RElem { // secret × public const → truncated
		if cf == 1 {
			return a.Copy()
		}
		if cf == math.Trunc(cf) {
			return a.Mul(rtype.FromInt(int(cf))) // integer coefficient preserves fixed-point scale
		}
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
	sumSquaresMat := func(a mpc_core.RMat) mpc_core.RElem {
		flat := make(mpc_core.RVec, 0, len(a)*len(a[0]))
		for i := range a {
			flat = append(flat, a[i]...)
		}
		sq := ast.ssSquare(flat)
		acc := rtype.Zero()
		for i := range sq {
			acc = acc.Add(sq[i])
		}
		return acc
	}
	traceSS := func(a mpc_core.RMat) mpc_core.RElem {
		acc := rtype.Zero()
		for i := range a {
			acc = acc.Add(a[i][i])
		}
		return acc
	}
	trProd := func(a, bb mpc_core.RMat) mpc_core.RElem { // tr(A·B) = Σ_ij A[i][j]B[j][i] = sumAll(A⊙Bᵀ)
		if len(a) == 0 || len(bb) == 0 {
			return rtype.Zero()
		}
		return sumAllMat(elemMulM(a, transpose(bb)))
	}
	vdot := ast.ssDot
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

	// ---- public-block gram/coupling. gg=GᵀpubGpub stays local plaintext, while Up=GᵀpubX/N is
	// entered once as an additive share; the private cache already holds the local A₁ contraction. ----
	var gg *mat.Dense
	Up := mpc_core.InitRMat(rtype.Zero(), m, c)
	sppDiag := mpc_core.InitRVec(rtype.Zero(), m)
	if pid > 0 && nsnps > 0 {
		gg = gl.gg
		for j := 0; j < m; j++ {
			sppDiag[j] = rtype.FromFloat64(gg.At(j, j)/N, fb)
			for l := 0; l < c; l++ {
				Up[j][l] = rtype.FromFloat64(gl.GtX.At(j, l)/N, fb)
			}
		}
	}
	tauProbe, tauMultiplier, exactPublic := skatTraceProbes(m, nProbes, int64(b)*1000003+1)
	psiProbe, psiMultiplier := tauProbe, tauMultiplier
	if !exactPublic {
		psiProbe, psiMultiplier, _ = skatTraceProbes(m, nProbes, int64(b)*1000003+7)
	}
	probeCount := 0
	if m > 0 {
		probeCount = len(tauProbe[0])
	}
	Omp := null.omp // c×c = N(XtX)⁻¹; cached once by the null model
	// m==0 (all-private gene): public objects are empty; only the pure-private Π/s terms contribute.
	Theta := mpc_core.InitRMat(rtype.Zero(), m, c)
	Dp2 := mpc_core.InitRVec(rtype.Zero(), m)
	if m > 0 {
		Theta = ssMulM(Up, Omp)
		Dp2 = ast.ssSquare(w)
	}
	upT := transpose(Up)

	// ---- τ₁ exact; τ₂,τ₃ exact-basis or Hutchinson over public probes (M_pp·A via chunked GᵀG). ----
	// Hutchinson-only progress over secure M actions. Exact-basis mode forms K directly and skips this
	// chunked path, so it has no synthetic percentage counter.
	nChunks := (m + gtgChunkRows - 1) / gtgChunkRows
	momTotal := 0.0
	if !exactPublic {
		momTotal = float64(nChunks * (3*probeCount + 2*c))
	}
	momDone := 0.0
	momStart := time.Now()
	momPct := 0
	mActionMat := func(A mpc_core.RMat) mpc_core.RMat { // M_pp·A = (GᵀG/N)A − Θ(UₚᵀA)
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
		// Θ and Uₚ already contain the one per-gene inverse solve. Reusing them here avoids solving
		// (XtX) against every probe/correction batch while computing the identical projection term.
		upTA := ssMulM(upT, A)
		gxsol := ssMulM(Theta, upTA)
		Mv := mpc_core.InitRMat(rtype.Zero(), m, R)
		for j := 0; j < m; j++ {
			for p := 0; p < R; p++ {
				Mv[j][p] = ga[j][p].Sub(gxsol[j][p])
			}
		}
		return Mv
	}
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

	ast.metricEnd("gene_moments_setup_gtx", setupMark)
	publicTraceMark := ast.metricMark()
	tTau := time.Now()
	tau1, tau2, tau3 := rtype.Zero(), rtype.Zero(), rtype.Zero()
	var Kexact mpc_core.RMat
	if m > 0 {
		if exactPublic {
			// The probe basis is I_m. Avoid the generic first M_pp·diag(w) dense multiply: form
			// K=½D(S_pp−ΘU_pᵀ)D with O(m²c)+elementwise work, then use one dense K² multiply.
			Spp := mpc_core.InitRMat(rtype.Zero(), m, m)
			for j := 0; j < m; j++ {
				if gg != nil {
					for k := 0; k < m; k++ {
						Spp[j][k] = rtype.FromFloat64(gg.At(j, k)/N, fb)
					}
				}
			}
			Mpp := subMat(Spp, ssMulM(Theta, upT))
			wCols := mpc_core.InitRMat(rtype.Zero(), m, m)
			for j := 0; j < m; j++ {
				for k := 0; k < m; k++ {
					wCols[j][k] = w[k].Copy()
				}
			}
			Kexact = vpmulMat(elemMulM(scaleRows(w, Mpp), wCols), 0.5)
			K2 := ssMulM(Kexact, Kexact)
			tau1 = traceSS(Kexact)
			tau2 = traceSS(K2)
			tau3 = trProd(K2, Kexact)
		} else {
			DvMat := mpc_core.InitRMat(rtype.Zero(), m, probeCount)
			for p := 0; p < probeCount; p++ {
				for j := 0; j < m; j++ {
					if tauProbe[j][p] > 0 {
						DvMat[j][p] = w[j].Copy()
					} else {
						DvMat[j][p] = w[j].Copy().Neg()
					}
				}
			}
			wMat := mpc_core.InitRMat(rtype.Zero(), m, probeCount)
			for j := 0; j < m; j++ {
				for p := 0; p < probeCount; p++ {
					wMat[j][p] = w[j].Copy()
				}
			}
			MvMat := mActionMat(DvMat)
			u1 := vpmulMat(elemMulM(wMat, MvMat), 0.5) // Kv
			Du1 := elemMulM(wMat, u1)
			Mu1Mat := mActionMat(Du1)
			u2 := vpmulMat(elemMulM(wMat, Mu1Mat), 0.5) // K²v
			tau2 = pmul(sumSquaresMat(u1), tauMultiplier)
			tau3 = pmul(sumAllMat(elemMulM(u1, u2)), tauMultiplier)
		}
	}
	tauSecs := time.Since(tTau).Seconds() // public τ₂/τ₃ phase (explicit K or two-pass Hutchinson)
	ast.metricEnd("gene_moments_public_trace_gtg_gtx", publicTraceMark)
	privateCrossMark := ast.metricMark()
	tRest := time.Now()

	// ---- private/cross corrections Δₖ, added to τₖ. Party B builds the per-gene tables in plaintext and
	// shares them: diagonals of Ψ₁,Ξ₁ (length m), Ψ₂,Ξ₂ (m×c), Π₀,Π₁,Π₂ (c×c), scalars sₖ. The one m×m
	// term tr(D_p²M_ppD_p²Ψ₁) uses the same exact/Hutchinson rule: B forms
	// Ψ₁·u = A₁D_v²(A₁ᵀu) locally for public probes u. ----
	// Party B's reduced tables (plaintext; nil on other parties → zero share).
	var psi1diagv, xi1diagv []float64
	var psi2v, xi2v, pi0v, pi1v, pi2v, psi1Uv [][]float64
	var s1v, s2v, s3v float64
	haveB := false
	if pid > 0 && priv != nil {
		mp := len(priv.Weight)
		if mp > 0 {
			haveB = true
			wv := priv.Weight
			dv2 := make([]float64, mp)
			for k := range wv {
				dv2[k] = wv[k] * wv[k]
			}
			// Reuse the cached /N contractions Svv=G_vᵀG_v/N and A₁=G_pᵀG_v/N.
			svv := priv.Svv
			if len(svv) != mp || (m > 0 && len(priv.A1) != m) {
				panic("skatMomentsSS: incomplete private gene cache")
			}
			uv := make([][]float64, c)
			for l := 0; l < c; l++ {
				uv[l] = make([]float64, mp)
				for k := 0; k < mp; k++ {
					uv[l][k] = priv.GtX.At(k, l) / N
				}
			}
			a1 := priv.A1
			// s_k = tr((D_v² S_vv)^k)
			dvSvv := make([][]float64, mp)
			for k := 0; k < mp; k++ {
				dvSvv[k] = make([]float64, mp)
				for k2 := 0; k2 < mp; k2++ {
					dvSvv[k][k2] = dv2[k] * svv[k][k2]
				}
			}
			dvSvv2 := matMul(dvSvv, dvSvv)
			s1v, s2v, s3v = matTrace(dvSvv), matTrace(dvSvv2), matTraceProduct(dvSvv2, dvSvv)
			// Π₀,Π₁,Π₂ (c×c)
			uvD2 := matMulD(uv, dv2)                       // Uv·D²
			pi0v = matMulT(uvD2, uv)                       // Uv D² Uvᵀ
			uvD2Svv := matMul(uvD2, svv)                   // Uv D² Svv
			uvD2SvvD2 := matMulD(uvD2Svv, dv2)             // Uv D² Svv D²
			pi1v = matMulT(uvD2SvvD2, uv)                  // Uv D² Svv D² Uvᵀ
			uvD2SvvD2Svv := matMul(uvD2SvvD2, svv)         // Uv D² Svv D² Svv
			pi2v = matMulT(matMulD(uvD2SvvD2Svv, dv2), uv) // ·D² Uvᵀ
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
				a1D2Svv := matMul(a1D2, svv)
				a1B := matMulD(a1D2Svv, dv2) // A₁ D² Svv D² (m×mp)
				xi1diagv = make([]float64, m)
				for j := 0; j < m; j++ {
					s := 0.0
					for k := 0; k < mp; k++ {
						s += a1[j][k] * a1B[j][k]
					}
					xi1diagv[j] = s
				}
				psi2v = matMulT(a1D2, uv) // Ψ₂ = A₁ D² Uvᵀ (m×c)
				xi2v = matMulT(a1B, uv)   // Ξ₂ = A₁ D² Svv D² Uvᵀ (m×c); reuse the matrix above
				// Ψ₁·U = A₁ D²(A₁ᵀU) (m×R), computed in B's plaintext (A₁ᵀU is mp×R).
				a1T := make([][]float64, mp)
				for k := 0; k < mp; k++ {
					a1T[k] = make([]float64, m)
					for j := 0; j < m; j++ {
						a1T[k][j] = a1[j][k]
					}
				}
				a1tUd := matMul(a1T, psiProbe) // A₁ᵀU (mp×probeCount)
				for k := 0; k < mp; k++ {
					for p := 0; p < probeCount; p++ {
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
	Psi1U := shareMat(psi1Uv, m, probeCount)
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
		// Hutchinson corrections share one wide M action; exact-basis corrections reuse explicit K.
		basis := mpc_core.InitRMat(rtype.Zero(), m, probeCount+2*c)
		for j := 0; j < m; j++ {
			for p := 0; p < probeCount; p++ {
				basis[j][p] = Psi1U[j][p].Copy()
			}
			for l := 0; l < c; l++ {
				basis[j][probeCount+l] = Psi2[j][l].Copy()
				basis[j][probeCount+c+l] = Theta[j][l].Copy()
			}
		}
		var acted mpc_core.RMat
		if exactPublic {
			// D²M_ppD²B = 2·D·K·D·B, reusing the explicit exact-basis K.
			acted = scaleRows(w, ssMulM(Kexact, scaleRows(w, basis)))
			for j := range acted {
				acted[j].MulScalar(rtype.FromInt(2))
			}
		} else {
			weightedBasis := scaleRows(Dp2, basis)
			// τ₁ = ½[Σ_j d_j²Spp_jj − Σ_jl(d_j²Θ_jl)Up_jl] exactly.
			tauLeft := make(mpc_core.RVec, m*(c+1))
			tauRight := make(mpc_core.RVec, m*(c+1))
			for j := 0; j < m; j++ {
				tauLeft[j], tauRight[j] = Dp2[j], sppDiag[j]
				for l := 0; l < c; l++ {
					idx := m + j*c + l
					tauLeft[idx] = weightedBasis[j][probeCount+c+l]
					tauRight[idx] = Up[j][l]
				}
			}
			tauTerms := ast.ssMul(tauLeft, tauRight)
			tau1Acc := rtype.Zero()
			for i := 0; i < m; i++ {
				tau1Acc = tau1Acc.Add(tauTerms[i])
			}
			for i := m; i < len(tauTerms); i++ {
				tau1Acc = tau1Acc.Sub(tauTerms[i])
			}
			tau1 = pmul(tau1Acc, 0.5)
			acted = scaleRows(Dp2, mActionMat(weightedBasis))
		}

		// Ψ₁ trace: basis sum if exact, otherwise (1/R)Σ uᵀD_p²M_ppD_p²Ψ₁u.
		psiAcc := rtype.Zero()
		for j := 0; j < m; j++ {
			for p := 0; p < probeCount; p++ {
				if psiProbe[j][p] > 0 {
					psiAcc = psiAcc.Add(acted[j][p])
				} else if psiProbe[j][p] < 0 {
					psiAcc = psiAcc.Sub(acted[j][p])
				}
			}
		}
		trPsi1 := pmul(psiAcc, psiMultiplier)
		// quad(Ψ₂) = tr(Θᵀ·D_p²M_ppD_p²·Ψ₂) = sumAll(Θ ⊙ D_p²M_pp(D_p²Ψ₂))
		DMV2 := mpc_core.InitRMat(rtype.Zero(), m, c)
		DMVt := mpc_core.InitRMat(rtype.Zero(), m, c)
		for j := 0; j < m; j++ {
			for l := 0; l < c; l++ {
				DMV2[j][l] = acted[j][probeCount+l]
				DMVt[j][l] = acted[j][probeCount+c+l]
			}
		}
		quadPsi2 := sumAllMat(elemMulM(Theta, DMV2))
		// quad(Θ) = Θᵀ·D_p²M_ppD_p²·Θ (c×c); tr(quad(Θ)·Π₀)
		quadTh := ssMulM(transpose(Theta), DMVt) // c×c
		trDp2MppDp2C = trPsi1.Sub(pmul(quadPsi2, 2.0)).Add(trProd(quadTh, Pi0))
	}

	// pure-private traces (c×c / scalar): tr(Ω'Πₖ), tr((Ω'Π₀)²), tr((Ω'Π₀)³), tr(Π₁Ω'Π₀Ω')
	OmPi1 := ssMulM(Omp, Pi1)
	M0sq := ssMulM(OmPi0, OmPi0)
	trOmPi0 := traceSS(OmPi0)
	trOmPi1 := traceSS(OmPi1)
	trOmPi2 := trProd(Omp, Pi2)
	trM0sq := traceSS(M0sq)
	trM0cu := trProd(M0sq, OmPi0)
	trP1OP0O := trProd(OmPi1, OmPi0) // cyclic: tr((OΠ₁)(OΠ₀)) = tr(Π₁OΠ₀O)

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
	ast.metricEnd("gene_moments_private_cross", privateCrossMark)
	if hub {
		mode := fmt.Sprintf("hutch(%d cols, 2 M actions)", probeCount)
		if exactPublic {
			mode = fmt.Sprintf("exact-basis(%d cols, explicit K+K²)", probeCount)
		}
		log.LLvl1(fmt.Sprintf("[skat_fed]   moments split: public-%s %.0fs | τ1/private/combine %.0fs",
			mode, tauSecs, time.Since(tRest).Seconds()))
	}
	return
}

func skatMomentFloor(fracBits int) float64 {
	return math.Ldexp(1.0, -fracBits)
}

func skatSkewFloor(fracBits int) float64 {
	// h=2*skew²/9 must also remain positive after fixed-point quantization; this gives h>=2 quanta.
	return math.Max(1e-4, 3.0*math.Sqrt(math.Ldexp(1.0, -fracBits)))
}

// skatZSSVec assembles the Wilson-Hilferty pivot z from (Q, S1, S2, S3) in secret shares (δ=0, S4 unused).
// Uses the algebraically reduced Liu form s=S3/S2^1.5, u=(Q−S1)s/√S2, h=2s²/9:
//
//	z = (∛(arg) − 1 + h)/√h.
//
// √S2,√h use SqrtAndSqrtInverse; ∛ uses secureCbrt. Reveal z ⇒ p = ½erfc(z/√2).
func (ast *AssocTest) skatZSSVec(Q, S1, S2, S3 mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	hub := mpcObj.GetPid() == mpcObj.GetHubPid()
	mul, square, pmul := ast.ssMul, ast.ssSquare, ast.ssPMul
	// Record the PSD/fixed-point domain guard before clamping. The bit remains secret; non-positive,
	// underflowed, or unrepresentably large genes use safe surrogate inputs and end at z=-9 (p≈1).
	momentFloor := skatMomentFloor(fb)
	momentCeil := math.Ldexp(1.0, db-fb-2) // one signed-integer headroom bit for sqrt/intermediates
	bv := mpcObj.GetBooleanShareFlag()
	s2Valid := mpcObj.NotLessThanPublic(S2, rtype.FromFloat64(momentFloor, fb), bv)
	s3Valid := mpcObj.NotLessThanPublic(S3, rtype.FromFloat64(momentFloor, fb), bv)
	s2BelowCeil := mpcObj.LessThanPublic(S2, rtype.FromFloat64(momentCeil, fb), bv)
	valid := mpcObj.SSMultElemVec(s2Valid, s3Valid)
	valid = mpcObj.SSMultElemVec(valid, s2BelowCeil) // scale-0 AND

	// Floor S2 off zero: an all-common window can underflow S2/N² to 0, while the normalizer-backed
	// square root requires a strictly positive input. The one-quantum floor is the smallest positive
	// value representable at fb bits; the ceil leaves one signed-integer headroom bit.
	S2 = ast.secureSelectOrPublicVec(s2Valid, S2, momentFloor)
	S2 = ast.secureSelectOrPublicVec(s2BelowCeil, S2, momentCeil)
	S3 = ast.secureSelectOrPublicVec(valid, S3, momentFloor)

	// Form the dimensionless ratios directly. The split skew path below and
	// arg=1+((Q−S1)/√S2)*(S3/S2^1.5) avoid both the S2²/S2³ intermediates and a
	// 1/s² -> 1/l inverse round trip (two extra secure inverse-square-root calls).
	_, si := mpcObj.SqrtAndSqrtInverse(S2, false) // 1/√S2
	// Compute min(S3/S2^1.5,1) without allowing an inconsistent Hutchinson pair to wrap before
	// the final clamp. For S2<1, cap after every increasing multiply by si; for S2>=1, si<=1 so
	// the three products are monotonically non-increasing. The inactive branch receives safe inputs.
	oneVec := ast.hubVec(1.0, len(S2))
	smallS2 := mpcObj.LessThanPublic(S2, rtype.FromFloat64(1.0, fb), bv)
	lowSi := ast.secureSelectOrPublicVec(smallS2, si, 1.0)
	lowSkew := ast.secureClampVec(S3, 0.0, 1.0)
	for it := 0; it < 3; it++ {
		lowSkew = mul(lowSkew, lowSi)
		lowSkew = ast.secureClampVec(lowSkew, 0.0, 1.0)
	}
	highSi := ast.secureSelectVec(smallS2, oneVec, si)
	highSkew := ast.secureSelectVec(smallS2, oneVec, S3)
	for it := 0; it < 3; it++ {
		highSkew = mul(highSkew, highSi)
	}
	skew := ast.secureSelectVec(smallS2, lowSkew, highSkew)
	// A PSD spectrum has 0 < S3/S2^(3/2) <= 1. Enforce that feasible range after noisy estimation;
	// the lower floor also keeps h and its inverse square root away from zero on invalid genes.
	skew = ast.secureClampVec(skew, skatSkewFloor(fb), 1.0)
	qmS1 := make(mpc_core.RVec, len(Q))
	for i := range Q {
		qmS1[i] = Q[i].Sub(S1[i])
	}
	arg := mul(mul(qmS1, si), skew)
	if hub {
		one := rtype.FromFloat64(1.0, fb)
		for i := range arg {
			arg[i] = arg[i].Add(one) // 1 + u
		}
	}
	h := pmul(square(skew), 2.0/9.0) // 2s²/9
	cr := ast.secureCbrtVec(arg)
	_, invSqrtH := mpcObj.SqrtAndSqrtInverse(h, false)
	num := make(mpc_core.RVec, len(cr))
	for i := range cr {
		num[i] = cr[i].Add(h[i]) // ∛(t/k) + h
		if hub {
			num[i] = num[i].Sub(rtype.FromFloat64(1.0, fb)) // − 1
		}
	}
	z := mul(num, invSqrtH) // (∛(t/k) − 1 + h)/√h
	return ast.secureSelectOrPublicVec(valid, z, -9.0)
}
