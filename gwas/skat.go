package gwas

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

// LocalContraction is one party's plaintext contraction of (G, X, y0): only c-/m-dim
// aggregates (never n), additive across parties so Σ_party == full-cohort contraction.
type LocalContraction struct {
	XtX       *mat.Dense // c×c
	Xty0      []float64  // c
	Y0ty0     float64
	GtX       *mat.Dense // m×c
	Gty0      []float64  // m
	DosageSum []float64  // m
}

func vecToSlice(v *mat.VecDense) []float64 {
	out := make([]float64, v.Len())
	for i := range out {
		out[i] = v.AtVec(i)
	}
	return out
}

// LocalContract contracts (G, X, y0) locally; y0 centered, X includes intercept.
func LocalContract(G, X *mat.Dense, y0 []float64) LocalContraction {
	n, _ := X.Dims()
	_, m := G.Dims()
	y0v := mat.NewVecDense(n, y0)

	var XtX mat.Dense
	XtX.Mul(X.T(), X)

	var Xty mat.VecDense
	Xty.MulVec(X.T(), y0v)

	var GtX mat.Dense
	GtX.Mul(G.T(), X)

	var Gty mat.VecDense
	Gty.MulVec(G.T(), y0v)

	dosage := make([]float64, m)
	for j := 0; j < m; j++ {
		for i := 0; i < n; i++ {
			dosage[j] += G.At(i, j)
		}
	}

	return LocalContraction{
		XtX:       &XtX,
		Xty0:      vecToSlice(&Xty),
		Y0ty0:     mat.Dot(y0v, y0v),
		GtX:       &GtX,
		Gty0:      vecToSlice(&Gty),
		DosageSum: dosage,
	}
}

// Add sums two contractions party-wise (all fields additive) — the n-independence invariant
// (Σ_party == full cohort) checked by TestSKATLocalPartyAdditivity.
func (a LocalContraction) Add(b LocalContraction) LocalContraction {
	var XtX, GtX mat.Dense
	XtX.Add(a.XtX, b.XtX)
	GtX.Add(a.GtX, b.GtX)
	addVec := func(x, y []float64) []float64 {
		out := make([]float64, len(x))
		for i := range x {
			out[i] = x[i] + y[i]
		}
		return out
	}
	return LocalContraction{
		XtX:       &XtX,
		Xty0:      addVec(a.Xty0, b.Xty0),
		Y0ty0:     a.Y0ty0 + b.Y0ty0,
		GtX:       &GtX,
		Gty0:      addVec(a.Gty0, b.Gty0),
		DosageSum: addVec(a.DosageSum, b.DosageSum),
	}
}

type geneLocal struct {
	LocalContraction            // the n-free aggregates (GᵀX, Gᵀy0, dosage), promoted for direct access
	Gloc             *mat.Dense // this party's aligned public genotype block (n×m — the n-dim data LocalContraction never holds)
	gg               *mat.Dense // GᵀpubGpub local Gram
}

func (ast *AssocTest) computeGeneLocal(b, nsnps int, X *mat.Dense, y0 []float64) *geneLocal {
	if ast.general.mpcObj[0].GetPid() == 0 || nsnps == 0 {
		return &geneLocal{}
	}
	Gloc := ast.readGenoBlockLocal(b)
	var gg mat.Dense
	gg.Mul(Gloc.T(), Gloc)
	return &geneLocal{LocalContraction: LocalContract(Gloc, X, y0), Gloc: Gloc, gg: &gg}
}

func (ast *AssocTest) localFor(b, nsnps int, X *mat.Dense, y0 []float64, gl *geneLocal) *geneLocal {
	if gl != nil {
		return gl
	}
	return ast.computeGeneLocal(b, nsnps, X, y0)
}

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

// --- secure low-rank null model (β̂, RSS from c-dim aggregates only; n never enters) ---

type localNull struct {
	X     *mat.Dense
	Y0    []float64
	XtX   *mat.Dense
	Xty0  []float64
	Y0ty0 float64
}

// localNullEquations builds this party's X=[1|cov], centered y0, and XᵀX/Xᵀy₀/y₀ᵀy₀.
// nil cov/pheno (pid-0) yields an empty result (zero aggregates).
func localNullEquations(cov, pheno *mat.Dense, center float64) localNull {
	if cov == nil || pheno == nil {
		return localNull{}
	}
	cn, ncov := cov.Dims()
	pn, _ := pheno.Dims()
	n := cn
	if pn < n {
		n = pn
	}
	c := ncov + 1

	X := mat.NewDense(n, c, nil)
	y0 := make([]float64, n)
	for i := 0; i < n; i++ {
		X.Set(i, 0, 1.0) // intercept
		for j := 0; j < ncov; j++ {
			X.Set(i, j+1, cov.At(i, j))
		}
		y0[i] = pheno.At(i, 0) - center
	}
	y0v := mat.NewVecDense(n, y0)

	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	var Xty mat.VecDense
	Xty.MulVec(X.T(), y0v)

	return localNull{X: X, Y0: y0, XtX: &XtX, Xty0: vecToSlice(&Xty), Y0ty0: mat.Dot(y0v, y0v)}
}

// matrices flattens the c-dim aggregates into the plaintext forms the encoders consume;
// the empty (pid-0) localNull yields c×c / c zeros so AggregateCMat sums align.
func (ln localNull) matrices(c int) (xtx [][]float64, xty []float64, y0ty0 float64) {
	xtx = make([][]float64, c)
	for j := range xtx {
		xtx[j] = make([]float64, c)
	}
	xty = make([]float64, c)
	if ln.XtX == nil {
		return xtx, xty, 0
	}
	for j := 0; j < c; j++ {
		for k := 0; k < c; k++ {
			xtx[j][k] = ln.XtX.At(j, k)
		}
		xty[j] = ln.Xty0[j]
	}
	return xtx, xty, ln.Y0ty0
}

// skatNull holds the secure null-model results from the c-dim aggregates.
type skatNull struct {
	betaHat  crypto.CipherVector // β̂ in slots 0..c-1
	betaRep  crypto.CipherVector // betaRep[ℓ] = Enc(β̂_ℓ) in every slot (PART B CKKS score)
	betaSS   mpc_core.RVec       // β̂ in secret shares (exact SS score)
	xtxSS    mpc_core.RMat       // XᵀX in secret shares (reused by the Burden-variance solve)
	xtxL     mpc_core.RMat       // cached Cholesky factor of xtxSS (factor once, choleskySolve per RHS)
	xtxDinv  mpc_core.RVec       // 1/diag(xtxL); pairs with xtxL for every moment/burden solve
	rssSS    mpc_core.RElem      // residual-norm RSS = y₀ᵀy₀ − 2Xᵀy₀·β̂ + β̂ᵀXᵀXβ̂ (robust σ̂²)
	xtyEnc   crypto.CipherVector // global Xᵀy₀ (aggregated; xtySS derives from it in the null)
	y0ty0Enc crypto.CipherVector // global y₀ᵀy₀ (for RSS)
	c        int
	center   float64 // public centering constant, reused by the score
}

// ridgeRel: tiny Tikhonov ridge on the XᵀX diagonal (hub, once, intercept excluded)
// to keep the inverse finite on a singular design; negligible on a well-conditioned one.
const ridgeRel = 1e-6

// gtgChunkRows bounds how many GᵀG rows enter one secret SSMultMat in the Burden-variance path, so
// a large-m gene's m×m gram is revealed in O(gtgChunkRows·m) pieces instead of one O(m²) blob (OOM).
const gtgChunkRows = 256

// choleskyFactor returns the Cholesky factor L and dInv=1/diag(L) of an SPD matrix A (A=L·Lᵀ) in
// secret shares. The factor is RHS-independent: a reused A (e.g. XᵀX) is factored once and every
// solve against it goes through choleskySolve / choleskySolveMat with the stored (L, dInv).
func (ast *AssocTest) choleskyFactor(A mpc_core.RMat) (mpc_core.RMat, mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	c := len(A)
	mul := func(x, y mpc_core.RElem) mpc_core.RElem { // SS scalar product, back to fracBits
		return mpcObj.TruncVec(mpcObj.SSMultElemVec(mpc_core.RVec{x}, mpc_core.RVec{y}), db, fb)[0]
	}
	L := mpc_core.InitRMat(rtype.Zero(), c, c) // lower-triangular Cholesky factor
	dInv := make(mpc_core.RVec, c)             // dInv[j] = 1/L[j][j]
	for j := 0; j < c; j++ {
		d := A[j][j].Copy() // d = A[j][j] − Σ_{k<j} L[j][k]²
		for k := 0; k < j; k++ {
			d = d.Sub(mul(L[j][k], L[j][k]))
		}
		sq, sqInv := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{d}, false)
		L[j][j], dInv[j] = sq[0], sqInv[0]
		for i := j + 1; i < c; i++ { // L[i][j] = (A[i][j] − Σ_{k<j} L[i][k]L[j][k])/L[j][j]
			v := A[i][j].Copy()
			for k := 0; k < j; k++ {
				v = v.Sub(mul(L[i][k], L[j][k]))
			}
			L[i][j] = mul(v, dInv[j])
		}
	}
	return L, dInv
}

// choleskySolve solves A·x=b given the precomputed factor (L,dInv) from choleskyFactor(A), via
// forward (L·z=b) then back (Lᵀ·x=z) substitution.
func (ast *AssocTest) choleskySolve(L mpc_core.RMat, dInv mpc_core.RVec, b mpc_core.RVec) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	c := len(b)
	mul := func(x, y mpc_core.RElem) mpc_core.RElem {
		return mpcObj.TruncVec(mpcObj.SSMultElemVec(mpc_core.RVec{x}, mpc_core.RVec{y}), db, fb)[0]
	}
	z := make(mpc_core.RVec, c) // forward:  L·z = b
	for i := 0; i < c; i++ {
		v := b[i].Copy()
		for k := 0; k < i; k++ {
			v = v.Sub(mul(L[i][k], z[k]))
		}
		z[i] = mul(v, dInv[i])
	}
	x := make(mpc_core.RVec, c) // back:  Lᵀ·x = z
	for i := c - 1; i >= 0; i-- {
		v := z[i].Copy()
		for k := i + 1; k < c; k++ {
			v = v.Sub(mul(L[k][i], x[k]))
		}
		x[i] = mul(v, dInv[i])
	}
	return x
}

// choleskySolveMat solves A·X=B for a c×R RHS matrix B given the factor (L,dInv), with the
// forward/back substitution running over all R columns together (one length-R round per (i,k) step).
func (ast *AssocTest) choleskySolveMat(L mpc_core.RMat, dInv mpc_core.RVec, B mpc_core.RMat) mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	c := len(dInv)
	R := 0
	if len(B) > 0 {
		R = len(B[0])
	}
	smul := func(s mpc_core.RElem, row mpc_core.RVec) mpc_core.RVec { // secret scalar × length-R row
		return mpcObj.TruncVec(mpcObj.SSMultElemVec(mpc_core.InitRVec(s, R), row), db, fb)
	}
	Z := mpc_core.InitRMat(rtype.Zero(), c, R) // forward:  L·Z = B
	for i := 0; i < c; i++ {
		v := make(mpc_core.RVec, R)
		for p := 0; p < R; p++ {
			v[p] = B[i][p].Copy()
		}
		for k := 0; k < i; k++ {
			t := smul(L[i][k], Z[k])
			for p := 0; p < R; p++ {
				v[p] = v[p].Sub(t[p])
			}
		}
		Z[i] = smul(dInv[i], v)
	}
	X := mpc_core.InitRMat(rtype.Zero(), c, R) // back:  Lᵀ·X = Z
	for i := c - 1; i >= 0; i-- {
		v := make(mpc_core.RVec, R)
		for p := 0; p < R; p++ {
			v[p] = Z[i][p].Copy()
		}
		for k := i + 1; k < c; k++ {
			t := smul(L[k][i], X[k])
			for p := 0; p < R; p++ {
				v[p] = v[p].Sub(t[p])
			}
		}
		X[i] = smul(dInv[i], v)
	}
	return X
}

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

// computeBetaHatEnc computes β̂ = (XᵀX)⁻¹Xᵀy₀ as an encrypted c-vector from the cross-party
// normal equations — only the c-dim aggregates cross the secure boundary, n never does.
func (ast *AssocTest) computeBetaHatEnc() skatNull {
	cps := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()

	c := ast.general.gwasParams.NumCov() + 1
	if c > cps.GetSlots() {
		panic(fmt.Sprintf("computeBetaHatEnc: c=%d exceeds slots=%d", c, cps.GetSlots()))
	}

	center := 0.0
	if ast.general.config.BinaryPheno {
		center = 0.5 // binary {0,1} phenotype; absorbed by the intercept (Q-invariant)
	}
	ln := localNullEquations(ast.general.cov, ast.general.pheno, center)
	xtxLocal, xtyLocal, y0ty0Local := ln.matrices(c)

	// Ridge on the hub only (AggregateCMat sums per-party XᵀX, so add ε once); skip intercept.
	if pid == mpcObj.GetHubPid() {
		var trace float64
		for k := 1; k < c; k++ {
			trace += xtxLocal[k][k]
		}
		eps := ridgeRel * (trace / float64(c))
		for k := 1; k < c; k++ {
			xtxLocal[k][k] += eps
		}
	}

	tNull := time.Now()
	xtxEnc, _, _, err := crypto.EncryptFloatMatrixRow(cps, xtxLocal)
	if err != nil {
		panic(err)
	}
	xtyEnc, _ := crypto.EncryptFloatVector(cps, xtyLocal)
	y0ty0Enc, _ := crypto.EncryptFloatVector(cps, []float64{y0ty0Local})

	xtxEnc = mpcObj.Network.AggregateCMat(cps, xtxEnc)
	if rows := mpcObj.Network.AggregateCMat(cps, crypto.CipherMatrix{xtyEnc}); len(rows) > 0 {
		xtyEnc = rows[0]
	} else {
		xtyEnc = nil
	}
	if rows := mpcObj.Network.AggregateCMat(cps, crypto.CipherMatrix{y0ty0Enc}); len(rows) > 0 {
		y0ty0Enc = rows[0]
	} else {
		y0ty0Enc = nil
	}
	xtxEnc = mpcObj.Network.CollectiveBootstrapMat(cps, xtxEnc, -1)
	fedTimings.nullAgg = time.Since(tNull)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: encrypt+AggregateCMat(XtX,Xty,y0ty0)+bootstrap %v", fedTimings.nullAgg.Round(time.Millisecond)))
	tNull = time.Now()

	// Solve (XᵀX)·β̂ = Xᵀy₀ for the small SPD system by Cholesky + forward/back substitution
	xtxSS := mpcObj.CMatToSS(cps, rtype, xtxEnc, -1, c, 1, c)
	var xtyCt *rlwe.Ciphertext
	if len(xtyEnc) > 0 {
		xtyCt = xtyEnc[0]
	}
	xtySS := mpcObj.CiphertextToSS(cps, rtype, xtyCt, mpcObj.GetHubPid(), c)
	xtxL, xtxDinv := ast.choleskyFactor(xtxSS) // factor once; reused by every moment/burden solve
	betaSS := ast.choleskySolve(xtxL, xtxDinv, xtySS)
	fedTimings.nullInv = time.Since(tNull)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: XtX->SS + Cholesky solve %v", fedTimings.nullInv.Round(time.Millisecond)))
	tNull = time.Now()

	// Residual-norm RSS (robust σ̂²): RSS = y₀ᵀy₀ − 2·(Xᵀy₀·β̂) + β̂ᵀ(XᵀX)β̂, 2nd-order in the β̂ error
	// (the plain identity y₀ᵀy₀ − Xᵀy₀·β̂ is 1st-order). All from the c-dim SS aggregates.
	sumRVec := func(v mpc_core.RVec) mpc_core.RElem {
		acc := rtype.Zero()
		for k := 0; k < len(v); k++ {
			acc = acc.Add(v[k])
		}
		return acc
	}
	betaCol2 := make(mpc_core.RMat, c)
	for k := 0; k < c; k++ {
		betaCol2[k] = mpc_core.RVec{betaSS[k]}
	}
	xtxBetaMat := mpcObj.SSMultMat(xtxSS, betaCol2) // (XᵀX)·β̂
	xtxBeta := make(mpc_core.RVec, c)
	for k := 0; k < c; k++ {
		xtxBeta[k] = xtxBetaMat[k][0]
	}
	xtxBeta = mpcObj.TruncVec(xtxBeta, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	xtyBeta := sumRVec(mpcObj.TruncVec(mpcObj.SSMultElemVec(xtySS, betaSS), mpcObj.GetDataBits(), mpcObj.GetFracBits()))
	betaXtxBeta := sumRVec(mpcObj.TruncVec(mpcObj.SSMultElemVec(betaSS, xtxBeta), mpcObj.GetDataBits(), mpcObj.GetFracBits()))
	var y0Ct *rlwe.Ciphertext
	if len(y0ty0Enc) > 0 {
		y0Ct = y0ty0Enc[0]
	}
	y0ty0SS := mpcObj.CiphertextToSS(cps, rtype, y0Ct, mpcObj.GetHubPid(), 1)[0]
	rssSS := y0ty0SS.Sub(xtyBeta).Sub(xtyBeta).Add(betaXtxBeta) // y₀ᵀy₀ − 2·Xᵀy₀·β̂ + β̂ᵀXᵀXβ̂

	betaHat := mpcObj.SSToCVec(cps, betaSS)

	// betaRep[ℓ] = β̂_ℓ replicated in every slot (for the score's CPMult). pid-0 → nil.
	slots := cps.GetSlots()
	betaRep := make(crypto.CipherVector, c)
	for j := 0; j < c; j++ {
		rep := mpc_core.InitRVec(rtype.Zero(), slots)
		if j < len(betaSS) {
			for k := 0; k < slots; k++ {
				rep[k] = betaSS[j].Copy()
			}
		}
		if cv := mpcObj.SSToCVec(cps, rep); len(cv) > 0 {
			betaRep[j] = cv[0]
		}
	}

	fedTimings.nullBeta = time.Since(tNull)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: Xty->SS + beta=inv*Xty + betaRep %v", fedTimings.nullBeta.Round(time.Millisecond)))
	return skatNull{betaHat: betaHat, betaRep: betaRep, betaSS: betaSS, xtxSS: xtxSS, xtxL: xtxL, xtxDinv: xtxDinv, rssSS: rssSS, xtyEnc: xtyEnc, y0ty0Enc: y0ty0Enc, c: c, center: center}
}

// computeNullRSSEnc wraps the residual-norm RSS (formed in SS as null.rssSS) as a 1-element
// CipherVector for rareVariantScaleShares. SSToCVec is collective, so all parties call it.
func (ast *AssocTest) computeNullRSSEnc(null skatNull) crypto.CipherVector {
	if null.rssSS == nil {
		return nil
	}
	return ast.general.mpcObj[0].SSToCVec(ast.general.cps, mpc_core.InitRVec(null.rssSS, 1))
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

// skatBlockNumSnps returns block b's SNP count, syncing pid 0 over the network.
func (ast *AssocTest) skatBlockNumSnps(block int) int {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	hubPid := mpcObj.GetHubPid()
	if pid == 0 {
		return mpcObj.Network.ReceiveInt(hubPid)
	}

	var nsnpsBlock int
	blockSize := ast.general.genoBlockSizes[block]
	shift := uint64(0)
	for i := 0; i < block; i++ {
		shift += uint64(ast.general.genoBlockSizes[i])
	}
	if ast.general.IsPgen() {
		if ast.general.gwasParams.snpFilt == nil {
			nsnpsBlock = blockSize
		} else {
			nsnpsBlock = SumBool(ast.general.gwasParams.snpFilt[shift : shift+uint64(blockSize)])
		}
	} else {
		nsnpsBlock = int(ast.general.genoBlocks[block].NumColsToKeep())
	}

	if pid == hubPid {
		mpcObj.Network.SendInt(nsnpsBlock, 0)
	}
	return nsnpsBlock
}

// openBlockGenoStream returns this party's genotype stream for block b (the open "blocks"
// stream, or a pgen block materialized to a temp file). nil for an empty block.
func (ast *AssocTest) openBlockGenoStream(b int) *GenoFileStream {
	if !ast.general.IsPgen() {
		if b < len(ast.general.genoBlocks) {
			return ast.general.genoBlocks[b]
		}
		return nil
	}

	gp := ast.general.gwasParams
	blockSize := ast.general.genoBlockSizes[b]
	shift := 0
	for i := 0; i < b; i++ {
		shift += ast.general.genoBlockSizes[i]
	}
	var snpFilt []bool
	if gp.snpFilt == nil {
		snpFilt = OnesBool(blockSize)
	} else {
		snpFilt = gp.snpFilt[shift : shift+blockSize]
	}
	nsnps := SumBool(snpFilt)
	if nsnps == 0 {
		return nil
	}
	numInd := ast.skatNumInds()[ast.general.mpcObj[0].GetPid()]
	pgenFile := fmt.Sprintf(ast.general.config.GenoFilePrefix, b+1)
	tmp := ast.general.CachePath(fmt.Sprintf("lowrank_pgen_gfs.%d.tmp", b))
	FilterMatrixFilePgen(pgenFile, numInd, nsnps, ast.general.config.SampleKeepFile,
		ast.general.config.SnpIdsFile, shift, snpFilt, tmp, false)
	return NewGenoFileStream(tmp, uint64(numInd), uint64(nsnps), true)
}

// denseFromStream reads a genotype stream into a dense matrix (samples × variants), missing
// (negative) dosages → 0. nil/empty stream yields a 0×0 matrix.
func denseFromStream(gfs *GenoFileStream) *mat.Dense {
	if gfs == nil {
		return mat.NewDense(0, 0, nil)
	}
	gfs.Reset()

	var rows [][]float64
	for {
		row := gfs.NextRow()
		if row == nil {
			break
		}
		fr := make([]float64, len(row))
		for j, v := range row {
			if v > 0 {
				fr[j] = float64(v)
			}
		}
		rows = append(rows, fr)
	}
	if len(rows) == 0 || len(rows[0]) == 0 {
		return nil // empty block (e.g. gene with no private variants); callers treat nil as empty
	}
	G := mat.NewDense(len(rows), len(rows[0]), nil)
	for i := range rows {
		G.SetRow(i, rows[i])
	}
	return G
}

// loadDenseBlocks reads nGenes per-gene int8 genotype block files into dense matrices for the
// federated-private private side. path b = fmt.Sprintf(prefix, b), row-major n×m_b with m_b
// inferred from file size / n (int8 = 1 byte/cell).
func loadDenseBlocks(prefix string, nGenes, n int) ([]*mat.Dense, error) {
	if n <= 0 {
		return nil, fmt.Errorf("loadDenseBlocks: n=%d must be positive", n)
	}
	blocks := make([]*mat.Dense, nGenes)
	for b := 0; b < nGenes; b++ {
		path := fmt.Sprintf(prefix, b)
		fi, err := os.Stat(path)
		if err != nil {
			return nil, err
		}
		if fi.Size()%int64(n) != 0 {
			return nil, fmt.Errorf("loadDenseBlocks: %s size %d not divisible by n=%d", path, fi.Size(), n)
		}
		m := int(fi.Size()) / n
		if m == 0 {
			blocks[b] = nil // gene with no private variants; privateQRaw(nil) -> Enc(0)
			continue
		}
		blocks[b] = denseFromStream(NewGenoFileStream(path, uint64(n), uint64(m), false))
	}
	return blocks, nil
}

// readGenoBlockLocal streams this party's genotype block b into a dense matrix
// (samples × variants), missing (negative) dosages → 0. Empty on pid 0 or an empty block.
func (ast *AssocTest) readGenoBlockLocal(b int) *mat.Dense {
	return denseFromStream(ast.openBlockGenoStream(b))
}

// nullSetup runs the secure null model (β̂ + RSS) and builds the block-independent
// per-party design X=[1|cov] and centered y₀ (empty on pid 0).
func (ast *AssocTest) nullSetup() (null skatNull, nullRSS crypto.CipherVector, X *mat.Dense, y0 []float64) {
	null = ast.computeBetaHatEnc()
	tRSS := time.Now()
	nullRSS = ast.computeNullRSSEnc(null)
	fedTimings.nullRSS = time.Since(tRSS)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: RSS = y0ty0 - Xty·beta %v", fedTimings.nullRSS.Round(time.Millisecond)))
	if ast.general.mpcObj[0].GetPid() > 0 {
		center := 0.0
		if ast.general.config.BinaryPheno {
			center = 0.5
		}
		ln := localNullEquations(ast.general.cov, ast.general.pheno, center)
		X, y0 = ln.X, ln.Y0
	}
	return
}

// blockStat computes one block's raw statistics as 1-elem RVecs: qRawSS = Σ ŵ²s²
// and bLinSS = Σ ŵ·s (burden linear term, squared by the caller). Local plaintext
// contraction → key-free score → AggregateCMat → oriented weights → ScoreCalculation → SS.
func (ast *AssocTest) blockStat(b, nsnps int, null skatNull, X *mat.Dense, y0 []float64, gl *geneLocal) (qRawSS, bLinSS, wSignedSS mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	pid := mpcObj.GetPid()

	var GtX *mat.Dense
	var Gty0 []float64
	dosage := make([]float64, nsnps)
	if pid > 0 {
		var lc LocalContraction
		if gl != nil && gl.Gloc != nil {
			lc = gl.LocalContraction
		} else {
			lc = LocalContract(ast.readGenoBlockLocal(b), X, y0)
		}
		GtX, Gty0, dosage = lc.GtX, lc.Gty0, lc.DosageSum
	}

	// Score in secret shares — the Gᵀy₀ − GᵀX·β̂ cancellation is exact in fixed-point (β̂ from the
	// Cholesky null solve is accurate). All parties (incl. pid 0, with zero shares) take part.
	sCVec := mpcObj.SSToCVec(cps, ast.scoreSS(GtX, Gty0, null.betaSS, nsnps))
	wSigned := ast.signedWeight(dosage, nsnps) // SS signed weight; also returned for the Burden-variance solve
	weightEnc := mpcObj.SSToCVec(cps, wSigned)

	var qRaw, bLin crypto.CipherVector
	if pid > 0 && len(sCVec) > 0 {
		qRaw, bLin, _, _, _, _ = ast.ScoreCalculation(sCVec, weightEnc)
	}
	return ast.scalarCiphertextToShares(qRaw), ast.scalarCiphertextToShares(bLin), wSigned
}

// burdenVarSS returns one gene's Burden-variance factor zᵀPz = z_fullᵀP z_full (P = I − X(XᵀX)⁻¹Xᵀ,
// z_full = Σⱼ ŵⱼ Gⱼ over the public list ∪ private variants) as an SS scalar. It never forms z (which
// is n-dim); instead it works in the small m/c space:
//
//	zᵀPz = ŵ_pubᵀ(GᵀG)_pub ŵ_pub  +  2·ŵ_pubᵀ d  +  z_privᵀz_priv  −  (Xᵀz)ᵀ(XtX)⁻¹(Xᵀz)
//
// where the public GᵀG/GᵀX are federated by the "local contraction = SS share → global sum" trick
// (same as scoreSS), and the private party contributes d = G_pubᵀz_priv (m_pub), z_privᵀz_priv, and
// Xᵀz_priv (all n-contracted locally, so the private variant count m_priv stays hidden).
func (ast *AssocTest) burdenVarSS(b, nsnps int, null skatNull, X *mat.Dense, y0 []float64, privG *mat.Dense, privatePid int, wPubIn mpc_core.RVec, gl *geneLocal) mpc_core.RElem {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	pid := mpcObj.GetPid()
	c := null.c

	// zᵀPz/N (else ŵᵀGᵀGŵ ~2^31 at AoU overflows the fixed-point wall): quadratic pieces /N, Xᵀz pieces
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
	var Gloc, gg *mat.Dense // Gloc = local aligned public genotype (reused by the cross term); gg = local GᵀG
	gtxSS := mpc_core.InitRMat(rtype.Zero(), nsnps, c)
	if pid > 0 && nsnps > 0 {
		g := ast.localFor(b, nsnps, X, y0, gl)
		Gloc, gg, dosage = g.Gloc, g.gg, g.DosageSum
		for j := 0; j < nsnps; j++ {
			for l := 0; l < c; l++ {
				gtxSS[j][l] = rtype.FromFloat64(g.GtX.At(j, l)/sqrtN, fb) // GᵀX/√N
			}
		}
	}
	wPub := wPubIn // signed weight reused from PART A's blockStat (same gene/dosage → same value)
	if wPub == nil {
		wPub = ast.signedWeight(dosage, nsnps) // fallback: compute here when not supplied
	}

	// pubZZ = ŵ_pubᵀ(GᵀG)ŵ_pub ; pubXtz = (XᵀG_pub)ŵ_pub (c-vector)
	pubZZ := rtype.Zero()
	pubXtz := mpc_core.InitRVec(rtype.Zero(), c)
	if nsnps > 0 {
		wCol := make(mpc_core.RMat, nsnps)
		for j := 0; j < nsnps; j++ {
			wCol[j] = mpc_core.RVec{wPub[j]}
		}
		// (GᵀG)·ŵ in row-chunks: forming the whole m×m secret gram (and its single Beaver reveal)
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
		// pubXtz = (GᵀX)ᵀ·ŵ_pub (c-vector) as one SSMultMat over the c×nsnps transpose.
		gtxT := mpc_core.InitRMat(rtype.Zero(), c, nsnps)
		for l := 0; l < c; l++ {
			for j := 0; j < nsnps; j++ {
				gtxT[l][j] = gtxSS[j][l]
			}
		}
		pxM := mpcObj.SSMultMat(gtxT, wCol) // c×1 (untruncated 2·fb)
		for l := 0; l < c; l++ {
			pubXtz[l] = pxM[l][0]
		}
		pubXtz = mpcObj.TruncVec(pubXtz, db, fb)
	}

	// --- private (privatePid only): z_priv = G_priv·ŵ_priv, contracted locally (count hidden) ---
	privZZ := rtype.Zero()
	crossD := mpc_core.InitRVec(rtype.Zero(), nsnps) // d = G_pubᵀ z_priv (m_pub)
	privXtz := mpc_core.InitRVec(rtype.Zero(), c)    // Xᵀ z_priv (c)
	if pid == privatePid && privG != nil {
		np, mp := privG.Dims()
		if mp > 0 {
			dpriv := make([]float64, mp)
			for k := 0; k < mp; k++ {
				for i := 0; i < np; i++ {
					dpriv[k] += privG.At(i, k)
				}
			}
			wpv := skatBetaWeightSigned(dpriv, ast.skatTotalNumInds())
			zpv := make([]float64, np)
			for i := 0; i < np; i++ {
				for k := 0; k < mp; k++ {
					zpv[i] += privG.At(i, k) * wpv[k]
				}
			}
			zz := 0.0
			for i := 0; i < np; i++ {
				zz += zpv[i] * zpv[i]
			}
			privZZ = rtype.FromFloat64(zz/N, fb) // priv-priv quadratic /N
			for l := 0; l < c; l++ {
				s := 0.0
				for i := 0; i < np; i++ {
					s += X.At(i, l) * zpv[i]
				}
				privXtz[l] = rtype.FromFloat64(s/sqrtN, fb) // Xᵀz_priv /√N
			}
			if nsnps > 0 { // cross term needs B's aligned public genotype (Gloc)
				for j := 0; j < nsnps; j++ {
					s := 0.0
					for i := 0; i < np; i++ {
						s += Gloc.At(i, j) * zpv[i]
					}
					crossD[j] = rtype.FromFloat64(s/N, fb) // G_pubᵀz_priv /N (cross quadratic)
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

	// zᵀPz = zz − (Xᵀz)ᵀ(XtX)⁻¹(Xᵀz), reusing the null Cholesky solve.
	a := ast.choleskySolve(null.xtxL, null.xtxDinv, xtz)
	return zz.Sub(vdot(xtz, a))
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

// ComputeSKATStatistics returns the whole-genome secure SKAT (Q) and Burden as
// 1-elem CipherVectors (Σ over all blocks, scaled by 1/(2σ̂²)). Caller collectively decrypts.
func (ast *AssocTest) ComputeSKATStatistics() (qStat, qBurden crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()

	null, nullRSS, X, y0 := ast.nullSetup()

	finalQSS := mpc_core.InitRVec(rtype.Zero(), 1)
	finalBurdenSS := mpc_core.InitRVec(rtype.Zero(), 1)
	for b := 0; b < ast.general.config.GenoNumBlocks; b++ {
		nsnps := ast.skatBlockNumSnps(b)
		if nsnps == 0 {
			continue
		}
		qRawSS, bLinSS, _ := ast.blockStat(b, nsnps, null, X, y0, nil)
		finalQSS.Add(qRawSS)
		finalBurdenSS.Add(bLinSS)
	}

	finalBurdenSS = mpcObj.SSSquareElemVec(finalBurdenSS)
	finalBurdenSS = mpcObj.TruncVec(finalBurdenSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	if scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS); ok {
		finalQSS = ast.general.scaleRareVariantShareStat(finalQSS, scaleSS)
		finalBurdenSS = ast.general.scaleRareVariantShareStat(finalBurdenSS, scaleSS)
	}
	return mpcObj.SSToCVec(cps, finalQSS), mpcObj.SSToCVec(cps, finalBurdenSS)
}

// ComputeSKATStatisticsPerBlock returns per-block Q and Burden (slot b = block b's
// statistic, scaled by the common 1/(2σ̂²)) — the per-gene secure SKAT statistics.
func (ast *AssocTest) ComputeSKATStatisticsPerBlock() (qPerBlock, burdenPerBlock crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()
	hub := mpcObj.GetPid() == mpcObj.GetHubPid()

	t0 := time.Now()
	null, nullRSS, X, y0 := ast.nullSetup()
	if hub {
		log.LLvl1(fmt.Sprintf(">>> [secure] null model (β̂, RSS) done (%.1fs)", time.Since(t0).Seconds()))
	}

	nB := ast.general.config.GenoNumBlocks
	qBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)
	bLinBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)
	done := 0
	for b := 0; b < nB; b++ {
		nsnps := ast.skatBlockNumSnps(b)
		if nsnps == 0 {
			continue
		}
		tb := time.Now()
		q, bl, _ := ast.blockStat(b, nsnps, null, X, y0, nil)
		qBlockSS[b] = q[0]
		bLinBlockSS[b] = bl[0]
		done++
		if hub {
			log.LLvl1(fmt.Sprintf(">>> [secure] block %d/%d done (%d snps, %.1fs; elapsed %.1fs)",
				done, nB, nsnps, time.Since(tb).Seconds(), time.Since(t0).Seconds()))
		}
	}

	burdenSqSS := mpcObj.SSSquareElemVec(bLinBlockSS)
	burdenSqSS = mpcObj.TruncVec(burdenSqSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	if scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS); ok {
		scaleVec := mpc_core.InitRVec(rtype.Zero(), nB)
		for b := 0; b < nB; b++ {
			scaleVec[b] = scaleSS[0].Copy()
		}
		qBlockSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(qBlockSS, scaleVec), mpcObj.GetDataBits(), mpcObj.GetFracBits())
		burdenSqSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(burdenSqSS, scaleVec), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	}
	return mpcObj.SSToCVec(cps, qBlockSS), mpcObj.SSToCVec(cps, burdenSqSS)
}

func (ast *AssocTest) skatNumInds() []int {
	filtNumInds := ast.general.gwasParams.FiltNumInds()
	if len(filtNumInds) == ast.general.config.NumMainParties+1 {
		return filtNumInds
	}

	return ast.general.config.NumInds
}

func (ast *AssocTest) skatTotalNumInds() int {
	nrows := ast.skatNumInds()
	total := 0
	for p := 1; p <= ast.general.config.NumMainParties; p++ {
		total += nrows[p]
	}
	return total
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

// --- federated SKAT with party-private variants ---
//
// Two parties hold partially-overlapping variants; one party's list is PUBLIC. Per gene a variant
// is shared (both), public_only (public-list party only — other contributes 0) or private (private
// party only, never shared). MVP+AoU is the motivating instance (MVP=public list, AoU=private).
//
//	per-gene Q = PART A (secure over the public list) + PART B (private, computed locally)
//
// Both parts are raw Σŵ²s²; the common 1/(2σ̂²) is applied once at the end. Equals the pooled
// single-cohort SKAT; plaintext oracle = SKATFederatedPrivate (skat_plain.go).

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

// skatBetaWeightSigned folds the minor-allele orientation (ŵ_j = −w_j iff p̄_j>½) into the weight, so
// ScoreCalculation yields both Q (sign²=1, same as unsigned) and the R::SKAT-oriented Burden linear
// term Σŵ·s. The private party's dosage is local, so its sign is computed in plaintext.
func skatBetaWeightSigned(dosageSum []float64, totalInds int) []float64 {
	w := skatBetaWeight(dosageSum, totalInds)
	twoN := float64(2 * totalInds)
	for j, d := range dosageSum {
		if d/twoN > 0.5 {
			w[j] = -w[j]
		}
	}
	return w
}

// privateRawStats returns the private party's local raw SKAT = Σ ŵ² s² and Burden linear term Σ ŵ·s
// (two ciphertext scalars in slot 0) for one gene's private variants: key-free score (β̂ stays
// encrypted) × plaintext signed weight, no AggregateCMat. The signed weight leaves SKAT unchanged
// (sign²=1) and makes the linear term match R::SKAT's orientation. Enc(0) for an empty gene. Only on
// the private party.
func (ast *AssocTest) privateRawStats(G *mat.Dense, null skatNull, X *mat.Dense, y0 []float64) (skat, burdenLin crypto.CipherVector) {
	cps := ast.general.cps

	var m int
	if G != nil {
		_, m = G.Dims()
	}
	if m == 0 {
		z, _ := crypto.EncryptFloatVector(cps, []float64{0})
		zb, _ := crypto.EncryptFloatVector(cps, []float64{0})
		return z, zb
	}

	lc := LocalContract(G, X, y0)
	s := ast.partyScore(lc.GtX, lc.Gty0, null)
	wEnc, _ := crypto.EncryptFloatVector(cps, skatBetaWeightSigned(lc.DosageSum, ast.skatTotalNumInds()))

	s, wEnc = alignCipherVectorLevels(cps, s, wEnc)
	skat, burdenLin, _, _, _, _ = ast.ScoreCalculation(s, wEnc)
	return skat, burdenLin
}

// privateBlockStat returns one gene's private-variant raw SKAT and Burden linear term as 1-elem SS
// RVecs. The private party computes them locally (privateRawStats); CiphertextToSS reshares each from
// privatePid (the other party gets only the ciphertext, pid 0 sits out). Must run for EVERY gene —
// the Enc(0) on empty genes keeps which genes have private variants indistinguishable.
func (ast *AssocTest) privateBlockStat(G *mat.Dense, null skatNull, X *mat.Dense, y0 []float64, privatePid int) (skatSS, burdenSS mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()

	var skatCt, burdenCt *rlwe.Ciphertext
	if mpcObj.GetPid() == privatePid {
		skat, burdenLin := ast.privateRawStats(G, null, X, y0)
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
//
// fedTimings accumulates skat_fed phase durations for the end-of-run tree (runFederatedPrivate
// prints it, combined with mpc.SetupTiming). One run per process; not concurrency-safe.
var fedTimings struct {
	initTotal, loadPriv, assocInit                                time.Duration // pre-compute phases
	nullAgg, nullInv, nullBeta, nullRSS, nullTotal, blocks, total time.Duration
	blockSecs                                                     []float64
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
	null, nullRSS, X, y0 := ast.nullSetup()
	fedTimings.nullTotal = time.Since(tStart)
	log.LLvl1(fmt.Sprintf("[skat_fed] null model: %v", fedTimings.nullTotal.Round(time.Millisecond)))

	nB := ast.general.config.GenoNumBlocks
	if mpcObj.GetPid() == privatePid && privateOnly != nil && len(privateOnly) != nB {
		panic(fmt.Sprintf("ComputeSKATFederatedPrivate: privateOnly has %d blocks, want %d", len(privateOnly), nB))
	}
	skatBlockSS := mpc_core.InitRVec(rtype.Zero(), nB) // Σ ŵ²s²   (SKAT raw, per gene)
	bLinBlockSS := mpc_core.InitRVec(rtype.Zero(), nB) // Σ ŵ·s    (Burden linear term, per gene)
	zpzBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)  // zᵀPz     (Burden variance, per gene; unscaled)
	nProbes := ast.general.config.SkatPValueProbes     // SKAT p-value (Hutchinson); 0 = disabled
	s1B := mpc_core.InitRVec(rtype.Zero(), nB)         // SKAT kernel moments per gene (N-normalized), if enabled
	s2B := mpc_core.InitRVec(rtype.Zero(), nB)
	s3B := mpc_core.InitRVec(rtype.Zero(), nB)
	tBlocks := time.Now()
	blockSecs := make([]float64, 0, nB)
	// ETA weight: the moment cost is O(m_pub²·probes), so weight each gene by nsnps² of its public block.
	var totalWork, doneWork float64
	for b := 0; b < nB; b++ {
		m := float64(ast.general.genoBlockSizes[b])
		totalWork += m * m
	}
	for b := 0; b < nB; b++ {
		tb := time.Now()
		accSkat := mpc_core.InitRVec(rtype.Zero(), 1)
		accBurden := mpc_core.InitRVec(rtype.Zero(), 1)
		nsnps := ast.skatBlockNumSnps(b) // collective; public-list size for gene b
		// Not hub-gated: run_fed.sh tees party 2 (not the hub) to the terminal, so every party logs this.
		if nProbes > 0 { // moments: tr(K_ppᵏ) Hutchinson (O(m_pub²·probes)) + block-contracted private
			log.LLvl1(fmt.Sprintf("[skat_fed] gene %d/%d start (m_pub=%d, probes=%d)", b+1, nB, nsnps, nProbes))
		} else {
			log.LLvl1(fmt.Sprintf("[skat_fed] gene %d/%d start (m_pub=%d)", b+1, nB, nsnps))
		}

		// Read + contract + Gram the public block ONCE per gene; blockStat/burdenVarSS/skatMomentsSS reuse it.
		gl := ast.computeGeneLocal(b, nsnps, X, y0)

		// PART A: secure SKAT over the public list (existing per-block path).
		var wSignedA mpc_core.RVec // signed weight from PART A, reused by burdenVarSS below
		var tBlk, tBur, tMom float64
		if nsnps > 0 {
			t := time.Now()
			skatA, burdenA, wA := ast.blockStat(b, nsnps, null, X, y0, gl)
			tBlk = time.Since(t).Seconds()
			accSkat.Add(skatA)
			accBurden.Add(burdenA)
			wSignedA = wA
		}

		// PART B: private variants for this gene (uniform across all genes).
		var G *mat.Dense
		if mpcObj.GetPid() == privatePid && b < len(privateOnly) {
			G = privateOnly[b]
		}
		skatB, burdenB := ast.privateBlockStat(G, null, X, y0, privatePid)
		accSkat.Add(skatB)
		accBurden.Add(burdenB)

		skatBlockSS[b] = accSkat[0]
		bLinBlockSS[b] = accBurden[0]
		tBurStart := time.Now()
		zpzBlockSS[b] = ast.burdenVarSS(b, nsnps, null, X, y0, G, privatePid, wSignedA, gl) // Burden variance zᵀPz
		tBur = time.Since(tBurStart).Seconds()
		if nProbes > 0 { // SKAT p-value kernel moments (N-normalized)
			t := time.Now()
			s1B[b], s2B[b], s3B[b] = ast.skatMomentsSS(b, nsnps, nProbes, null, X, y0, G, gl)
			tMom = time.Since(t).Seconds()
		}
		dt := time.Since(tb).Seconds()
		blockSecs = append(blockSecs, dt)
		// ETA weights remaining genes by m² (the O(m_pub²·R) moment cost dominates).
		doneWork += float64(nsnps) * float64(nsnps)
		elapsed := time.Since(tBlocks).Seconds()
		eta := 0.0
		if doneWork > 0 {
			eta = elapsed * (totalWork - doneWork) / doneWork
		}
		log.LLvl1(fmt.Sprintf("[skat_fed] gene %d/%d done %.0fs [blockStat %.0f | burden %.0f | moments %.0f] | elapsed %.0fs  ETA ~%.0fm (%.0f%%)",
			b+1, nB, dt, tBlk, tBur, tMom, elapsed, eta/60, 100*doneWork/totalWork))
	}
	fedTimings.blocks = time.Since(tBlocks)
	fedTimings.blockSecs = blockSecs
	log.LLvl1(fmt.Sprintf("[skat_fed] %d blocks: %v total | per-block %s",
		nB, fedTimings.blocks.Round(time.Millisecond), blockSecStats(blockSecs)))

	// Burden = (Σ ŵ·s)² — square the accumulated linear term (A+B) per gene before scaling.
	// Q_b/N to match zᵀPz/N: scale the score by 1/√N BEFORE squaring (b² hits the fixed-point wall at the square).
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	invSqrtN := rtype.FromFloat64(1.0/math.Sqrt(float64(ast.skatTotalNumInds())), fb)
	bLinNorm := make(mpc_core.RVec, nB)
	for b := 0; b < nB; b++ {
		bLinNorm[b] = bLinBlockSS[b].Mul(invSqrtN)
	}
	bLinNorm = mpcObj.TruncVec(bLinNorm, db, fb)
	burdenBlockSS := mpcObj.TruncVec(mpcObj.SSSquareElemVec(bLinNorm), db, fb)

	// Common 1/(2σ̂²) applied once to both stats (linear, distributes over A+B).
	scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS)
	if !ok {
		panic("ComputeSKATFederatedPrivate: div undefined (dof = n-c <= 0); need more samples than covariates")
	}
	scaleVec := mpc_core.InitRVec(rtype.Zero(), nB)
	for b := 0; b < nB; b++ {
		scaleVec[b] = scaleSS[0].Copy()
	}
	skatBlockSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(skatBlockSS, scaleVec), db, fb)
	burdenBlockSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(burdenBlockSS, scaleVec), db, fb)

	// Burden p-value: reveal ONLY √(T/2) = √(Burden/zᵀPz) = √Burden·(1/√zᵀPz). Both the Burden
	// statistic and zᵀPz stay secret-shared (never decrypted), so neither leaves and zᵀPz=2·Burden/T
	// is not derivable. √Burden and 1/√zᵀPz come from one SqrtAndSqrtInverse each; driver applies erfc.
	// SqrtAndSqrtInverse assumes strictly-positive inputs; Burden=B²/(2σ̂²)≥0 and zᵀPz≥0 hold with a
	// wide margin for any polymorphic gene (MAF-filtered upstream). A degenerate (monomorphic) gene
	// would give an unreliable — but bounded, non-crashing — p; the driver's √(T/2)>0 guard bounds output.
	sqrtBurden, _ := mpcObj.SqrtAndSqrtInverse(burdenBlockSS, false)
	_, invSqrtZpz := mpcObj.SqrtAndSqrtInverse(zpzBlockSS, false)
	sqrtT2 := mpcObj.TruncVec(mpcObj.SSMultElemVec(sqrtBurden, invSqrtZpz), db, fb)

	// SKAT p-value pivot z per gene: Q_norm = (scaled SKAT stat)/N, then WH z = skatZSS(Q,S1,S2,S3).
	// Reveal only z (Q_skat, moments stay secret); driver applies p = ½erfc(z/√2).
	if nProbes > 0 {
		invN := rtype.FromFloat64(1.0/float64(ast.skatTotalNumInds()), fb)
		qNorm := make(mpc_core.RVec, nB)
		for b := 0; b < nB; b++ {
			qNorm[b] = skatBlockSS[b].Mul(invN)
		}
		qNorm = mpcObj.TruncVec(qNorm, db, fb)
		zB := ast.skatZSSVec(qNorm, s1B, s2B, s3B)
		skatZStat = mpcObj.SSToCVec(cps, zB)
	}

	fedTimings.total = time.Since(tStart)
	log.LLvl1(fmt.Sprintf("[skat_fed] total compute: %v", fedTimings.total.Round(time.Millisecond)))
	// SKAT-p mode reveals ONLY z (the equivalent p) — the raw SKAT statistic Q_skat must not leave.
	// Suppress it (nil ⇒ driver skips decrypt/save); the Burden path (nProbes==0) still returns Q_skat.
	var skatStatOut crypto.CipherVector
	if nProbes == 0 {
		skatStatOut = mpcObj.SSToCVec(cps, skatBlockSS)
	}
	return skatStatOut, mpcObj.SSToCVec(cps, sqrtT2), skatZStat
}
