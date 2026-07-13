package gwas

import (
	"fmt"
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
