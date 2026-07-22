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

// localGenotypeContract computes only the gene-dependent fields. Per-gene SKAT callers use this
// instead of recomputing the null-only XtX, Xty0, and y0ty0 for every public/private block.
func localGenotypeContract(G, X *mat.Dense, y0 []float64) LocalContraction {
	n, _ := X.Dims()
	gn, m := G.Dims()
	if gn != n || len(y0) != n {
		panic(fmt.Sprintf("localGenotypeContract: G/X/y0 rows differ (%d/%d/%d)", gn, n, len(y0)))
	}
	y0v := mat.NewVecDense(n, y0)

	var GtX mat.Dense
	GtX.Mul(G.T(), X)

	var Gty mat.VecDense
	Gty.MulVec(G.T(), y0v)

	dosage := make([]float64, m)
	for j := 0; j < m; j++ {
		dosage[j] = GtX.At(j, 0) // X[:,0] is the intercept, so this is already Σ_i G[i,j]
	}

	return LocalContraction{
		GtX:       &GtX,
		Gty0:      vecToSlice(&Gty),
		DosageSum: dosage,
	}
}

// LocalContract contracts (G, X, y0) locally; y0 centered, X includes intercept.
func LocalContract(G, X *mat.Dense, y0 []float64) LocalContraction {
	lc := localGenotypeContract(G, X, y0)
	n, _ := X.Dims()
	y0v := mat.NewVecDense(n, y0)

	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	var Xty mat.VecDense
	Xty.MulVec(X.T(), y0v)
	lc.XtX = &XtX
	lc.Xty0 = vecToSlice(&Xty)
	lc.Y0ty0 = mat.Dot(y0v, y0v)
	return lc
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
	Gloc := orientGenotypeLocal(ast.readGenoBlockLocal(b))
	var gg mat.Dense
	gg.Mul(Gloc.T(), Gloc)
	return &geneLocal{LocalContraction: localGenotypeContract(Gloc, X, y0), Gloc: Gloc, gg: &gg}
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
		if cov != nil || pheno != nil {
			panic("localNullEquations: covariates and phenotype must both be present")
		}
		return localNull{}
	}
	cn, ncov := cov.Dims()
	pn, pc := pheno.Dims()
	if cn != pn || pc == 0 {
		panic(fmt.Sprintf("localNullEquations: covariate/phenotype dimensions differ (%dx%d/%dx%d)", cn, ncov, pn, pc))
	}
	n := cn
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
	betaRep crypto.CipherVector // betaRep[ℓ] = Enc(β̂_ℓ) in every slot (PART B CKKS score)
	betaSS  mpc_core.RVec       // β̂ in secret shares (exact SS score)
	omp     mpc_core.RMat       // Ω'=N(XᵀX)⁻¹, solved once and reused by every gene/action
	rssSS   mpc_core.RElem      // residual-norm RSS = y₀ᵀy₀ − 2Xᵀy₀·β̂ + β̂ᵀXᵀXβ̂ (robust σ̂²)
	c       int
}

// ridgeRel: relative Tikhonov ridge on the covariate diagonal, keeps XᵀX⁻¹ finite on a singular design.
const ridgeRel = 1e-6

// choleskyFactor returns the Cholesky factor L and dInv=1/diag(L) of an SPD matrix A (A=L·Lᵀ) in
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

// nullSetup solves the cross-party c-dimensional normal equations and returns the local X/y₀
// alongside the secure result so per-gene work can reuse them without rebuilding XᵀX.
func (ast *AssocTest) nullSetup() (null skatNull, X *mat.Dense, y0 []float64) {
	cps := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()

	c := ast.general.gwasParams.NumCov() + 1
	if c > cps.GetSlots() {
		panic(fmt.Sprintf("nullSetup: c=%d exceeds slots=%d", c, cps.GetSlots()))
	}

	center := 0.0
	if ast.general.config.BinaryPheno {
		center = 0.5 // binary {0,1} phenotype; absorbed by the intercept (Q-invariant)
	}
	localMark := ast.metricMark()
	ln := localNullEquations(ast.general.cov, ast.general.pheno, center)
	xtxLocal, xtyLocal, y0ty0Local := ln.matrices(c)

	if pid > 0 {
		var trace float64
		for k := 1; k < c; k++ {
			trace += xtxLocal[k][k]
		}
		eps := ridgeRel * (trace / float64(c))
		for k := 1; k < c; k++ {
			xtxLocal[k][k] += eps
		}
	}
	ast.metricEnd("null_local_xtx_xty_yty", localMark)

	xtxMark := ast.metricMark()
	xtxEnc, _, _, err := crypto.EncryptFloatMatrixRow(cps, xtxLocal)
	if err != nil {
		panic(err)
	}
	aggVec := func(v crypto.CipherVector) crypto.CipherVector {
		if rows := mpcObj.Network.AggregateCMat(cps, crypto.CipherMatrix{v}); len(rows) > 0 {
			return rows[0]
		}
		return nil
	}
	xtxEnc = mpcObj.Network.AggregateCMat(cps, xtxEnc)
	xtxDuration := ast.metricEnd("null_aggregate_xtx", xtxMark)

	xtyMark := ast.metricMark()
	xtyEnc, _ := crypto.EncryptFloatVector(cps, xtyLocal)
	xtyEnc = aggVec(xtyEnc)
	xtyDuration := ast.metricEnd("null_aggregate_xty", xtyMark)

	ytyMark := ast.metricMark()
	y0ty0Enc, _ := crypto.EncryptFloatVector(cps, []float64{y0ty0Local})
	y0ty0Enc = aggVec(y0ty0Enc)
	ytyDuration := ast.metricEnd("null_aggregate_yty", ytyMark)
	nullAgg := xtxDuration + xtyDuration + ytyDuration
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: encrypt+AggregateCMat(XtX,Xty,y0ty0) %v", nullAgg.Round(time.Millisecond)))
	solveMark := ast.metricMark()

	xtxSS := mpcObj.CMatToSS(cps, rtype, xtxEnc, -1, c, 1, c)
	var xtyCt *rlwe.Ciphertext
	if len(xtyEnc) > 0 {
		xtyCt = xtyEnc[0]
	}
	xtySS := mpcObj.CiphertextToSS(cps, rtype, xtyCt, mpcObj.GetHubPid(), c)
	xtxL, xtxDinv := ast.choleskyFactor(xtxSS)
	// Solve [β̂ | Ω']=(XᵀX)⁻¹[Xᵀy₀ | N·I] as one c×(c+1) RHS. Batching removes a duplicate
	// forward/back-substitution schedule, and cached Ω'=N(XᵀX)⁻¹ removes every per-gene solve.
	rhs := mpc_core.InitRMat(rtype.Zero(), c, c+1)
	for i := 0; i < c; i++ {
		rhs[i][0] = xtySS[i]
	}
	if pid == mpcObj.GetHubPid() {
		nE := rtype.FromFloat64(float64(ast.skatTotalNumInds()), mpcObj.GetFracBits())
		for j := 0; j < c; j++ {
			rhs[j][j+1] = nE
		}
	}
	solvedNull := ast.choleskySolveMat(xtxL, xtxDinv, rhs)
	betaSS := make(mpc_core.RVec, c)
	omp := mpc_core.InitRMat(rtype.Zero(), c, c)
	for i := 0; i < c; i++ {
		betaSS[i] = solvedNull[i][0]
		copy(omp[i], solvedNull[i][1:])
	}
	solveDuration := ast.metricEnd("null_solve", solveMark)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: XtX->SS + Cholesky solve %v", solveDuration.Round(time.Millisecond)))
	betaMark := ast.metricMark()

	// betaRep[ℓ] = β̂_ℓ replicated in every slot (for the score's CPMult). Convert all
	// c rows in one masked SS→CKKS schedule instead of c separate collective conversions.
	var betaRepSS mpc_core.RMat
	if pid > 0 { // pid 0 sits out SSToCMat and needs no c×slots temporary
		slots := cps.GetSlots()
		betaRepSS = make(mpc_core.RMat, c)
		for j := 0; j < c; j++ {
			betaRepSS[j] = mpc_core.InitRVec(betaSS[j], slots)
		}
	}
	betaRep := make(crypto.CipherVector, c)
	if betaRepCM := mpcObj.SSToCMat(cps, betaRepSS); pid > 0 {
		for j := 0; j < c; j++ {
			betaRep[j] = betaRepCM[j][0]
		}
	}
	betaDuration := ast.metricEnd("null_beta_pack", betaMark)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: betaRep SS->CKKS %v", betaDuration.Round(time.Millisecond)))
	rssMark := ast.metricMark()

	// Residual-norm RSS (robust σ̂²): RSS = y₀ᵀy₀ − 2·(Xᵀy₀·β̂) + β̂ᵀ(XᵀX)β̂, 2nd-order in the β̂ error
	// (the plain identity y₀ᵀy₀ − Xᵀy₀·β̂ is 1st-order). All from the c-dim SS aggregates.
	sumRVec := func(v mpc_core.RVec) mpc_core.RElem {
		acc := rtype.Zero()
		for _, x := range v {
			acc = acc.Add(x)
		}
		return acc
	}
	betaCol := make(mpc_core.RMat, c)
	for k := 0; k < c; k++ {
		betaCol[k] = mpc_core.RVec{betaSS[k]}
	}
	xtxBetaMat := mpcObj.SSMultMat(xtxSS, betaCol) // (XᵀX)·β̂
	xtxBeta := make(mpc_core.RVec, c)
	for k := 0; k < c; k++ {
		xtxBeta[k] = xtxBetaMat[k][0]
	}
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	xtxBeta = mpcObj.TruncVec(xtxBeta, db, fb)
	xtyBeta := sumRVec(mpcObj.TruncVec(mpcObj.SSMultElemVec(xtySS, betaSS), db, fb))
	betaXtxBeta := sumRVec(mpcObj.TruncVec(mpcObj.SSMultElemVec(betaSS, xtxBeta), db, fb))
	var y0Ct *rlwe.Ciphertext
	if len(y0ty0Enc) > 0 {
		y0Ct = y0ty0Enc[0]
	}
	y0ty0SS := mpcObj.CiphertextToSS(cps, rtype, y0Ct, mpcObj.GetHubPid(), 1)[0]
	rssSS := y0ty0SS.Sub(xtyBeta).Sub(xtyBeta).Add(betaXtxBeta)
	rssDuration := ast.metricEnd("null_rss", rssMark)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: RSS %v", rssDuration.Round(time.Millisecond)))

	null = skatNull{betaRep: betaRep, betaSS: betaSS, omp: omp, rssSS: rssSS, c: c}
	if ln.X != nil {
		X, y0 = ln.X, ln.Y0
	}
	return
}
