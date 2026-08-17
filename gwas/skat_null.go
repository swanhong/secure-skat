package gwas

import (
	"fmt"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
	"gonum.org/v1/gonum/mat"
)

// LocalContraction holds one party's n-free, gene-dependent aggregates of (G, X, y0).
// Its fields add component-wise across parties to the pooled contraction.
type LocalContraction struct {
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

// localGenotypeContract computes the gene-dependent contraction used by per-gene SKAT.
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

func localGenotypeTerms(G, X, y0 *mat.Dense) (LocalContraction, *mat.Dense) {
	n, _ := X.Dims()
	gn, m := G.Dims()
	if gn != n {
		panic(fmt.Sprintf("localGenotypeTerms: G/X rows differ (%d/%d)", gn, n))
	}
	var gtx mat.Dense
	gtx.Mul(G.T(), X)
	dosage := make([]float64, m)
	for j := range dosage {
		dosage[j] = gtx.At(j, 0)
	}
	terms := LocalContraction{GtX: &gtx, DosageSum: dosage}
	if y0 == nil {
		return terms, nil
	}
	yn, _ := y0.Dims()
	if yn != n {
		panic(fmt.Sprintf("localGenotypeTerms: G/Y0 rows differ (%d/%d)", gn, yn))
	}
	var gty mat.Dense
	gty.Mul(G.T(), y0)
	return terms, &gty
}

type windowLocalContraction struct {
	LocalContraction
	Gty0All *mat.Dense // m×q, selected one column at a time during the phenotype pass
	Gamma   *mat.Dense
	Private *privateGeneLocal
}

func normalizedGram(G *mat.Dense, invN float64) *mat.Dense {
	var gamma mat.Dense
	gamma.Mul(G.T(), G)
	data := gamma.RawMatrix().Data
	for i := range data {
		data[i] *= invN
	}
	return &gamma
}

// computeWindowLocal contracts one manifest window without retaining its genotype matrices.
func (ast *AssocTest) computeWindowLocal(window GeneBatchWindow, X, y0 *mat.Dense, privateOnly []*mat.Dense, privatePid int, needMoments bool) []windowLocalContraction {
	local := make([]windowLocalContraction, len(window.Tiles))
	pid := ast.general.mpcObj[0].GetPid()
	invN := 1.0 / float64(ast.skatTotalNumInds())
	for i, tile := range window.Tiles {
		var G *mat.Dense
		if pid > 0 && tile.Variants > 0 {
			G = orientGenotypeLocal(ast.readGenoBlockLocal(tile.Gene))
			local[i].LocalContraction, local[i].Gty0All = localGenotypeTerms(G, X, y0)
			local[i].Gamma = normalizedGram(G, invN)
		}
		var privateG *mat.Dense
		if pid == privatePid && tile.Gene < len(privateOnly) {
			privateG = orientedGenotypeLocalCopy(privateOnly[tile.Gene])
		}
		gl := &geneLocal{LocalContraction: local[i].LocalContraction, Gloc: G}
		local[i].Private = ast.computePrivateGeneLocalMulti(privateG, X, y0, gl, needMoments)
	}
	return local
}

func phenotypeWindowLocal(local []windowLocalContraction, phenotype int) []windowLocalContraction {
	out := make([]windowLocalContraction, len(local))
	for i := range local {
		out[i] = local[i]
		if local[i].Gty0All != nil {
			out[i].Gty0 = mat.Col(nil, phenotype, local[i].Gty0All)
		}
		if local[i].Private != nil {
			private := *local[i].Private
			if private.Gty0All != nil {
				private.Gty0 = mat.Col(nil, phenotype, private.Gty0All)
			}
			out[i].Private = &private
		}
	}
	return out
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

// --- secure low-rank null model (β̂, RSS from c-dim aggregates only; n never enters) ---

type localNull struct {
	X     *mat.Dense
	Y0    []float64
	XtX   *mat.Dense
	Xty0  []float64
	Y0ty0 float64
}

type localMultiNull struct {
	X    *mat.Dense // n×c
	Y0   *mat.Dense // n×q
	XtX  *mat.Dense // c×c
	Xty0 *mat.Dense // c×q
	Yty0 []float64  // q
}

func localMultiNullEquations(cov, pheno *mat.Dense, center float64, q int) localMultiNull {
	if cov == nil || pheno == nil {
		if cov != nil || pheno != nil {
			panic("localMultiNullEquations: covariates and phenotype must both be present")
		}
		return localMultiNull{}
	}
	cn, ncov := cov.Dims()
	pn, pc := pheno.Dims()
	if cn != pn || pc != q || q < 1 {
		panic(fmt.Sprintf("localMultiNullEquations: covariate/phenotype dimensions differ (%dx%d/%dx%d, q=%d)", cn, ncov, pn, pc, q))
	}

	X := mat.NewDense(cn, ncov+1, nil)
	Y0 := mat.NewDense(cn, q, nil)
	for i := 0; i < cn; i++ {
		X.Set(i, 0, 1)
		for j := 0; j < ncov; j++ {
			X.Set(i, j+1, cov.At(i, j))
		}
		for t := 0; t < q; t++ {
			Y0.Set(i, t, pheno.At(i, t)-center)
		}
	}
	var xtx, xty mat.Dense
	xtx.Mul(X.T(), X)
	xty.Mul(X.T(), Y0)
	yty := make([]float64, q)
	for t := 0; t < q; t++ {
		column := mat.Col(nil, t, Y0)
		yty[t] = mat.Dot(mat.NewVecDense(cn, column), mat.NewVecDense(cn, column))
	}
	return localMultiNull{X: X, Y0: Y0, XtX: &xtx, Xty0: &xty, Yty0: yty}
}

func denseRows(matrix *mat.Dense, rows, columns int) [][]float64 {
	out := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		out[i] = make([]float64, columns)
		if matrix != nil {
			for j := 0; j < columns; j++ {
				out[i][j] = matrix.At(i, j)
			}
		}
	}
	return out
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
	c := len(A)
	mul := func(x, y mpc_core.RElem) mpc_core.RElem { // SS scalar product, back to fracBits
		return ast.ssMul(mpc_core.RVec{x}, mpc_core.RVec{y})[0]
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
	c := len(dInv)
	R := 0
	if len(B) > 0 {
		R = len(B[0])
	}
	smul := func(s mpc_core.RElem, row mpc_core.RVec) mpc_core.RVec { // secret scalar × length-R row
		return ast.ssMul(mpc_core.InitRVec(s, R), row)
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

// nullSetupMulti factors pooled XᵀX once and solves all phenotype RHS columns together.
func (ast *AssocTest) nullSetupMulti() (nulls []skatNull, X, y0 *mat.Dense) {
	cps := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()

	c := ast.general.gwasParams.NumCov() + 1
	q := ast.general.config.NumPhenos
	if q == 0 {
		q = 1
	}
	if c > cps.GetSlots() {
		panic(fmt.Sprintf("nullSetupMulti: c=%d exceeds slots=%d", c, cps.GetSlots()))
	}

	center := 0.0
	if ast.general.config.BinaryPheno {
		center = 0.5 // binary {0,1} phenotype; absorbed by the intercept (Q-invariant)
	}
	localMark := ast.metricMark()
	ln := localMultiNullEquations(ast.general.cov, ast.general.pheno, center, q)
	xtxLocal := denseRows(ln.XtX, c, c)
	xtyLocal := denseRows(ln.Xty0, c, q)
	ytyLocal := make([]float64, q)
	copy(ytyLocal, ln.Yty0)

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
	xtyEnc, _, _, err := crypto.EncryptFloatMatrixRow(cps, xtyLocal)
	if err != nil {
		panic(err)
	}
	xtyEnc = mpcObj.Network.AggregateCMat(cps, xtyEnc)
	xtyDuration := ast.metricEnd("null_aggregate_xty", xtyMark)

	ytyMark := ast.metricMark()
	ytyEnc, _ := crypto.EncryptFloatVector(cps, ytyLocal)
	ytyEnc = aggVec(ytyEnc)
	ytyDuration := ast.metricEnd("null_aggregate_yty", ytyMark)
	nullAgg := xtxDuration + xtyDuration + ytyDuration
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: encrypt+AggregateCMat(XtX,Xty,y0ty0) %v", nullAgg.Round(time.Millisecond)))
	solveMark := ast.metricMark()

	xtxSS := mpcObj.CMatToSS(cps, rtype, xtxEnc, -1, c, 1, c)
	xtySS := mpcObj.CMatToSS(cps, rtype, xtyEnc, -1, c, 1, q)
	xtxL, xtxDinv := ast.choleskyFactor(xtxSS)
	// Solve [β̂₀..β̂q₋₁ | Ω']=(XᵀX)⁻¹[XᵀY₀ | N·I] in one factorization.
	rhs := mpc_core.InitRMat(rtype.Zero(), c, q+c)
	for i := 0; i < c; i++ {
		copy(rhs[i][:q], xtySS[i])
	}
	if pid == mpcObj.GetHubPid() {
		nE := rtype.FromFloat64(float64(ast.skatTotalNumInds()), mpcObj.GetFracBits())
		for j := 0; j < c; j++ {
			rhs[j][q+j] = nE
		}
	}
	solvedNull := ast.choleskySolveMat(xtxL, xtxDinv, rhs)
	betas := make([]mpc_core.RVec, q)
	for t := range betas {
		betas[t] = make(mpc_core.RVec, c)
	}
	omp := mpc_core.InitRMat(rtype.Zero(), c, c)
	for i := 0; i < c; i++ {
		for t := 0; t < q; t++ {
			betas[t][i] = solvedNull[i][t]
		}
		copy(omp[i], solvedNull[i][q:])
	}
	solveDuration := ast.metricEnd("null_solve", solveMark)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: XtX->SS + Cholesky solve %v", solveDuration.Round(time.Millisecond)))
	betaMark := ast.metricMark()

	// Convert all q*c replicated coefficients in one masked SS→CKKS schedule.
	var betaRepSS mpc_core.RMat
	if pid > 0 {
		slots := cps.GetSlots()
		betaRepSS = make(mpc_core.RMat, q*c)
		for t := 0; t < q; t++ {
			for j := 0; j < c; j++ {
				betaRepSS[t*c+j] = mpc_core.InitRVec(betas[t][j], slots)
			}
		}
	}
	betaReps := make([]crypto.CipherVector, q)
	for t := range betaReps {
		betaReps[t] = make(crypto.CipherVector, c)
	}
	if betaRepCM := mpcObj.SSToCMat(cps, betaRepSS); pid > 0 {
		for t := 0; t < q; t++ {
			for j := 0; j < c; j++ {
				betaReps[t][j] = betaRepCM[t*c+j][0]
			}
		}
	}
	betaDuration := ast.metricEnd("null_beta_pack", betaMark)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: betaRep SS->CKKS %v", betaDuration.Round(time.Millisecond)))
	rssMark := ast.metricMark()

	// Residual RSS for every phenotype; pooled XᵀX and its factorization are shared.
	db, fb := mpcObj.GetDataBits(), mpcObj.GetFracBits()
	ytySS := mpcObj.CVecToSS(cps, rtype, ytyEnc, -1, len(ytyEnc), q, 1)
	rss := make(mpc_core.RVec, q)
	for t := 0; t < q; t++ {
		phenotypeStarted := time.Now()
		xtyColumn := make(mpc_core.RVec, c)
		for i := 0; i < c; i++ {
			xtyColumn[i] = xtySS[i][t]
		}
		xtxBeta := mpcObj.TruncVec(ast.ssMatVec(xtxSS, betas[t]), db, fb)
		xtyBeta := ast.ssDot(xtyColumn, betas[t])
		betaXtxBeta := ast.ssDot(betas[t], xtxBeta)
		rss[t] = ytySS[t].Sub(xtyBeta).Sub(xtyBeta).Add(betaXtxBeta)
		fedTimings.phenotypeNullRSS[t] += time.Since(phenotypeStarted)
	}
	rssDuration := ast.metricEnd("null_rss", rssMark)
	log.LLvl1(fmt.Sprintf("[skat_fed]   null: RSS %v", rssDuration.Round(time.Millisecond)))

	nulls = make([]skatNull, q)
	for t := 0; t < q; t++ {
		nulls[t] = skatNull{betaRep: betaReps[t], betaSS: betas[t], omp: omp, rssSS: rss[t], c: c}
	}
	if ln.X != nil {
		X, y0 = ln.X, ln.Y0
	}
	return
}

// nullSetup preserves the scalar API for the ordinary q=1 implementation.
func (ast *AssocTest) nullSetup() (null skatNull, X *mat.Dense, y0 []float64) {
	nulls, X, Y0 := ast.nullSetupMulti()
	if len(nulls) != 1 {
		panic("ordinary SKAT path supports exactly one phenotype")
	}
	if Y0 != nil {
		y0 = mat.Col(nil, 0, Y0)
	}
	return nulls[0], X, y0
}
