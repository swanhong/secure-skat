package gwas

import (
	"fmt"
	"math"
	"os"
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

// Add sums two contractions party-wise (all fields additive).
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
	betaRep  crypto.CipherVector // betaRep[ℓ] = Enc(β̂_ℓ) in every slot (for the score)
	xtyEnc   crypto.CipherVector // global Xᵀy₀ (kept for RSS)
	y0ty0Enc crypto.CipherVector // global y₀ᵀy₀ (kept for RSS)
	c        int
	center   float64 // public centering constant, reused by the score
}

// ridgeRel: tiny Tikhonov ridge on the XᵀX diagonal (hub, once, intercept excluded)
// to keep the inverse finite on a singular design; negligible on a well-conditioned one.
const ridgeRel = 1e-6

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

	// Invert XᵀX in SS (c rows × 1 ct, c data slots), then β̂ = (XᵀX)⁻¹·Xᵀy₀.
	xtxSS := mpcObj.CMatToSS(cps, rtype, xtxEnc, -1, c, 1, c)
	invSS, _ := mpcObj.MatrixInverseSymPos(xtxSS)

	var xtyCt *rlwe.Ciphertext
	if len(xtyEnc) > 0 {
		xtyCt = xtyEnc[0]
	}
	xtySS := mpcObj.CiphertextToSS(cps, rtype, xtyCt, mpcObj.GetHubPid(), c)
	xtyCol := make(mpc_core.RMat, c)
	for i := 0; i < c; i++ {
		xtyCol[i] = mpc_core.RVec{xtySS[i]}
	}
	betaMat := mpcObj.SSMultMat(invSS, xtyCol)
	betaSS := make(mpc_core.RVec, c)
	for i := 0; i < c; i++ {
		betaSS[i] = betaMat[i][0]
	}
	betaSS = mpcObj.TruncVec(betaSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())

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

	return skatNull{betaHat: betaHat, betaRep: betaRep, xtyEnc: xtyEnc, y0ty0Enc: y0ty0Enc, c: c, center: center}
}

// maskFirstSlots zeros slots >= c, keeping 0..c-1 — used before InnerSumAll so the
// (Xᵀy₀)ᵀβ̂ inner product never sums CKKS garbage in the padding slots.
func (ast *AssocTest) maskFirstSlots(v crypto.CipherVector, c int) crypto.CipherVector {
	cps := ast.general.cps
	out := make(crypto.CipherVector, len(v))
	for i := range v {
		if v[i] == nil {
			continue
		}
		out[i] = crypto.MaskTrunc(cps, v[i], c)
	}
	return out
}

// computeNullRSSEnc returns RSS = y₀ᵀy₀ − (Xᵀy₀)ᵀβ̂ as a 1-element CipherVector (the
// orthogonality identity makes this exact without touching n). pid 0 returns nil.
func (ast *AssocTest) computeNullRSSEnc(null skatNull) crypto.CipherVector {
	cps := ast.general.cps
	pid := ast.general.mpcObj[0].GetPid()

	if pid == 0 || len(null.xtyEnc) == 0 || len(null.y0ty0Enc) == 0 ||
		len(null.betaHat) == 0 || null.xtyEnc[0] == nil || null.y0ty0Enc[0] == nil {
		return nil
	}

	xtyEnc, betaHat := alignCipherVectorLevels(cps, null.xtyEnc, null.betaHat)
	prod := crypto.CMult(cps, xtyEnc, betaHat)
	prod = ast.maskFirstSlots(prod, null.c)
	dotCt := crypto.InnerSumAll(cps, prod)

	y0ty0Enc, dotVec := alignCipherVectorLevels(cps, null.y0ty0Enc, crypto.CipherVector{dotCt})
	return crypto.CSub(cps, y0ty0Enc, dotVec)
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

// signedWeight returns the minor-allele-oriented weight ŵ_j = t_j·w_j (t_j=−1 iff p̄_j>½),
// folding the orientation into the weight so ScoreCalculation gives both Q (sign²=1) and the
// R::SKAT-oriented Burden from one vector. Returned in SS.
func (ast *AssocTest) signedWeight(dosageSum []float64, nsnps int) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	fracBits := mpcObj.GetFracBits()

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
// inferred from file size / n (int8 = 1 byte/cell). ponytail: dense load — see warning.md (full-N
// needs streaming; #4 subsamples n so this is fine).
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
	nullRSS = ast.computeNullRSSEnc(null)
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
func (ast *AssocTest) blockStat(b, nsnps int, null skatNull, X *mat.Dense, y0 []float64) (qRawSS, bLinSS mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	pid := mpcObj.GetPid()

	var SBlock crypto.CipherMatrix
	dosage := make([]float64, nsnps)
	if pid > 0 {
		lc := LocalContract(ast.readGenoBlockLocal(b), X, y0)
		SBlock = crypto.CipherMatrix{ast.partyScore(lc.GtX, lc.Gty0, null)}
		dosage = lc.DosageSum
	}

	sAggr := mpcObj.Network.AggregateCMat(cps, SBlock)
	if pid > 0 && len(sAggr) > 0 && len(sAggr[0]) > 0 &&
		mpcObj.Network.CanCollectiveBootstrap(cps, sAggr[0][0].Level()) {
		sAggr = mpcObj.Network.CollectiveBootstrapMat(cps, sAggr, -1)
	}
	weightEnc := mpcObj.SSToCVec(cps, ast.signedWeight(dosage, nsnps))

	var qRaw, bLin crypto.CipherVector
	if pid > 0 && len(sAggr) > 0 {
		qRaw, bLin, _, _, _, _ = ast.ScoreCalculation(sAggr[0], weightEnc)
	}
	return ast.scalarCiphertextToShares(qRaw), ast.scalarCiphertextToShares(bLin)
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
		qRawSS, bLinSS := ast.blockStat(b, nsnps, null, X, y0)
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
		q, bl := ast.blockStat(b, nsnps, null, X, y0)
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

func (ast *AssocTest) zeroPlainVectorForNonDataParty() crypto.PlainVector {
	zeros := make([]float64, ast.general.cps.GetSlots())
	pv, _ := crypto.EncodeFloatVector(ast.general.cps, zeros)
	return pv
}

// weightsCalculation computes the shared beta-density term used by the SKAT score.
// weightsCalculation returns (pVec=1−p̄, pBarVec=p̄, w24=25(1−MAF)^24), all in secret shares.
func (ast *AssocTest) weightsCalculation(dosageSum []float64, nsnps_block int) (mpc_core.RVec, mpc_core.RVec, mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()

	xSumRVec := mpc_core.InitRVec(rtype.Zero(), nsnps_block)

	if pid > 0 {
		for j := 0; j < nsnps_block; j++ {
			xSumRVec[j] = rtype.FromFloat64(dosageSum[j], 0) // exact integer sums
		}
	}

	// Calculate MAF securely
	// Total number of individuals N is perfectly public across all parties, so 2N is public.
	// We can compute p_j = xSum / (2N) securely by multiplying the secret-shared xSumRVec by the plaintext scalar 1/(2N).
	totalIndivs := ast.skatTotalNumInds()
	inv2N := rtype.FromFloat64(1.0/float64(2*totalIndivs), mpcObj.GetFracBits())

	p_j := xSumRVec.Copy()
	p_j.MulScalar(inv2N)
	// p_j inherently possesses mpcObj.GetFracBits() precision now (0 bits * frac bits = frac bits)

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

// ScoreCalculation calculates the final SKAT Score statistic iteratively
// ScoreCalculation returns Q=Σw²s² and Burden=Σw·s (1-elem CipherVectors), plus the
// intermediates S2/w2/w2S2/wS for debug.
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

// privateQRaw returns the private party's local raw Q = Σ w² s² (ciphertext, scalar in slot 0) for
// one gene's private variants: key-free score (β̂ stays encrypted) × plaintext weight, no
// AggregateCMat. Enc(0) for an empty gene. Called only on the private party.
func (ast *AssocTest) privateQRaw(G *mat.Dense, null skatNull, X *mat.Dense, y0 []float64) crypto.CipherVector {
	cps := ast.general.cps

	var m int
	if G != nil {
		_, m = G.Dims()
	}
	if m == 0 {
		z, _ := crypto.EncryptFloatVector(cps, []float64{0})
		return z
	}

	lc := LocalContract(G, X, y0)
	s := ast.partyScore(lc.GtX, lc.Gty0, null)
	wEnc, _ := crypto.EncryptFloatVector(cps, skatBetaWeight(lc.DosageSum, ast.skatTotalNumInds()))

	s, wEnc = alignCipherVectorLevels(cps, s, wEnc)
	q, _, _, _, _, _ := ast.ScoreCalculation(s, wEnc)
	return q
}

// privateBlockStat returns one gene's private-variant raw Q as a 1-elem SS RVec. The private party
// computes it locally (privateQRaw); CiphertextToSS reshares from privatePid (the other party gets
// only the ciphertext, pid 0 sits out). Must run for EVERY gene — the Enc(0) on empty genes keeps
// which genes have private variants indistinguishable.
func (ast *AssocTest) privateBlockStat(G *mat.Dense, null skatNull, X *mat.Dense, y0 []float64, privatePid int) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()

	var ct *rlwe.Ciphertext
	if mpcObj.GetPid() == privatePid {
		if q := ast.privateQRaw(G, null, X, y0); len(q) > 0 {
			ct = q[0]
		}
	}
	return mpcObj.CiphertextToSS(ast.general.cps, rtype, ct, privatePid, 1)
}

// ComputeSKATFederatedPrivate returns per-gene Q (slot b = gene/block b). genoBlocks hold the
// public-list genotypes (the private party's aligned, public_only cols = 0); privateOnly[b] is the
// private party's private-variant genotypes for gene b (only read on privatePid, must be nil or
// length GenoNumBlocks). Caller collectively decrypts the result.
//
// ponytail: privateOnly injected as dense per-gene matrices (real deployment would stream them);
// fine until the private set outgrows memory. Genes absent from the public list (no block) are out
// of scope — the public gene/block set is defined by the public-list party.
func (ast *AssocTest) ComputeSKATFederatedPrivate(privateOnly []*mat.Dense, privatePid int) crypto.CipherVector {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()

	tStart := time.Now()
	null, nullRSS, X, y0 := ast.nullSetup()
	log.LLvl1(fmt.Sprintf("[skat_fed] null model: %v", time.Since(tStart).Round(time.Millisecond)))

	nB := ast.general.config.GenoNumBlocks
	if mpcObj.GetPid() == privatePid && privateOnly != nil && len(privateOnly) != nB {
		panic(fmt.Sprintf("ComputeSKATFederatedPrivate: privateOnly has %d blocks, want %d", len(privateOnly), nB))
	}
	qBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)
	tBlocks := time.Now()
	for b := 0; b < nB; b++ {
		tb := time.Now()
		acc := mpc_core.InitRVec(rtype.Zero(), 1)

		// PART A: secure SKAT over the public list (existing per-block path).
		if nsnps := ast.skatBlockNumSnps(b); nsnps > 0 {
			qA, _ := ast.blockStat(b, nsnps, null, X, y0)
			acc.Add(qA)
		}

		// PART B: private variants for this gene (uniform across all genes).
		var G *mat.Dense
		if mpcObj.GetPid() == privatePid && b < len(privateOnly) {
			G = privateOnly[b]
		}
		acc.Add(ast.privateBlockStat(G, null, X, y0, privatePid))

		qBlockSS[b] = acc[0]
		log.LLvl1(fmt.Sprintf("[skat_fed] block %d/%d: %v", b+1, nB, time.Since(tb).Round(time.Millisecond)))
	}
	log.LLvl1(fmt.Sprintf("[skat_fed] all %d blocks: %v", nB, time.Since(tBlocks).Round(time.Millisecond)))

	// Common 1/(2σ̂²) applied once (linear, distributes over A+B).
	if scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS); ok {
		scaleVec := mpc_core.InitRVec(rtype.Zero(), nB)
		for b := 0; b < nB; b++ {
			scaleVec[b] = scaleSS[0].Copy()
		}
		qBlockSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(qBlockSS, scaleVec), mpcObj.GetDataBits(), mpcObj.GetFracBits())
	}
	log.LLvl1(fmt.Sprintf("[skat_fed] total compute: %v", time.Since(tStart).Round(time.Millisecond)))
	return mpcObj.SSToCVec(cps, qBlockSS)
}
