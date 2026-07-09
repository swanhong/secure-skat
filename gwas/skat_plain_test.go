package gwas

import (
	"math"
	"math/rand"
	"testing"

	"gonum.org/v1/gonum/mat"
)

// Go/gonum plaintext L1 oracle for the low-rank SKAT/Burden statistic.
//
// The oracle (skat_plain.go) computes the score via the LOW-RANK identities
//   s   = Gᵀy − (GᵀX)β̂          (not Gᵀ(y − Xβ̂))
//   RSS = y₀ᵀy₀ − (Xᵀy₀)ᵀβ̂       (not ‖y − Xβ̂‖²)
// This test recomputes Q/Burden/RSS the DIRECT (residual) way with gonum — an
// independent code path — and asserts they agree, pinning those identities and
// the 1/(2σ̂²) scaling. (Weight vs R::SKAT is pinned later, L2.)

// small deterministic fixture: n=6, m=3, c=2 (intercept + 1 covariate)
func plainFixture() (G, X *mat.Dense, y []float64) {
	// col0 is ALT-major (Σ=10, p̄=0.83 > 0.5) so it exercises the minor-allele
	// orientation flip; col1 (p̄≈0.42) and col2 (p̄=0.5) are not flipped.
	gRows := [][]float64{
		{2, 1, 2},
		{1, 0, 2},
		{2, 1, 0},
		{2, 0, 1},
		{1, 2, 1},
		{2, 1, 0},
	}
	xCov := []float64{0.5, -1.2, 0.3, 2.0, -0.7, 1.1} // one covariate
	y = []float64{1.0, -0.5, 0.3, 2.2, -1.1, 0.8}

	n := len(y)
	G = mat.NewDense(n, 3, nil)
	X = mat.NewDense(n, 2, nil)
	for i := 0; i < n; i++ {
		for j := 0; j < 3; j++ {
			G.Set(i, j, gRows[i][j])
		}
		X.Set(i, 0, 1.0) // intercept
		X.Set(i, 1, xCov[i])
	}
	return G, X, y
}

// directSKAT computes Q/Burden/RSS the residual way, independent of skat_plain.go.
func directSKAT(G, X *mat.Dense, y []float64) (q, burden, rss float64) {
	n, c := X.Dims()
	_, m := G.Dims()
	yv := mat.NewVecDense(n, y)

	// β̂ = (XᵀX)⁻¹Xᵀy via SolveVec
	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	var Xty mat.VecDense
	Xty.MulVec(X.T(), yv)
	var beta mat.VecDense
	if err := beta.SolveVec(&XtX, &Xty); err != nil {
		panic(err)
	}

	// resid = y − Xβ̂
	var fitted mat.VecDense
	fitted.MulVec(X, &beta)
	resid := make([]float64, n)
	rss = 0
	for i := 0; i < n; i++ {
		resid[i] = y[i] - fitted.AtVec(i)
		rss += resid[i] * resid[i]
	}
	dof := n - c
	sigma2 := rss / float64(dof)

	// s_j = Σ_i G_ij resid_i ; weights from column sums
	residV := mat.NewVecDense(n, resid)
	var s mat.VecDense
	s.MulVec(G.T(), residV)

	var sumWS, sumW2S2 float64
	for j := 0; j < m; j++ {
		colSum := 0.0
		for i := 0; i < n; i++ {
			colSum += G.At(i, j)
		}
		pbar := colSum / (2 * float64(n))
		maf := math.Min(pbar, 1-pbar)
		w := 25 * math.Pow(1-maf, 24)
		sj := s.AtVec(j)
		if pbar > 0.5 { // orient the score to the minor allele (R::SKAT convention)
			sj = -sj
		}
		sumWS += w * sj
		sumW2S2 += w * w * sj * sj
	}
	scale := 1.0 / (2 * sigma2)
	return sumW2S2 * scale, sumWS * sumWS * scale, rss
}

func approxEqual(a, b, tol float64) bool { return math.Abs(a-b) <= tol*(1+math.Abs(b)) }

func TestSKATPlainMatchesDirect(t *testing.T) {
	G, X, y := plainFixture()
	wantQ, wantB, wantRSS := directSKAT(G, X, y)

	got := SKATPlain(G, X, y)

	if !approxEqual(got.RSS, wantRSS, 1e-9) {
		t.Errorf("RSS: got %.12g want %.12g", got.RSS, wantRSS)
	}
	if !approxEqual(got.Q, wantQ, 1e-9) {
		t.Errorf("Q: got %.12g want %.12g", got.Q, wantQ)
	}
	if !approxEqual(got.Burden, wantB, 1e-9) {
		t.Errorf("Burden: got %.12g want %.12g", got.Burden, wantB)
	}
}

func TestSKATPlainDof(t *testing.T) {
	G, X, y := plainFixture()
	got := SKATPlain(G, X, y)
	if got.Dof != 6-2 { // n − c, intercept counted in c
		t.Errorf("dof: got %d want %d", got.Dof, 4)
	}
}

// The federated per-gene Q (secure over the public list + private variants computed locally) must
// equal the pooled single-cohort SKAT Q on the union genotype matrix (0-filled for party-unique).

func fedFixture() (pub, priv FedParty) {
	r := rand.New(rand.NewSource(5))
	const c = 2
	makeX := func(n int) *mat.Dense {
		X := mat.NewDense(n, c, nil)
		for i := 0; i < n; i++ {
			X.Set(i, 0, 1)
			X.Set(i, 1, r.NormFloat64())
		}
		return X
	}
	makeG := func(n int, maf []float64) *mat.Dense {
		G := mat.NewDense(n, len(maf), nil)
		for i := 0; i < n; i++ {
			for j := range maf {
				d := 0.0
				if r.Float64() < maf[j] {
					d++
				}
				if r.Float64() < maf[j] {
					d++
				}
				G.Set(i, j, d)
			}
		}
		return G
	}
	makeY := func(X *mat.Dense) []float64 {
		n, _ := X.Dims()
		y := make([]float64, n)
		for i := range y {
			y[i] = 0.5*X.At(i, 1) + r.NormFloat64()
		}
		return y
	}
	// genes G1 {a:shared, b:public_only, c:private}, G2 {d,e:shared}
	pub = FedParty{
		X: makeX(8), G: makeG(8, []float64{0.30, 0.20, 0.25, 0.15}),
		ID: []string{"a", "b", "d", "e"}, Gene: []string{"G1", "G1", "G2", "G2"},
		Role: []string{"shared", "public_only", "shared", "shared"},
	}
	pub.Y = makeY(pub.X)
	priv = FedParty{
		X: makeX(6), G: makeG(6, []float64{0.30, 0.18, 0.25, 0.15}),
		ID: []string{"a", "c", "d", "e"}, Gene: []string{"G1", "G1", "G2", "G2"},
		Role: []string{"shared", "private", "shared", "shared"},
	}
	priv.Y = makeY(priv.X)
	return
}

// pooledPerGeneQ — the ground truth: build the union genotype (0-fill party-unique), run
// pooled single-cohort SKAT per gene.
func pooledPerGene(pub, priv FedParty) (skat, burden, burdenP map[string]float64) {
	np, c := pub.X.Dims()
	nq, _ := priv.X.Dims()
	N := np + nq
	type uv struct {
		gene       string
		pcol, qcol int
	}
	order := []string{}
	u := map[string]*uv{}
	put := func(id, gene string) *uv {
		if _, ok := u[id]; !ok {
			u[id] = &uv{gene, -1, -1}
			order = append(order, id)
		}
		return u[id]
	}
	for k, id := range pub.ID {
		put(id, pub.Gene[k]).pcol = k
	}
	for j, id := range priv.ID {
		put(id, priv.Gene[j]).qcol = j
	}
	M := len(order)
	G := mat.NewDense(N, M, nil)
	for kk, id := range order {
		e := u[id]
		if e.pcol >= 0 {
			for i := 0; i < np; i++ {
				G.Set(i, kk, pub.G.At(i, e.pcol))
			}
		}
		if e.qcol >= 0 {
			for i := 0; i < nq; i++ {
				G.Set(np+i, kk, priv.G.At(i, e.qcol))
			}
		}
	}
	X := mat.NewDense(N, c, nil)
	y := make([]float64, N)
	for i := 0; i < np; i++ {
		for j := 0; j < c; j++ {
			X.Set(i, j, pub.X.At(i, j))
		}
		y[i] = pub.Y[i]
	}
	for i := 0; i < nq; i++ {
		for j := 0; j < c; j++ {
			X.Set(np+i, j, priv.X.At(i, j))
		}
		y[np+i] = priv.Y[i]
	}
	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	var xty mat.VecDense
	xty.MulVec(X.T(), mat.NewVecDense(N, y))
	var beta mat.VecDense
	if err := beta.SolveVec(&XtX, &xty); err != nil {
		panic(err)
	}
	res := fedResidual(X, y, &beta)
	sigma2 := (fedDot(y, y) - mat.Dot(&xty, &beta)) / float64(N-c)
	skat = map[string]float64{}
	bLin := map[string]float64{}
	zGene := map[string]*mat.VecDense{} // gene -> z = Σ ŵ G_col (N-dim), for the Burden variance
	for kk, id := range order {
		score := fedColDot(G, kk, res)
		count := fedColSum(G, kk)
		p := count / (2 * float64(N))
		w := 25 * math.Pow(1-math.Min(p, 1-p), 24)
		sw := w
		if p > 0.5 {
			score = -score
			sw = -w
		}
		gene := u[id].gene
		skat[gene] += w * w * score * score / (2 * sigma2)
		bLin[gene] += w * score
		z, ok := zGene[gene]
		if !ok {
			z = mat.NewVecDense(N, nil)
			zGene[gene] = z
		}
		z.AddScaledVec(z, sw, G.ColView(kk))
	}
	burden = map[string]float64{}
	burdenP = map[string]float64{}
	for g, l := range bLin {
		burden[g] = l * l / (2 * sigma2)
		z := zGene[g]
		var xtz, sol mat.VecDense
		xtz.MulVec(X.T(), z)
		if err := sol.SolveVec(&XtX, &xtz); err != nil {
			panic(err)
		}
		zPz := mat.Dot(z, z) - mat.Dot(&xtz, &sol)
		burdenP[g] = 1.0
		if zPz > 0 {
			burdenP[g] = math.Erfc(math.Sqrt(burden[g] / zPz))
		}
	}
	return skat, burden, burdenP
}

func TestSKATFederatedPrivateMatchesPooled(t *testing.T) {
	pub, priv := fedFixture()
	skatFed, burdenFed, burdenPFed := SKATFederatedPrivate(pub, priv)
	skatPool, burdenPool, burdenPPool := pooledPerGene(pub, priv)
	if len(skatFed) != len(skatPool) {
		t.Fatalf("gene count: fed %d pooled %d", len(skatFed), len(skatPool))
	}
	for g, sp := range skatPool {
		if !approxEqual(skatFed[g], sp, 1e-9) {
			t.Errorf("gene %s SKAT: federated=%.12g pooled=%.12g", g, skatFed[g], sp)
		}
		if !approxEqual(burdenFed[g], burdenPool[g], 1e-9) {
			t.Errorf("gene %s Burden: federated=%.12g pooled=%.12g", g, burdenFed[g], burdenPool[g])
		}
		if !approxEqual(burdenPFed[g], burdenPPool[g], 1e-9) {
			t.Errorf("gene %s BurdenP: federated=%.12g pooled=%.12g", g, burdenPFed[g], burdenPPool[g])
		}
	}
}
