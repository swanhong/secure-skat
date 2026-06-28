package gwas

import (
	"math"
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

func TestSKATPlainLowRankMatchesDirect(t *testing.T) {
	G, X, y := plainFixture()
	wantQ, wantB, wantRSS := directSKAT(G, X, y)

	got := SKATPlainLowRank(G, X, y)

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

func TestSKATPlainLowRankDof(t *testing.T) {
	G, X, y := plainFixture()
	got := SKATPlainLowRank(G, X, y)
	if got.Dof != 6-2 { // n − c, intercept counted in c
		t.Errorf("dof: got %d want %d", got.Dof, 4)
	}
}
