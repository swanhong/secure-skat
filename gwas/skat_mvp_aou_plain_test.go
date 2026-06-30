package gwas

import (
	"math"
	"math/rand"
	"testing"

	"gonum.org/v1/gonum/mat"
)

// The federated MVP+AoU per-gene Q (secure over MVP's public list + aou_only local) must equal
// the pooled single-cohort SKAT Q on the union genotype matrix (0-filled for party-unique).

func fedFixture() (mvp, aou FedParty) {
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
	// genes G1 {a:intersection, b:mvp_only, c:aou_only}, G2 {d,e:intersection}
	mvp = FedParty{
		X: makeX(8), G: makeG(8, []float64{0.30, 0.20, 0.25, 0.15}),
		ID: []string{"a", "b", "d", "e"}, Gene: []string{"G1", "G1", "G2", "G2"},
		Role: []string{"intersection", "mvp_only", "intersection", "intersection"},
	}
	mvp.Y = makeY(mvp.X)
	aou = FedParty{
		X: makeX(6), G: makeG(6, []float64{0.30, 0.18, 0.25, 0.15}),
		ID: []string{"a", "c", "d", "e"}, Gene: []string{"G1", "G1", "G2", "G2"},
		Role: []string{"intersection", "aou_only", "intersection", "intersection"},
	}
	aou.Y = makeY(aou.X)
	return
}

// pooledPerGeneQ — the ground truth: build the union genotype (0-fill party-unique), run
// pooled single-cohort SKAT per gene.
func pooledPerGeneQ(mvp, aou FedParty) map[string]float64 {
	nm, c := mvp.X.Dims()
	na, _ := aou.X.Dims()
	N := nm + na
	type uv struct {
		gene             string
		mcol, acol       int
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
	for k, id := range mvp.ID {
		put(id, mvp.Gene[k]).mcol = k
	}
	for j, id := range aou.ID {
		put(id, aou.Gene[j]).acol = j
	}
	M := len(order)
	G := mat.NewDense(N, M, nil)
	for kk, id := range order {
		e := u[id]
		if e.mcol >= 0 {
			for i := 0; i < nm; i++ {
				G.Set(i, kk, mvp.G.At(i, e.mcol))
			}
		}
		if e.acol >= 0 {
			for i := 0; i < na; i++ {
				G.Set(nm+i, kk, aou.G.At(i, e.acol))
			}
		}
	}
	X := mat.NewDense(N, c, nil)
	y := make([]float64, N)
	for i := 0; i < nm; i++ {
		for j := 0; j < c; j++ {
			X.Set(i, j, mvp.X.At(i, j))
		}
		y[i] = mvp.Y[i]
	}
	for i := 0; i < na; i++ {
		for j := 0; j < c; j++ {
			X.Set(nm+i, j, aou.X.At(i, j))
		}
		y[nm+i] = aou.Y[i]
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
	Q := map[string]float64{}
	for kk, id := range order {
		score := fedColDot(G, kk, res)
		count := fedColSum(G, kk)
		p := count / (2 * float64(N))
		w := 25 * math.Pow(1-math.Min(p, 1-p), 24)
		Q[u[id].gene] += w * w * score * score / (2 * sigma2)
	}
	return Q
}

func TestSKATFederatedMVPAoUMatchesPooled(t *testing.T) {
	mvp, aou := fedFixture()
	Qfed := SKATFederatedMVPAoU(mvp, aou)
	Qpool := pooledPerGeneQ(mvp, aou)
	if len(Qfed) != len(Qpool) {
		t.Fatalf("gene count: fed %d pooled %d", len(Qfed), len(Qpool))
	}
	for g, qp := range Qpool {
		if !approxEqual(Qfed[g], qp, 1e-9) {
			t.Errorf("gene %s: federated=%.12g pooled=%.12g", g, Qfed[g], qp)
		}
	}
}
