package gwas

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// skat_mvp_aou_plain.go — plaintext oracle for the MVP+AoU FEDERATED per-gene SKAT Q.
// Mirrors the secure design (and .local/mvp_aou/federated_q.py); the secure implementation
// is tested against this, and this must equal the pooled single-cohort SKAT.
//
// Scenario: MVP publishes its variant list. Per gene the variants split into
//   intersection (MVP∩AoU)  — both parties have them → secure federated score;
//   mvp_only                — only MVP (AoU contributes 0, but they're on MVP's public list);
//   aou_only                — only AoU (computed locally by AoU; MVP never sees them).
// β̂/σ̂² are the pooled null model (covariates only, variant-independent).
//
//	per-gene Q = PART A (secure over MVP's public list) + PART B (aou_only, local to AoU)

// FedParty is one party's plaintext data. Variants are columns of G, each labelled by Gene and
// Role ∈ {"intersection","mvp_only","aou_only"}; ID aligns shared variants across parties.
type FedParty struct {
	G    *mat.Dense // n × m genotypes (dosages)
	X    *mat.Dense // n × c design (intercept + covariates)
	Y    []float64  // n phenotype
	ID   []string   // m variant ids
	Gene []string   // m gene per variant
	Role []string   // m role
}

// SKATFederatedMVPAoU returns per-gene Q. mvp.ID is the public list (intersection + mvp_only);
// aou.ID is intersection + aou_only. Intersection variants share the same ID across parties.
func SKATFederatedMVPAoU(mvp, aou FedParty) map[string]float64 {
	nm, c := mvp.X.Dims()
	na, _ := aou.X.Dims()
	N := float64(nm + na)

	// ---- shared null model: pooled aggregates (covariates only) ----
	var XtX, tmp mat.Dense
	XtX.Mul(mvp.X.T(), mvp.X)
	tmp.Mul(aou.X.T(), aou.X)
	XtX.Add(&XtX, &tmp) // Xm'Xm + Xa'Xa

	var xtyM, xtyA mat.VecDense
	xtyM.MulVec(mvp.X.T(), mat.NewVecDense(nm, mvp.Y))
	xtyA.MulVec(aou.X.T(), mat.NewVecDense(na, aou.Y))
	xty := mat.NewVecDense(c, nil)
	xty.AddVec(&xtyM, &xtyA) // Xm'ym + Xa'ya

	var beta mat.VecDense
	if err := beta.SolveVec(&XtX, xty); err != nil {
		panic(err)
	}
	rm := fedResidual(mvp.X, mvp.Y, &beta)
	ra := fedResidual(aou.X, aou.Y, &beta)
	rss := fedDot(mvp.Y, mvp.Y) + fedDot(aou.Y, aou.Y) - mat.Dot(xty, &beta)
	sigma2 := rss / float64(nm+na-c)

	aCol := make(map[string]int, len(aou.ID))
	for j, id := range aou.ID {
		aCol[id] = j
	}

	Q := make(map[string]float64)
	add := func(gene string, score, count float64) {
		p := count / (2 * N)
		w := 25 * math.Pow(1-math.Min(p, 1-p), 24)
		Q[gene] += w * w * score * score / (2 * sigma2)
	}

	// ---- PART A: secure over MVP's public list (AoU adds its intersection cols) ----
	for k := range mvp.ID {
		score := fedColDot(mvp.G, k, rm)
		count := fedColSum(mvp.G, k)
		if mvp.Role[k] == "intersection" {
			j := aCol[mvp.ID[k]]
			score += fedColDot(aou.G, j, ra) // AoU local contribution (federated sum)
			count += fedColSum(aou.G, j)
		}
		add(mvp.Gene[k], score, count)
	}
	// ---- PART B: aou_only variants, local to AoU ----
	for j := range aou.ID {
		if aou.Role[j] == "aou_only" {
			add(aou.Gene[j], fedColDot(aou.G, j, ra), fedColSum(aou.G, j))
		}
	}
	return Q
}

func fedResidual(X *mat.Dense, y []float64, beta *mat.VecDense) []float64 {
	n, _ := X.Dims()
	var fit mat.VecDense
	fit.MulVec(X, beta)
	r := make([]float64, n)
	for i := range r {
		r[i] = y[i] - fit.AtVec(i)
	}
	return r
}

func fedColDot(G *mat.Dense, k int, v []float64) float64 {
	n, _ := G.Dims()
	s := 0.0
	for i := 0; i < n; i++ {
		s += G.At(i, k) * v[i]
	}
	return s
}

func fedColSum(G *mat.Dense, k int) float64 {
	n, _ := G.Dims()
	s := 0.0
	for i := 0; i < n; i++ {
		s += G.At(i, k)
	}
	return s
}

func fedDot(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s += a[i] * b[i]
	}
	return s
}
