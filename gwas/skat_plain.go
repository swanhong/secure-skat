package gwas

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

// skat_plain.go — the plaintext L1 oracle for the low-rank SKAT/Burden statistic
// (quantitative / linear trait). Pure gonum, no crypto. This is the reference the
// secure Go pipeline is tested against (Go-secure ≈ L1 ≈ R::SKAT).
//
// Low-rank identities:
//
//	β̂   = (XᵀX)⁻¹ Xᵀy
//	RSS = yᵀy − (Xᵀy)ᵀ β̂          (orthogonality identity — no n-dim residual)
//	s   = Gᵀy − (GᵀX) β̂           (residualized genotype·phenotype score)
//	σ̂²  = RSS / dof,  dof = n − c  (intercept counted in c = NumCov + 1)
//	w_j = 25·max(p_j, 1−p_j)^24,  p_j = colSum(G)_j / (2n)
//	Q       = Σ_j w_j² s_j² / (2σ̂²)
//	Burden  = (Σ_j w_j s_j)²  / (2σ̂²)
//
// Scores are oriented to the minor allele (s_j→−s_j when p_j>0.5) so Burden
// matches R::SKAT's convention; Q is orientation-invariant. X must already
// contain the intercept column.

type SKATPlainResult struct {
	Beta   []float64 // c, null-model coefficients
	RSS    float64
	Dof    int
	Sigma2 float64
	Score  []float64 // m, per-variant score s
	Weight []float64 // m, per-variant weight w
	Q      float64
	Burden float64
}

func SKATPlain(G, X *mat.Dense, y []float64) SKATPlainResult {
	n, c := X.Dims()
	_, m := G.Dims()
	yv := mat.NewVecDense(n, y)

	// β̂ = (XᵀX)⁻¹ Xᵀy
	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	var Xty mat.VecDense
	Xty.MulVec(X.T(), yv)
	var beta mat.VecDense
	if err := beta.SolveVec(&XtX, &Xty); err != nil {
		panic(err)
	}

	// RSS = yᵀy − (Xᵀy)ᵀ β̂  (orthogonality identity)
	yty := mat.Dot(yv, yv)
	rss := yty - mat.Dot(&Xty, &beta)
	dof := n - c
	sigma2 := rss / float64(dof)

	// s = Gᵀy − (GᵀX) β̂
	var Gty mat.VecDense
	Gty.MulVec(G.T(), yv) // m
	var GtX mat.Dense
	GtX.Mul(G.T(), X) // m × c
	var GtXbeta mat.VecDense
	GtXbeta.MulVec(&GtX, &beta) // m
	score := make([]float64, m)

	// weights from genotype column sums
	weight := make([]float64, m)
	var sumWS, sumW2S2 float64
	for j := 0; j < m; j++ {
		score[j] = Gty.AtVec(j) - GtXbeta.AtVec(j)

		colSum := 0.0
		for i := 0; i < n; i++ {
			colSum += G.At(i, j)
		}
		pbar := colSum / (2 * float64(n))
		maf := math.Min(pbar, 1-pbar)
		weight[j] = 25 * math.Pow(1-maf, 24)

		// Orient the score to the minor allele (R::SKAT convention): flipping a
		// variant sends G_ij→2−G_ij, hence s_j→−s_j (Σ residuals = 0). Q (Σw²s²)
		// is unaffected; Burden (Σw s)² is sign-sensitive, so this matters there.
		if pbar > 0.5 {
			score[j] = -score[j]
		}

		sumWS += weight[j] * score[j]
		sumW2S2 += weight[j] * weight[j] * score[j] * score[j]
	}

	scale := 1.0 / (2 * sigma2)
	betaOut := make([]float64, c)
	for i := 0; i < c; i++ {
		betaOut[i] = beta.AtVec(i)
	}

	return SKATPlainResult{
		Beta:   betaOut,
		RSS:    rss,
		Dof:    dof,
		Sigma2: sigma2,
		Score:  score,
		Weight: weight,
		Q:      sumW2S2 * scale,
		Burden: sumWS * sumWS * scale,
	}
}

// --- federated SKAT with party-private variants: plaintext oracle ---
//
// Plaintext mirror of the secure ComputeSKATFederatedPrivate (skat.go); equals the pooled
// single-cohort SKAT Q. See that doc for the design.

// FedParty is one party's plaintext data. Variants are columns of G, each labelled by Gene and
// Role ∈ {"shared","public_only","private"}; ID aligns shared variants across parties.
type FedParty struct {
	G    *mat.Dense // n × m genotypes (dosages)
	X    *mat.Dense // n × c design (intercept + covariates)
	Y    []float64  // n phenotype
	ID   []string   // m variant ids
	Gene []string   // m gene per variant
	Role []string   // m role
}

// SKATFederatedPrivate returns per-gene Q. pub.ID is the public list (shared + public_only);
// priv.ID is shared + private. Shared variants share the same ID across parties.
func SKATFederatedPrivate(pub, priv FedParty) map[string]float64 {
	np, c := pub.X.Dims()
	nq, _ := priv.X.Dims()
	N := float64(np + nq)

	// ---- shared null model: pooled aggregates (covariates only) ----
	var XtX, tmp mat.Dense
	XtX.Mul(pub.X.T(), pub.X)
	tmp.Mul(priv.X.T(), priv.X)
	XtX.Add(&XtX, &tmp) // Xp'Xp + Xq'Xq

	var xtyP, xtyQ mat.VecDense
	xtyP.MulVec(pub.X.T(), mat.NewVecDense(np, pub.Y))
	xtyQ.MulVec(priv.X.T(), mat.NewVecDense(nq, priv.Y))
	xty := mat.NewVecDense(c, nil)
	xty.AddVec(&xtyP, &xtyQ) // Xp'yp + Xq'yq

	var beta mat.VecDense
	if err := beta.SolveVec(&XtX, xty); err != nil {
		panic(err)
	}
	rp := fedResidual(pub.X, pub.Y, &beta)
	rq := fedResidual(priv.X, priv.Y, &beta)
	rss := fedDot(pub.Y, pub.Y) + fedDot(priv.Y, priv.Y) - mat.Dot(xty, &beta)
	sigma2 := rss / float64(np+nq-c)

	qCol := make(map[string]int, len(priv.ID))
	for j, id := range priv.ID {
		qCol[id] = j
	}

	Q := make(map[string]float64)
	add := func(gene string, score, count float64) {
		p := count / (2 * N)
		w := 25 * math.Pow(1-math.Min(p, 1-p), 24)
		Q[gene] += w * w * score * score / (2 * sigma2)
	}

	// ---- PART A: secure over the public list (the private party adds its shared cols) ----
	for k := range pub.ID {
		score := fedColDot(pub.G, k, rp)
		count := fedColSum(pub.G, k)
		if pub.Role[k] == "shared" {
			j := qCol[pub.ID[k]]
			score += fedColDot(priv.G, j, rq) // private party's local contribution (federated sum)
			count += fedColSum(priv.G, j)
		}
		add(pub.Gene[k], score, count)
	}
	// ---- PART B: private variants, local to the private party ----
	for j := range priv.ID {
		if priv.Role[j] == "private" {
			add(priv.Gene[j], fedColDot(priv.G, j, rq), fedColSum(priv.G, j))
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
