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

func SKATPlainLowRank(G, X *mat.Dense, y []float64) SKATPlainResult {
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
