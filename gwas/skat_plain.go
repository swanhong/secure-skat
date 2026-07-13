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
	Beta      []float64 // c, null-model coefficients
	RSS       float64
	Dof       int
	Sigma2    float64
	Score     []float64 // m, per-variant score s
	Weight    []float64 // m, per-variant weight w
	Q         float64
	Burden    float64
	BurdenVar float64 // zᵀPz = ŵᵀ(GᵀG − GᵀX(XtX)⁻¹XᵀG)ŵ, the Burden variance factor
	BurdenP   float64 // Burden p-value = erfc(√(T/2)), T = B²/(σ̂²·zᵀPz) ~ χ²₁ (R::SKAT r.corr=1)
	SkatZ     float64 // Wilson-Hilferty pivot z (SKAT p-value; screening); p = ½erfc(z/√2)
	SkatP     float64 // SKAT p-value = ½erfc(z/√2), Q ~ Σλχ²₁ via WH moment match (S4 unneeded, δ=0)
}

// whCleanZ is the Wilson-Hilferty pivot for a PSD mixture Σλχ²₁ using only the first three power sums
// (S4 drops out because δ=0 always for PSD kernels: S3²≤S2·S4 by Cauchy-Schwarz). z ↔ p is monotone;
// p = ½erfc(z/√2). Returns +Inf (⇒ p=0) guard handled by the caller for degenerate S2,S3.
func whCleanZ(Q, S1, S2, S3 float64) float64 {
	u := (Q - S1) * S3 / (S2 * S2)
	h := 2 * S3 * S3 / (9 * S2 * S2 * S2)
	return (math.Cbrt(1+u) - 1 + h) / math.Sqrt(h)
}

// skatMoments returns tr(K), tr(K²), tr(K³) of the SKAT kernel K = ½ D(GᵀPG)D, D=diag(w). Exact
// (plaintext oracle); the secure path estimates tr(K³) by Hutchinson.
func skatMoments(GtPG *mat.Dense, w []float64) (S1, S2, S3 float64) {
	m := len(w)
	K := mat.NewDense(m, m, nil)
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			K.Set(i, j, 0.5*w[i]*w[j]*GtPG.At(i, j))
		}
	}
	S1 = mat.Trace(K)
	var K2 mat.Dense
	K2.Mul(K, K)
	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			kij := K.At(i, j)
			S2 += kij * kij         // ‖K‖_F² = tr(K²)
			S3 += K2.At(i, j) * kij // Σ(K²)_ij K_ji = tr(K³) (K symmetric)
		}
	}
	return
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
	signedW := make([]float64, m) // ŵ = minor-allele-oriented weight, for z = Gŵ (Burden variance)
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
		signedW[j] = weight[j]

		// Orient the score to the minor allele (R::SKAT convention): flipping a
		// variant sends G_ij→2−G_ij, hence s_j→−s_j (Σ residuals = 0). Q (Σw²s²)
		// is unaffected; Burden (Σw s)² is sign-sensitive, so this matters there.
		if pbar > 0.5 {
			score[j] = -score[j]
			signedW[j] = -weight[j]
		}

		sumWS += weight[j] * score[j]
		sumW2S2 += weight[j] * weight[j] * score[j] * score[j]
	}

	scale := 1.0 / (2 * sigma2)
	burden := sumWS * sumWS * scale

	// Burden p-value: zᵀPz = ŵᵀ(GᵀG − GᵀX(XtX)⁻¹XᵀG)ŵ; T = B²/(σ̂²·zᵀPz) = 2·Burden/zᵀPz ~ χ²₁.
	var GtG mat.Dense
	GtG.Mul(G.T(), G) // m×m
	var solved mat.Dense
	if err := solved.Solve(&XtX, GtX.T()); err != nil { // (XtX)⁻¹ XᵀG : c×m
		panic(err)
	}
	var corr mat.Dense
	corr.Mul(&GtX, &solved) // GᵀX(XtX)⁻¹XᵀG : m×m
	var GtPG mat.Dense
	GtPG.Sub(&GtG, &corr)
	sw := mat.NewVecDense(m, signedW)
	var Pz mat.VecDense
	Pz.MulVec(&GtPG, sw)
	burdenVar := mat.Dot(sw, &Pz)
	burdenP := 1.0
	if burdenVar > 0 {
		T := 2 * burden / burdenVar
		burdenP = math.Erfc(math.Sqrt(T / 2))
	}

	// SKAT p-value (Wilson-Hilferty, screening): moments of K=½D(GᵀPG)D → z → p.
	skatQ := sumW2S2 * scale
	S1, S2, S3 := skatMoments(&GtPG, weight)
	skatZ, skatP := 0.0, 1.0 // degenerate gene (S2 or S3 ≈ 0) → p=1
	if S2 > 0 && S3 > 0 {
		skatZ = whCleanZ(skatQ, S1, S2, S3)
		skatP = 0.5 * math.Erfc(skatZ/math.Sqrt2)
	}

	betaOut := make([]float64, c)
	for i := 0; i < c; i++ {
		betaOut[i] = beta.AtVec(i)
	}

	return SKATPlainResult{
		Beta:      betaOut,
		RSS:       rss,
		Dof:       dof,
		Sigma2:    sigma2,
		Score:     score,
		Weight:    weight,
		Q:         skatQ,
		Burden:    burden,
		BurdenVar: burdenVar,
		BurdenP:   burdenP,
		SkatZ:     skatZ,
		SkatP:     skatP,
	}
}

// --- federated SKAT with party-private variants: plaintext oracle ---
//
// Plaintext mirror of the secure ComputeSKATFederatedPrivate; equals the pooled
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

// SKATFederatedPrivate returns per-gene SKAT, Burden, and the Burden p-value. pub.ID is the public
// list (shared + public_only); priv.ID is shared + private. Shared variants share the same ID.
func SKATFederatedPrivate(pub, priv FedParty) (skat, burden, burdenP map[string]float64) {
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

	skat = make(map[string]float64)
	bLin := make(map[string]float64)
	// z = Σⱼ ŵⱼ Gⱼ, the weighted burden collapse, kept per gene per cohort (each cohort accumulates
	// its own rows; zᵀPz is then additive across cohorts). Needed only for the Burden p-value.
	zA := map[string]*mat.VecDense{}
	zB := map[string]*mat.VecDense{}
	getZ := func(mp map[string]*mat.VecDense, gene string, n int) *mat.VecDense {
		v, ok := mp[gene]
		if !ok {
			v = mat.NewVecDense(n, nil)
			mp[gene] = v
		}
		return v
	}
	// add folds one variant into SKAT/Burden and returns its minor-allele-oriented weight ŵ (for z).
	add := func(gene string, score, count float64) float64 {
		p := count / (2 * N)
		w := 25 * math.Pow(1-math.Min(p, 1-p), 24)
		sw := w
		if p > 0.5 { // orient to minor allele: SKAT invariant, Burden (Σw s)² sign-sensitive
			score = -score
			sw = -w
		}
		skat[gene] += w * w * score * score / (2 * sigma2)
		bLin[gene] += w * score
		return sw
	}

	// ---- PART A: secure over the public list (the private party adds its shared cols) ----
	for k := range pub.ID {
		score := fedColDot(pub.G, k, rp)
		count := fedColSum(pub.G, k)
		jShared := -1
		if pub.Role[k] == "shared" {
			jShared = qCol[pub.ID[k]]
			score += fedColDot(priv.G, jShared, rq) // private party's local contribution (federated sum)
			count += fedColSum(priv.G, jShared)
		}
		sw := add(pub.Gene[k], score, count)
		getZ(zA, pub.Gene[k], np).AddScaledVec(getZ(zA, pub.Gene[k], np), sw, pub.G.ColView(k))
		if jShared >= 0 {
			getZ(zB, pub.Gene[k], nq).AddScaledVec(getZ(zB, pub.Gene[k], nq), sw, priv.G.ColView(jShared))
		}
	}
	// ---- PART B: private variants, local to the private party ----
	for j := range priv.ID {
		if priv.Role[j] == "private" {
			sw := add(priv.Gene[j], fedColDot(priv.G, j, rq), fedColSum(priv.G, j))
			getZ(zB, priv.Gene[j], nq).AddScaledVec(getZ(zB, priv.Gene[j], nq), sw, priv.G.ColView(j))
		}
	}

	burden = make(map[string]float64, len(bLin))
	burdenP = make(map[string]float64, len(bLin))
	for g, l := range bLin {
		burden[g] = l * l / (2 * sigma2) // Burden = (Σ w s)² / (2σ̂²)
		// zᵀPz = zᵀz − (Xᵀz)ᵀ(XtX)⁻¹(Xᵀz); zᵀz and Xᵀz are additive across cohorts.
		zz := 0.0
		xtz := mat.NewVecDense(c, nil)
		if za := zA[g]; za != nil {
			zz += mat.Dot(za, za)
			var t mat.VecDense
			t.MulVec(pub.X.T(), za)
			xtz.AddVec(xtz, &t)
		}
		if zb := zB[g]; zb != nil {
			zz += mat.Dot(zb, zb)
			var t mat.VecDense
			t.MulVec(priv.X.T(), zb)
			xtz.AddVec(xtz, &t)
		}
		var sol mat.VecDense
		if err := sol.SolveVec(&XtX, xtz); err != nil {
			panic(err)
		}
		zPz := zz - mat.Dot(xtz, &sol)
		burdenP[g] = 1.0
		if zPz > 0 {
			burdenP[g] = math.Erfc(math.Sqrt(burden[g] / zPz)) // T=2·Burden/zᵀPz, p=erfc(√(T/2))
		}
	}
	return skat, burden, burdenP
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
