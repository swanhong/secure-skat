package gwas

import "gonum.org/v1/gonum/mat"

// LocalContract and LocalContraction.Add are test-only helpers: they back
// TestSKATLocalPartyAdditivity, which verifies the Σ_party == full-cohort n-independence
// invariant central to the federated design. The shipped pipeline uses only
// localGenotypeContract and the embedded LocalContraction struct.

// LocalContract contracts (G, X, y0) locally; y0 centered, X includes intercept.
func LocalContract(G, X *mat.Dense, y0 []float64) LocalContraction {
	lc := localGenotypeContract(G, X, y0)
	n, _ := X.Dims()
	lc.XtX, lc.Xty0, lc.Y0ty0 = normalEqs(X, mat.NewVecDense(n, y0))
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
