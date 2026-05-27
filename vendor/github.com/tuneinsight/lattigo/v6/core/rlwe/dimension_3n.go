package rlwe

// Specification of dimension N = 2^Order2 * 2^Order3.
type DimSpec struct {
	Order2 int
	Order3 int
}

// Corresponding dimension of the specification.
func (spec DimSpec) Dimension() int {
	return Compute2a3b(spec.Order2, spec.Order3)
}
