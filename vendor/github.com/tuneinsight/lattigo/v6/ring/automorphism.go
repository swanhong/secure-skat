package ring

import (
	"fmt"
	"math/bits"
	"unsafe"

	"github.com/tuneinsight/lattigo/v6/utils"
)

// Import rlwe for CRT functions (will need to handle circular dependency)
// For now, we'll reimplement minimal versions here to avoid circular import

// AutomorphismNTTIndex computes the look-up table for the automorphism X^{i} -> X^{i*k mod NthRoot}.
func AutomorphismNTTIndex(N int, NthRoot, GalEl uint64) (index []uint64, err error) {

	if N&(N-1) != 0 {
		return nil, fmt.Errorf("N must be a power of two")
	}

	if NthRoot&(NthRoot-1) != 0 {
		return nil, fmt.Errorf("NthRoot must be w power of two")
	}

	var mask, tmp1, tmp2 uint64
	logNthRoot := int(bits.Len64(NthRoot-1) - 1)
	mask = NthRoot - 1
	index = make([]uint64, N)

	for i := 0; i < N; i++ {
		tmp1 = 2*utils.BitReverse64(i, logNthRoot) + 1
		tmp2 = ((GalEl * tmp1 & mask) - 1) >> 1
		index[i] = utils.BitReverse64(tmp2, logNthRoot)
	}

	return
}

// automorphismNTTIndex3N computes the permutation table for 3N ring automorphism in NTT domain.
// NTT form stores evaluations at roots ω^(exponent[i]).
// Automorphism X → X^k maps ω^e → ω^(k*e), which permutes the evaluation points.
func automorphismNTTIndex3N(N int, NthRoot, GalEl uint64) (index []uint64, err error) {
	M := NthRoot // For 3N rings, NthRoot = 3N

	// Get the NTT exponent vector (which roots are used in NTT)
	// For N = 2^a * 3^b, we need the exponent vector from the FFT
	exponents, err := getNTTExponents3N(N)
	if err != nil {
		return nil, err
	}

	// Create exponent-to-index map for fast lookup
	expToIdx := make(map[uint64]int)
	for i, exp := range exponents {
		expToIdx[exp] = i
	}

	// Compute permutation table
	index = make([]uint64, N)
	for i := 0; i < N; i++ {
		// Original exponent at position i
		origExp := exponents[i]

		// After automorphism X → X^GalEl, exponent becomes GalEl * origExp
		newExp := (GalEl * origExp) % M

		// Find which position has this new exponent
		newIdx, ok := expToIdx[newExp]
		if !ok {
			return nil, fmt.Errorf("automorphism permutation error: exponent %d not found", newExp)
		}

		index[i] = uint64(newIdx)
	}

	return index, nil
}

// getNTTExponents3N returns the exponent vector for 3N-NTT.
// For N = 2^a * 3^b, returns the exponents [e₀, e₁, ..., e_{N-1}] such that
// NTT evaluates polynomial at roots [ω^e₀, ω^e₁, ..., ω^e_{N-1}] where ω = e^(2πi/3N).
func getNTTExponents3N(N int) ([]uint64, error) {
	// Factor N as 2^a * 3^b
	a := 0
	tempN := N
	for tempN%2 == 0 {
		a++
		tempN /= 2
	}

	b := 0
	for tempN%3 == 0 {
		b++
		tempN /= 3
	}

	if tempN != 1 {
		return nil, fmt.Errorf("N must be of form 2^a * 3^b")
	}

	M := 3 * N

	// Compute g1, g2 using CRT (same as FFT)
	m1 := 1 << a      // 2^a
	m2 := Pow3(b + 1) // 3^(b+1)

	// g1 ≡ 5 (mod 2^a), g1 ≡ 1 (mod 3^(b+1))
	g1 := CRTPairSimple(5%m1, m1, 1%m2, m2) % M

	// g2 ≡ -1 (mod 2^a), g2 ≡ 4 (mod 3^(b+1))
	g2 := CRTPairSimple((m1-1)%m1, m1, 4%m2, m2) % M

	// Build exponent vector for NTT
	// NTT uses N evaluation points (first N elements of the full exponent matrix)
	rows := 1 << (a + 1) // 2^(a+1)
	cols := 2 * Pow3(b)  // 2 * 3^b
	// Note: rows * cols = 2N (full exponent matrix size)

	// We only use the first N exponents for NTT
	exponents := make([]uint64, N)
	idx := 0
	for i := 0; i < rows && idx < N; i++ {
		g1_pow := ModExp(uint64(g1), uint64(i), uint64(M))
		for j := 0; j < cols && idx < N; j++ {
			g2_pow := ModExp(uint64(g2), uint64(j), uint64(M))
			exp := (g1_pow * g2_pow) % uint64(M)
			exponents[idx] = exp
			idx++
		}
	}

	return exponents, nil
}

// AutomorphismNTT applies the automorphism X^{i} -> X^{i*gen} on a polynomial in the NTT domain.
// Uses NTT domain permutation to preserve Montgomery form.
// It must be noted that the result cannot be in-place.
func (r Ring) AutomorphismNTT(polIn Poly, gen uint64, polOut Poly) {
	// Use NTT domain permutation for both 3N and power-of-2 rings
	// This preserves Montgomery form!
	var index []uint64
	var err error

	if r.Type() == Matrix {
		// 3N ring: use AutomorphismNTTIndex3N
		index, err = r.AutomorphismNTTIndex3N(gen)
	} else {
		// Power-of-2: use AutomorphismNTTIndex
		index, err = AutomorphismNTTIndex(r.N(), r.NthRoot(), gen)
	}

	if err != nil {
		panic(err)
	}

	r.AutomorphismNTTWithIndex(polIn, index, polOut)
}

// AutomorphismNTTIndex3N computes the permutation table for 3N ring automorphism
// using empirical method to match factorized NTT ordering.
//
// Complexity: O(N log N) per galEl (dominated by NTT operations)
// This is necessary because factorized NTT uses its own slot ordering.
func (r Ring) AutomorphismNTTIndex3N(galEl uint64) ([]uint64, error) {
	N := r.N()

	// Get NTT transformer for level 0
	ntt0, ok := r.SubRings[0].GetNTT3N()
	if !ok {
		return nil, fmt.Errorf("ring is not using NumberTheoreticTransformer3N")
	}

	q := ntt0.p

	// Method: Use X^1 polynomial to probe factorized NTT ordering (Python approach)
	// This avoids relying on a separate canonical evaluation-point table.

	// Step 1: Precompute factorized NTT index ordering (O(N log N), once per galEl)
	// Use X^1 = [0, 1, 0, 0, ...] to get evaluation points directly
	coeffsX := r.NewPoly()
	coeffsX.Coeffs[0][1] = 1 // X^1

	nttX := r.NewPoly()
	r.NTT(coeffsX, nttX)

	// nttX[i] is the evaluation point at slot i in factorized ordering.
	evalToSlot := make(map[uint64]int, N)
	for slot := 0; slot < N; slot++ {
		evalToSlot[nttX.Coeffs[0][slot]%q] = slot
	}

	// Step 2: Compute automorphism permutation in O(N)
	index := make([]uint64, N)
	for outSlot := 0; outSlot < N; outSlot++ {
		evalPt := nttX.Coeffs[0][outSlot] % q

		// After automorphism: evalPt -> evalPt^galEl
		evalPtK := ModExp(evalPt, galEl, q)

		inSlot, ok := evalToSlot[evalPtK]
		if !ok {
			return nil, fmt.Errorf("automorphism permutation error: evaluation point %d not found", evalPtK)
		}

		// nttAfter[outSlot] = nttBefore[inSlot]
		index[outSlot] = uint64(inSlot)
	}

	return index, nil
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// AutomorphismNTTWithIndex applies the automorphism X^{i} -> X^{i*gen} on a polynomial in the NTT domain.
// `index` is the lookup table storing the mapping of the automorphism.
// It must be noted that the result cannot be in-place.
func (r Ring) AutomorphismNTTWithIndex(polIn Poly, index []uint64, polOut Poly) {

	level := r.level

	N := r.N()

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(index)%8 != 0  */
		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := 0; i < level+1; i++ {

			/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(polOut.Coeffs)%8 != 0 */
			z := (*[8]uint64)(unsafe.Pointer(&polOut.Coeffs[i][j]))
			y := polIn.Coeffs[i]

			z[0] = y[x[0]]
			z[1] = y[x[1]]
			z[2] = y[x[2]]
			z[3] = y[x[3]]
			z[4] = y[x[4]]
			z[5] = y[x[5]]
			z[6] = y[x[6]]
			z[7] = y[x[7]]
		}
	}
}

// AutomorphismNTTWithIndexThenAddLazy applies the automorphism X^{i} -> X^{i*gen} on a polynomial in the NTT domain .
// `index` is the lookup table storing the mapping of the automorphism.
// The result of the automorphism is added on polOut.
func (r Ring) AutomorphismNTTWithIndexThenAddLazy(polIn Poly, index []uint64, polOut Poly) {

	level := r.level

	N := r.N()

	for j := 0; j < N; j = j + 8 {

		/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(index)%8 != 0 */
		x := (*[8]uint64)(unsafe.Pointer(&index[j]))

		for i := 0; i < level+1; i++ {

			/* #nosec G103 -- behavior and consequences well understood, possible buffer overflow if len(polOut.Coeffs)%8 != 0 */
			z := (*[8]uint64)(unsafe.Pointer(&polOut.Coeffs[i][j]))
			y := polIn.Coeffs[i]

			z[0] += y[x[0]]
			z[1] += y[x[1]]
			z[2] += y[x[2]]
			z[3] += y[x[3]]
			z[4] += y[x[4]]
			z[5] += y[x[5]]
			z[6] += y[x[6]]
			z[7] += y[x[7]]
		}
	}
}

// Automorphism applies the automorphism X^{i} -> X^{i*gen} on a polynomial outside of the NTT domain.
// It must be noted that the result cannot be in-place.
func (r Ring) Automorphism(polIn Poly, gen uint64, polOut Poly) {

	var mask, index, indexRaw, logN, tmp uint64

	/* #nosec G115 -- N cannot be negative */
	N := uint64(r.N())

	level := r.level

	if r.Type() == Matrix {
		r.Automorphism3N(polIn, gen, polOut)
		return
	}

	if r.Type() == ConjugateInvariant {

		mask = 2*N - 1

		/* #nosec G115 -- bitsize cannot be negative */
		logN = uint64(bits.Len64(mask))

		// TODO: find a more efficient way to do
		// the automorphism on Z[X+X^-1]
		for i := uint64(0); i < 2*N; i++ {

			indexRaw = i * gen

			index = indexRaw & mask

			tmp = (indexRaw >> logN) & 1

			// Only consider i -> index if within [0, N-1]
			if index < N {

				idx := i

				// If the starting index is within [N, 2N-1]
				if idx >= N {
					idx = 2*N - idx // Wrap back between [0, N-1]
					tmp ^= 1        // Negate
				}

				for j, s := range r.SubRings[:level+1] {
					polOut.Coeffs[j][index] = polIn.Coeffs[j][idx]*(tmp^1) | (s.Modulus-polIn.Coeffs[j][idx])*tmp
				}
			}
		}

	} else {

		mask = N - 1

		/* #nosec G115 -- bitsize cannot be negative */
		logN = uint64(bits.Len64(mask))

		for i := uint64(0); i < N; i++ {

			indexRaw = i * gen

			index = indexRaw & mask

			tmp = (indexRaw >> logN) & 1

			for j, s := range r.SubRings[:level+1] {
				polOut.Coeffs[j][index] = polIn.Coeffs[j][i]*(tmp^1) | (s.Modulus-polIn.Coeffs[j][i])*tmp
			}
		}
	}
}

// Automorphism3N applies the automorphism X^{i} -> X^{i*gen} on a polynomial for 3N rings.
// For Matrix CKKS with modular polynomial X^N - X^(N/2) + 1 = 0.
// The result cannot be in-place.
func (r Ring) Automorphism3N(polIn Poly, gen uint64, polOut Poly) {

	N := r.N()
	M := 3 * N // For 3N rings, M = 3N
	half := N / 2
	level := r.level

	// Zero out the output polynomial
	for i := 0; i <= level; i++ {
		for j := 0; j < N; j++ {
			polOut.Coeffs[i][j] = 0
		}
	}

	// Apply automorphism: p(X) -> p(X^gen)
	for srcIdx := 0; srcIdx < N; srcIdx++ {
		// Compute (gen * srcIdx) mod (3N)
		newIdxOrig := (uint64(srcIdx) * gen) % uint64(M)

		for i := 0; i <= level; i++ {
			coeff := polIn.Coeffs[i][srcIdx]
			if coeff == 0 {
				continue
			}

			modulus := r.SubRings[i].Modulus
			newIdx := newIdxOrig

			// Reduce using modular polynomial X^N - X^(N/2) + 1 = 0
			// Step 1: Reduce using X^(3N/2) = -1
			if newIdx >= uint64(M/2) {
				newIdx -= uint64(M / 2)
				coeff = modulus - coeff // negate
			}

			// Step 2: Reduce using X^N = X^(N/2) - 1
			q := newIdx / uint64(N)
			rem := newIdx % uint64(N)

			bredconst := r.SubRings[i].BRedConstant

			// Precompute -coeff to avoid overflow in (polOut + modulus - coeff)
			negCoeff := modulus - coeff

			switch q {
			case 0:
				// X^rem - stays as is
				polOut.Coeffs[i][rem] = BRedAdd(polOut.Coeffs[i][rem]+coeff, modulus, bredconst)

			case 1:
				// X^(N+rem) = X^rem * (X^(N/2) - 1) = X^(rem+N/2) - X^rem
				if rem < uint64(half) {
					// X^(rem+N/2) - X^rem
					polOut.Coeffs[i][rem+uint64(half)] = BRedAdd(polOut.Coeffs[i][rem+uint64(half)]+coeff, modulus, bredconst)
					polOut.Coeffs[i][rem] = BRedAdd(polOut.Coeffs[i][rem]+negCoeff, modulus, bredconst)
				} else {
					// rem >= N/2
					// X^(N+rem) = X^(rem-N/2) * X^N = X^(rem-N/2) * (X^(N/2) - 1) = -X^(rem-N/2)
					polOut.Coeffs[i][rem-uint64(half)] = BRedAdd(polOut.Coeffs[i][rem-uint64(half)]+negCoeff, modulus, bredconst)
				}

			case 2:
				// X^(2N+rem) = (X^N)^2 * X^rem = (X^(N/2) - 1)^2 * X^rem
				//          = (X^N - 2X^(N/2) + 1) * X^rem = -X^(N/2) * X^rem = -X^(rem+N/2)
				if rem < uint64(half) {
					polOut.Coeffs[i][rem+uint64(half)] = BRedAdd(polOut.Coeffs[i][rem+uint64(half)]+negCoeff, modulus, bredconst)
				} else {
					// -X^(rem+N/2) with rem >= N/2
					// = -X^(rem-N/2) * X^N = -X^(rem-N/2) * (X^(N/2) - 1) = -X^rem + X^(rem-N/2)
					polOut.Coeffs[i][rem] = BRedAdd(polOut.Coeffs[i][rem]+negCoeff, modulus, bredconst)
					polOut.Coeffs[i][rem-uint64(half)] = BRedAdd(polOut.Coeffs[i][rem-uint64(half)]+coeff, modulus, bredconst)
				}
			}
		}
	}
}
