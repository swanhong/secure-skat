package ring

import "fmt"

type Dimensions struct {
	Rows, Cols int
}

// EvalPolyModP evaluates y = sum poly[i] * x^{i} mod p.
func EvalPolyModP(x uint64, poly []uint64, p uint64) (y uint64) {
	brc := GenBRedConstant(p)
	y = poly[len(poly)-1]
	for i := len(poly) - 2; i >= 0; i-- {
		y = BRed(y, x, p, brc)
		y = CRed(y+poly[i], p)
	}

	return
}

// Min returns the minimum between to int
func Min(x, y int) int {
	if x > y {
		return y
	}

	return x
}

// ModExp performs the modular exponentiation x^e mod p,
// x and p are required to be at most 64 bits to avoid an overflow.
func ModExp(x, e, p uint64) (result uint64) {
	brc := GenBRedConstant(p)
	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = BRed(result, x, p, brc)
		}
		x = BRed(x, x, p, brc)
	}
	return result
}

// ModExpPow2 performs the modular exponentiation x^e mod p, where p is a power of two,
// x and p are required to be at most 64 bits to avoid an overflow.
func ModExpPow2(x, e, p uint64) (result uint64) {

	result = 1
	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result *= x
		}
		x *= x
	}
	return result & (p - 1)
}

// ModexpMontgomery performs the modular exponentiation x^e mod p,
// where x is in Montgomery form, and returns x^e in Montgomery form.
func ModexpMontgomery(x uint64, e int, q, mredconstant uint64, bredconstant [2]uint64) (result uint64) {

	result = MForm(1, q, bredconstant)

	for i := e; i > 0; i >>= 1 {
		if i&1 == 1 {
			result = MRed(result, x, q, mredconstant)
		}
		x = MRed(x, x, q, mredconstant)
	}
	return result
}

// Pow3 computes 3^b
func Pow3(b int) int {
	result := 1
	for i := 0; i < b; i++ {
		result *= 3
	}
	return result
}

// CRTPairSimple solves the Chinese Remainder Theorem for two coprime moduli.
// Returns x such that x ≡ a1 (mod m1) and x ≡ a2 (mod m2).
func CRTPairSimple(a1, m1, a2, m2 int) int {
	// Extended Euclidean algorithm to find inverse
	inv := invModSimple(m1%m2, m2)
	t := ((a2 - a1) % m2) * inv % m2
	if t < 0 {
		t += m2
	}
	return a1 + m1*t
}

// invModSimple computes modular inverse using Extended Euclidean algorithm
func invModSimple(a, m int) int {
	if a < 0 {
		a = a%m + m
	}
	a = a % m

	t, newT := 0, 1
	r, newR := m, a

	for newR != 0 {
		quotient := r / newR
		t, newT = newT, t-quotient*newT
		r, newR = newR, r-quotient*newR
	}

	if r > 1 {
		panic(fmt.Sprintf("%d is not invertible mod %d", a, m))
	}

	if t < 0 {
		t += m
	}

	return t
}
