// Package utils implements various helper functions.
package utils

import (
	"math/bits"
	"reflect"

	"golang.org/x/exp/constraints"
)

// Min returns the minimum value of the two inputs.
func Min[V constraints.Ordered](a, b V) (r V) {
	if a <= b {
		return a
	}
	return b
}

// Max returns the maximum value of the two inputs.
func Max[V constraints.Ordered](a, b V) (r V) {
	if a >= b {
		return a
	}
	return b
}

// IsNil returns true either type or value are nil.
// Only interfaces or pointers to objects should be passed as argument.
func IsNil(i interface{}) bool {
	return i == nil || reflect.ValueOf(i).IsNil()
}

// BitReverse64 returns the bit-reverse value of the input value, within a context of 2^bitLen.
func BitReverse64[V uint64 | uint32 | int | int64](index V, bitLen int) uint64 {
	return bits.Reverse64(uint64(index)) >> (64 - bitLen)
}

// HammingWeight64 returns the hamming weight if the input value.
func HammingWeight64[V uint64 | uint32 | int | int64](x V) V {
	y := uint64(x)
	y -= (y >> 1) & 0x5555555555555555
	y = (y & 0x3333333333333333) + ((y >> 2) & 0x3333333333333333)
	y = (y + (y >> 4)) & 0x0f0f0f0f0f0f0f0f
	return V(((y * 0x0101010101010101) & 0xffffffffffffffff) >> 56)
}

// AllDistinct returns true if all elements in s are distinct, and false otherwise.
func AllDistinct[V comparable](s []V) bool {
	m := make(map[V]struct{}, len(s))
	for _, si := range s {
		if _, exists := m[si]; exists {
			return false
		}
		m[si] = struct{}{}
	}
	return true
}

// GCD computes the greatest common divisor between a and b.
func GCD[V uint64 | uint32 | int | int64](a, b V) V {
	if a == 0 || b == 0 {
		return 0
	}
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

// ModInv computes the modular inverse of a modulo m using the extended Euclidean algorithm.
// Returns x such that (a * x) % m == 1, or 0 if no inverse exists.
func ModInv(a, m uint64) uint64 {
	if m == 1 {
		return 0
	}

	m0 := m
	x0, x1 := uint64(0), uint64(1)

	a = a % m
	for a > 1 {
		// q is quotient
		q := a / m

		m, a = a%m, m
		x0, x1 = x1-q*x0, x0
	}

	// Make x1 positive
	if int64(x1) < 0 {
		x1 += m0
	}

	return x1
}

// SolveTwoCRTs solves the Chinese Remainder Theorem for two moduli.
// Returns x such that x ≡ a1 (mod m1) and x ≡ a2 (mod m2).
func SolveTwoCRTs(a1, m1, a2, m2 uint64) uint64 {
	// Ensure inputs are reduced
	a1 = a1 % m1
	a2 = a2 % m2

	// Find modular inverse of m1 modulo m2
	inv := ModInv(m1%m2, m2)

	// Calculate t = ((a2 - a1) % m2) * inv % m2
	var t uint64
	if a2 >= a1 {
		t = ((a2 - a1) % m2) * inv % m2
	} else {
		// Handle underflow: a2 - a1 is negative
		t = ((m2 - ((a1 - a2) % m2)) % m2) * inv % m2
	}

	// x = a1 + m1 * t
	x := a1 + m1*t

	// Return x mod (m1 * m2)
	return x % (m1 * m2)
}
