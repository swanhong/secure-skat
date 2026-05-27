package rlwe

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v6/ring"
)

// CheckModuli checks that the provided q and p correspond to a valid moduli chain.
func CheckModuli3N(q, p []uint64) error {

	for i, qi := range q {
		/* #nosec G115 -- error is returned if integer overflow conversion */
		if uint64(bits.Len64(qi)-1) > MaxModuliSize+1 {
			return fmt.Errorf("a Qi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
		}
	}

	for i, qi := range q {
		if !ring.IsPrime(qi) {
			return fmt.Errorf("a Qi (i=%d) is not a prime", i)
		}
	}

	if p != nil {

		for i, pi := range p {
			/* #nosec G115 -- error is triggered if integer overflow conversion */
			if uint64(bits.Len64(pi)-1) > MaxModuliSize+2 {
				return fmt.Errorf("a Pi bit-size (i=%d) is larger than %d", i, MaxModuliSize)
			}
		}

		for i, pi := range p {
			if !ring.IsPrime(pi) {
				return fmt.Errorf("a Pi (i=%d) is not a prime", i)
			}
		}
	}

	return nil
}

// GenModuli generates a valid moduli chain from the provided moduli sizes.
func GenModuli3N(NthRoot int, logQ, logP []int) (q, p []uint64, err error) {
	// Skip size params check for now
	if err = checkModuliLogSize(logQ, logP); err != nil {
		return
	}

	// Extracts all the different primes bit size and maps their number
	primesbitlen := make(map[int]int)
	for _, qi := range logQ {
		primesbitlen[qi]++
	}

	for _, pj := range logP {
		primesbitlen[pj]++
	}

	// For each bit-size, finds that many 3N-friendly primes
	primes := make(map[int][]uint64)
	for bitsize, value := range primesbitlen {
		// /* #nosec G115 -- bitsize cannot be negative */
		g := ring.NewNTTFriendlyPrimesGenerator3N(uint64(bitsize), uint64(NthRoot))

		if primes[bitsize], err = g.NextAlternatingPrimes(value); err != nil {
			return q, p, fmt.Errorf("cannot GenModuli: failed to generate %d primes of bit-size=%d for NthRoot=%d: %w", value, bitsize, NthRoot, err)
		}

		// Use Find3NRNSPrimes to generate 3N-friendly primes
		// primes[bitsize], err = ring.Find3NRNSPrimes(NthRoot/3, bitsize, value, 1000)
		// if err != nil {
		// 	return q, p, fmt.Errorf("cannot GenModuli3N: %w", err)
		// }
	}

	// Assigns the primes to the moduli chain
	for _, qi := range logQ {
		q = append(q, primes[qi][0])
		primes[qi] = primes[qi][1:]
	}

	// Assigns the primes to the special primes list for the extended ring
	for _, pj := range logP {
		p = append(p, primes[pj][0])
		primes[pj] = primes[pj][1:]
	}

	return
}
