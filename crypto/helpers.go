package crypto

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func Max(x, y int) int {
	if x <= y {
		return y
	}
	return x
}

func Min(a int, b int) int {
	if a > b {
		return b
	}
	return a
}

func Mod(n int, modulus int) int {
	n = n % modulus
	if n < 0 {
		n += modulus
	}
	return n
}

func ConvertVectorComplexToFloat64(v []complex128) []float64 {
	res := make([]float64, len(v))
	for i, el := range v {
		res[i] = real(el)
	}
	return res
}

func PadVector(v []float64, slots int) []float64 {
	if len(v) >= slots {
		return append([]float64(nil), v[:slots]...)
	}
	out := make([]float64, slots)
	copy(out, v)
	return out
}

func FindClosestPow2(n int) int {
	bigPower2 := 1
	for bigPower2 < n {
		bigPower2 *= 2
	}
	return bigPower2
}

func intCeilLog2(n int) int {
	if n <= 1 {
		return 0
	}
	x := 0
	v := 1
	for v < n {
		v <<= 1
		x++
	}
	return x
}

func intCeilSqrt(n int) int {
	if n <= 0 {
		return 0
	}
	x := 1
	for x*x < n {
		x++
	}
	return x
}

func binaryCipherVectorOp(cryptoParams *CryptoParams, X CipherVector, Y CipherVector, op func(eval *ckks.Evaluator, a *rlwe.Ciphertext, b rlwe.Operand) (*rlwe.Ciphertext, error)) CipherVector {
	if len(X) == 0 || len(Y) == 0 {
		return nil
	}
	size := broadcastSize(len(X), len(Y))
	res := make(CipherVector, size)
	if err := cryptoParams.WithEvaluator(func(eval *ckks.Evaluator) error {
		for i := 0; i < size; i++ {
			var err error
			res[i], err = op(eval, X[minIndex(i, len(X))], Y[minIndex(i, len(Y))])
			if err != nil {
				return err
			}
		}
		return nil
	}); err != nil {
		panic(err)
	}
	return res
}

func minIndex(i, length int) int {
	if length <= 1 {
		return 0
	}
	return i
}

func broadcastSize(xlen, ylen int) int {
	if xlen == ylen {
		return xlen
	}
	if xlen == 1 {
		return ylen
	}
	if ylen == 1 {
		return xlen
	}
	panic(fmt.Sprintf("incompatible vector lengths for broadcast: %d and %d", xlen, ylen))
}
