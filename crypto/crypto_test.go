package crypto

import (
	"fmt"
	"math"
	"path/filepath"
	"runtime"
	"strings"
	"testing"

	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/mat"
)

func newTestCryptoParams(t *testing.T) *CryptoParams {
	t.Helper()

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
	})
	if err != nil {
		t.Fatalf("failed to create test parameters: %v", err)
	}

	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	rlk := kgen.GenRelinearizationKeyNew(sk)

	cp := NewCryptoParams(params, sk, sk.CopyNew(), pk, rlk, 53, runtime.GOMAXPROCS(0))
	cp.SetRotKeys([]RotationType{
		{Value: 1, Side: SideLeft},
		{Value: 1, Side: SideRight},
	})
	return cp
}

func requireApproxSlice(t *testing.T, want, got []float64, tol float64) {
	t.Helper()
	if len(got) < len(want) {
		t.Fatalf("short slice: got %d, want at least %d", len(got), len(want))
	}
	for i := range want {
		if math.Abs(want[i]-got[i]) > tol {
			t.Fatalf("mismatch at %d: want=%f got=%f tol=%f", i, want[i], got[i], tol)
		}
	}
}

func requirePanicContains(t *testing.T, want string, fn func()) {
	t.Helper()

	defer func() {
		r := recover()
		if r == nil {
			t.Fatalf("expected panic containing %q", want)
		}
		msg := ""
		switch v := r.(type) {
		case error:
			msg = v.Error()
		default:
			msg = fmt.Sprint(v)
		}
		if !strings.Contains(msg, want) {
			t.Fatalf("panic mismatch: got %q want substring %q", msg, want)
		}
	}()

	fn()
}

// Crypto basics.

func TestEncryptDecryptFloatVector(t *testing.T) {
	cp := newTestCryptoParams(t)
	values := []float64{1.25, -2.5, 3.75, 4.5, -5.25}

	ct, n := EncryptFloatVector(cp, values)
	if n != len(values) {
		t.Fatalf("encrypted count mismatch: got %d want %d", n, len(values))
	}

	got := DecryptFloatVector(cp, ct, len(values))
	requireApproxSlice(t, values, got, 1e-6)
}

func TestBasicOps(t *testing.T) {
	cp := newTestCryptoParams(t)

	x, _ := EncryptFloatVector(cp, []float64{1, 2, 3, 4})
	y, _ := EncryptFloatVector(cp, []float64{5, 6, 7, 8})

	requireApproxSlice(t, []float64{6, 8, 10, 12}, DecryptFloatVector(cp, CAdd(cp, x, y), 4), 1e-6)
	requireApproxSlice(t, []float64{5, 12, 21, 32}, DecryptFloatVector(cp, CMult(cp, x, y), 4), 1e-4)
	requireApproxSlice(t, []float64{0.5, 1, 1.5, 2}, DecryptFloatVector(cp, CMultConst(cp, x, 0.5, false), 4), 1e-5)

	doubled := CMultConst(cp, x, 2.0, false)
	requireApproxSlice(t, []float64{2, 4, 6, 8}, DecryptFloatVector(cp, doubled, 4), 1e-6)
	if doubled[0].Level() != x[0].Level() {
		t.Fatalf("integer constant rescale consumed a level: got %d want %d", doubled[0].Level(), x[0].Level())
	}

	rotInput := make([]float64, cp.GetSlots())
	rotInput[0] = 1
	rotInput[1] = 2
	rotInput[2] = 3
	rotInput[cp.GetSlots()-1] = 4
	rotCt, _ := EncryptFloatVector(cp, rotInput)

	requireApproxSlice(t, []float64{4, 1, 2, 3}, DecryptMultipleFloat(cp, RotateRight(cp, rotCt[0], 1), 4), 1e-6)
}

// Small functions and invariants.

func TestCiphertextNegationPreservesScaleForPlainSub(t *testing.T) {
	cp := newTestCryptoParams(t)

	plain, _ := EncodeFloatVector(cp, []float64{10, 20, 30, 40})
	y, _ := EncryptFloatVector(cp, []float64{1.5, -2, 3.25, -4.5})

	got := CPSubOther(cp, plain, y)
	if got[0].Scale.Cmp(y[0].Scale) != 0 {
		t.Fatalf("negation changed scale: got %f want %f", got[0].Scale.Float64(), y[0].Scale.Float64())
	}
	requireApproxSlice(t, []float64{8.5, 22, 26.75, 44.5}, DecryptFloatVector(cp, got, 4), 1e-6)
}

func TestPlainSubCiphertextAtLowerLevel(t *testing.T) {
	cp := newTestCryptoParams(t)

	plain, _ := EncodeFloatVector(cp, []float64{10, 20, 30, 40})
	plainCt := EncryptPlaintextMatrix(cp, PlainMatrix{plain})[0]
	y, _ := EncryptFloatVector(cp, []float64{1.5, -2, 3.25, -4.5})
	lower := CMultConst(cp, y, 0.5, false)
	requireApproxSlice(t, []float64{0.75, -1, 1.625, -2.25}, DecryptFloatVector(cp, lower, 4), 1e-5)

	got := CPSubOther(cp, plain, lower)
	requireApproxSlice(t, []float64{9.25, 21, 28.375, 42.25}, DecryptFloatVector(cp, got, 4), 1e-5)

	plainCt = DropLevel(cp, CipherMatrix{plainCt}, lower[0].Level())[0]
	gotCt := CSub(cp, plainCt, lower)
	requireApproxSlice(t, []float64{9.25, 21, 28.375, 42.25}, DecryptFloatVector(cp, gotCt, 4), 1e-5)
}

func TestConstOpsDoNotMutateInputWhenNotInPlace(t *testing.T) {
	cp := newTestCryptoParams(t)

	x, _ := EncryptFloatVector(cp, []float64{1, 2, 3, 4})
	xBefore := DecryptFloatVector(cp, x, 4)

	scaled := CMultConst(cp, x, 0.5, false)
	requireApproxSlice(t, []float64{0.5, 1, 1.5, 2}, DecryptFloatVector(cp, scaled, 4), 1e-5)
	requireApproxSlice(t, xBefore, DecryptFloatVector(cp, x, 4), 1e-6)

	added := CAddConst(cp, x, 1.5)
	requireApproxSlice(t, []float64{2.5, 3.5, 4.5, 5.5}, DecryptFloatVector(cp, added, 4), 1e-5)
	requireApproxSlice(t, xBefore, DecryptFloatVector(cp, x, 4), 1e-6)
}

func TestMaskBoundsChecks(t *testing.T) {
	cp := newTestCryptoParams(t)
	cv, _ := EncryptFloatVector(cp, []float64{1, 2, 3, 4})

	requirePanicContains(t, "MaskTrunc:", func() {
		MaskTrunc(cp, cv[0], cp.GetSlots()+1)
	})
	requirePanicContains(t, "Mask:", func() {
		Mask(cp, cv[0], cp.GetSlots(), false)
	})
	requirePanicContains(t, "MaskWithScaling:", func() {
		MaskWithScaling(cp, cv[0], -1, false, 1.0)
	})
	requirePanicContains(t, "CMask:", func() {
		CMask(cp, cv, cp.GetSlots()*len(cv), false)
	})
}

func TestConcatCipherMatrixShapeChecks(t *testing.T) {
	cp := newTestCryptoParams(t)

	requirePanicContains(t, "ConcatCipherMatrix: matrix 1 has 2 rows, expected 1", func() {
		ConcatCipherMatrix([]CipherMatrix{
			{CZeros(cp, 1)},
			{CZeros(cp, 1), CZeros(cp, 1)},
		})
	})

	requirePanicContains(t, "ConcatCipherMatrix: matrix 0 row 1 has 2 columns, expected 1", func() {
		ConcatCipherMatrix([]CipherMatrix{
			{CZeros(cp, 1), CZeros(cp, 2)},
		})
	})
}

func TestFlattenLevelsAndDropLevel(t *testing.T) {
	cp := newTestCryptoParams(t)

	base, _ := EncryptFloatVector(cp, []float64{1, 2, 3, 4})
	scaled := CMultConst(cp, base, 0.5, false)

	uniform := CipherMatrix{
		CipherVector{base[0]},
		CipherVector{base[0].CopyNew()},
	}

	flatUniform, uniformLevel := FlattenLevels(cp, uniform)
	if uniformLevel != base[0].Level() {
		t.Fatalf("uniform level mismatch: got %d want %d", uniformLevel, base[0].Level())
	}
	if flatUniform[0][0] != uniform[0][0] {
		t.Fatalf("expected FlattenLevels to return original matrix when levels are already flat")
	}

	droppedSame := DropLevel(cp, uniform, base[0].Level())
	if droppedSame[0][0] == uniform[0][0] {
		t.Fatalf("expected DropLevel fast path to return copied ciphertexts")
	}
	requireApproxSlice(t, []float64{1, 2, 3, 4}, DecryptMultipleFloat(cp, droppedSame[0][0], 4), 1e-6)

	mixed := CipherMatrix{
		CipherVector{base[0]},
		CipherVector{scaled[0]},
	}
	flattenedMixed, mixedLevel := FlattenLevels(cp, mixed)
	if mixedLevel != scaled[0].Level() {
		t.Fatalf("mixed level mismatch: got %d want %d", mixedLevel, scaled[0].Level())
	}
	for i := range flattenedMixed {
		for j := range flattenedMixed[i] {
			if flattenedMixed[i][j].Level() != mixedLevel {
				t.Fatalf("flattened level mismatch at (%d,%d): got %d want %d", i, j, flattenedMixed[i][j].Level(), mixedLevel)
			}
		}
	}

	requirePanicContains(t, "DropLevel:", func() {
		DropLevel(cp, mixed, base[0].Level()+1)
	})
}

// Round trips.

func TestDenseRoundTrip(t *testing.T) {
	cp := newTestCryptoParams(t)
	input := mat.NewDense(3, 2, []float64{
		1, 2,
		3, 4,
		5, 6,
	})
	decoded := PlaintextToDense(cp, EncodeDense(cp, input), 3)
	if !mat.EqualApprox(input, decoded, 1e-6) {
		t.Fatalf("dense mismatch:\nwant:\n%v\ngot:\n%v", mat.Formatted(input), mat.Formatted(decoded))
	}
}

func TestCipherMatrixMarshalRoundTrip(t *testing.T) {
	cp := newTestCryptoParams(t)
	input := [][]float64{
		{1, 2, 3},
		{4, 5, 6},
	}

	cm, _, _, err := EncryptFloatMatrixRow(cp, input)
	if err != nil {
		t.Fatalf("encrypt matrix: %v", err)
	}

	rawSizes, rawData := MarshalCM(cm)
	got := UnmarshalCM(cp, len(cm), len(cm[0]), rawSizes, rawData)
	gotMatrix := DecryptFloatMatrix(cp, got, 3)

	for i := range input {
		requireApproxSlice(t, input[i], gotMatrix[i], 1e-6)
	}
}

func TestCipherMatrixFileRoundTrip(t *testing.T) {
	cp := newTestCryptoParams(t)
	input := [][]float64{
		{1.5, 2.5},
		{3.5, 4.5},
	}
	cm, _, _, err := EncryptFloatMatrixRow(cp, input)
	if err != nil {
		t.Fatalf("encrypt matrix: %v", err)
	}

	dir := t.TempDir()
	path := filepath.Join(dir, "cm.bin")
	SaveCipherMatrixToFile(cm, path)

	got := LoadCipherMatrixFromFile(cp, path)
	gotMatrix := DecryptFloatMatrix(cp, got, 2)
	for i := range input {
		requireApproxSlice(t, input[i], gotMatrix[i], 1e-6)
	}
}

func TestCryptoParamsMarshalRoundTrip(t *testing.T) {
	cp := newTestCryptoParams(t)
	data, err := cp.MarshalBinary()
	if err != nil {
		t.Fatalf("marshal crypto params: %v", err)
	}

	var roundTrip CryptoParams
	if err := roundTrip.UnmarshalBinary(data); err != nil {
		t.Fatalf("unmarshal crypto params: %v", err)
	}

	if !cp.Sk.Equal(roundTrip.Sk) {
		t.Fatalf("secret key mismatch after round-trip")
	}
	if !cp.AggregateSk.Equal(roundTrip.AggregateSk) {
		t.Fatalf("aggregate secret key mismatch after round-trip")
	}
	if !cp.Pk.Equal(roundTrip.Pk) {
		t.Fatalf("public key mismatch after round-trip")
	}

	want, _ := EncryptFloatVector(&roundTrip, []float64{9, 8, 7})
	got := DecryptFloatVector(&roundTrip, want, 3)
	requireApproxSlice(t, []float64{9, 8, 7}, got, 1e-6)
}
