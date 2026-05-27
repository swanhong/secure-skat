package mpc_test

import (
	"math"
	"testing"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/multiparty"
	"github.com/tuneinsight/lattigo/v6/multiparty/mpckks"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

// Educational tests for the exact Lattigo multiparty functions that this repo
// uses. Each test is meant to be runnable on its own:
//
//   go test ./mpc -run TestLattigoMultiparty_PublicKeyGen_Example -v
//   go test ./mpc -run TestLattigoMultiparty_RelinearizationKeyGen_Example -v
//   go test ./mpc -run TestLattigoMultiparty_GaloisKeyGen_Example -v
//   go test ./mpc -run TestLattigoMultiparty_KeySwitchToZero_Example -v
//   go test ./mpc -run TestLattigoMultiparty_Refresh_Example -v

type lattigoMultipartyExampleContext struct {
	params   ckks.Parameters
	parties  int
	kgen     *rlwe.KeyGenerator
	encoder  *ckks.Encoder
	skShares []*rlwe.SecretKey
	skIdeal  *rlwe.SecretKey
	crs      sampling.PRNG
}

func newLattigoMultipartyExampleContext(t *testing.T) *lattigoMultipartyExampleContext {
	t.Helper()

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            13,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{61},
		LogDefaultScale: 45,
	})
	if err != nil {
		t.Fatalf("new parameters: %v", err)
	}

	const parties = 3

	kgen := ckks.NewKeyGenerator(params)
	skShares := make([]*rlwe.SecretKey, parties)
	skIdeal := rlwe.NewSecretKey(params)
	for i := range skShares {
		skShares[i] = kgen.GenSecretKeyNew()
		params.RingQP().Add(skIdeal.Value, skShares[i].Value, skIdeal.Value)
	}

	crs, err := sampling.NewKeyedPRNG([]byte("lattigo-multiparty-example-seed-000000000000"))
	if err != nil {
		t.Fatalf("new keyed prng: %v", err)
	}

	return &lattigoMultipartyExampleContext{
		params:   params,
		parties:  parties,
		kgen:     kgen,
		encoder:  ckks.NewEncoder(params),
		skShares: skShares,
		skIdeal:  skIdeal,
		crs:      crs,
	}
}

func encodeFloatValues(t *testing.T, tc *lattigoMultipartyExampleContext, values []float64, level int) *rlwe.Plaintext {
	t.Helper()

	pt := ckks.NewPlaintext(tc.params, level)
	pt.LogDimensions.Cols = 2 // 4 slots for easy mental tracing.
	if err := tc.encoder.Encode(values, pt); err != nil {
		t.Fatalf("encode values: %v", err)
	}
	return pt
}

func encryptFloatValues(t *testing.T, tc *lattigoMultipartyExampleContext, enc *rlwe.Encryptor, values []float64) *rlwe.Ciphertext {
	t.Helper()

	pt := encodeFloatValues(t, tc, values, tc.params.MaxLevel())
	ct, err := enc.EncryptNew(pt)
	if err != nil {
		t.Fatalf("encrypt values: %v", err)
	}
	return ct
}

func decryptFloatValues(t *testing.T, tc *lattigoMultipartyExampleContext, dec *rlwe.Decryptor, ct *rlwe.Ciphertext, wantLen int) []float64 {
	t.Helper()

	pt := ckks.NewPlaintext(tc.params, ct.Level())
	dec.Decrypt(ct, pt)

	got := make([]float64, 1<<pt.LogDimensions.Cols)
	if err := tc.encoder.Decode(pt, got); err != nil {
		t.Fatalf("decode plaintext: %v", err)
	}

	return append([]float64(nil), got[:wantLen]...)
}

func decodeZeroKeyCiphertext(t *testing.T, tc *lattigoMultipartyExampleContext, ct *rlwe.Ciphertext, wantLen int) []float64 {
	t.Helper()

	pt := ct.Plaintext()
	got := make([]float64, 1<<pt.LogDimensions.Cols)
	if err := tc.encoder.Decode(pt, got); err != nil {
		t.Fatalf("decode zero-key plaintext: %v", err)
	}

	return append([]float64(nil), got[:wantLen]...)
}

func requireApproxSlice(t *testing.T, want, got []float64, tol float64) {
	t.Helper()

	if len(want) != len(got) {
		t.Fatalf("length mismatch: want %d, got %d", len(want), len(got))
	}

	for i := range want {
		if math.Abs(want[i]-got[i]) > tol {
			t.Fatalf("slot %d mismatch: want %.6f, got %.6f, tol %.6f", i, want[i], got[i], tol)
		}
	}
}

func TestLattigoMultiparty_PublicKeyGen_Example(t *testing.T) {
	tc := newLattigoMultipartyExampleContext(t)
	input := []float64{0.25, -1.5, 2.75, -0.125}

	ckg := make([]multiparty.PublicKeyGenProtocol, tc.parties)
	for i := range ckg {
		if i == 0 {
			ckg[i] = multiparty.NewPublicKeyGenProtocol(tc.params)
		} else {
			ckg[i] = ckg[0].ShallowCopy()
		}
	}

	shares := make([]multiparty.PublicKeyGenShare, tc.parties)
	for i := range shares {
		shares[i] = ckg[i].AllocateShare()
	}

	crp := ckg[0].SampleCRP(tc.crs)
	for i := range shares {
		ckg[i].GenShare(tc.skShares[i], crp, &shares[i])
	}

	for i := 1; i < tc.parties; i++ {
		ckg[0].AggregateShares(shares[0], shares[i], &shares[0])
	}

	shareBytes, err := shares[0].MarshalBinary()
	if err != nil {
		t.Fatalf("marshal public-key share: %v", err)
	}

	pk := rlwe.NewPublicKey(tc.params)
	ckg[0].GenPublicKey(shares[0], crp, pk)

	encryptor := ckks.NewEncryptor(tc.params, pk)
	decryptor := ckks.NewDecryptor(tc.params, tc.skIdeal)
	ct := encryptFloatValues(t, tc, encryptor, input)
	got := decryptFloatValues(t, tc, decryptor, ct, len(input))

	t.Logf("protocol: AllocateShare -> SampleCRP -> GenShare -> AggregateShares -> GenPublicKey")
	t.Logf("input values:  %v", input)
	t.Logf("output values: %v", got)
	t.Logf("aggregate share bytes: %d", len(shareBytes))

	requireApproxSlice(t, input, got, 1e-3)
}

func TestLattigoMultiparty_RelinearizationKeyGen_Example(t *testing.T) {
	tc := newLattigoMultipartyExampleContext(t)
	left := []float64{1.0, 2.0, 3.0, 4.0}
	right := []float64{0.5, -1.0, 1.5, -2.0}
	want := []float64{0.5, -2.0, 4.5, -8.0}

	// First create the collective public key used to encrypt example inputs.
	ckg := multiparty.NewPublicKeyGenProtocol(tc.params)
	crpPk := ckg.SampleCRP(tc.crs)
	pkShare := make([]multiparty.PublicKeyGenShare, tc.parties)
	for i := range pkShare {
		pkShare[i] = ckg.AllocateShare()
		ckg.GenShare(tc.skShares[i], crpPk, &pkShare[i])
		if i > 0 {
			ckg.AggregateShares(pkShare[0], pkShare[i], &pkShare[0])
		}
	}
	pk := rlwe.NewPublicKey(tc.params)
	ckg.GenPublicKey(pkShare[0], crpPk, pk)

	rkg := make([]multiparty.RelinearizationKeyGenProtocol, tc.parties)
	for i := range rkg {
		if i == 0 {
			rkg[i] = multiparty.NewRelinearizationKeyGenProtocol(tc.params)
		} else {
			rkg[i] = rkg[0].ShallowCopy()
		}
	}

	ephSk := make([]*rlwe.SecretKey, tc.parties)
	share1 := make([]multiparty.RelinearizationKeyGenShare, tc.parties)
	share2 := make([]multiparty.RelinearizationKeyGenShare, tc.parties)
	for i := range rkg {
		ephSk[i], share1[i], share2[i] = rkg[i].AllocateShare()
	}

	crpRlk := rkg[0].SampleCRP(tc.crs)
	for i := range rkg {
		rkg[i].GenShareRoundOne(tc.skShares[i], crpRlk, ephSk[i], &share1[i])
	}
	for i := 1; i < tc.parties; i++ {
		rkg[0].AggregateShares(share1[0], share1[i], &share1[0])
	}

	for i := range rkg {
		rkg[i].GenShareRoundTwo(ephSk[i], tc.skShares[i], share1[0], &share2[i])
	}
	for i := 1; i < tc.parties; i++ {
		rkg[0].AggregateShares(share2[0], share2[i], &share2[0])
	}

	rlk := rlwe.NewRelinearizationKey(tc.params)
	rkg[0].GenRelinearizationKey(share1[0], share2[0], rlk)

	encryptor := ckks.NewEncryptor(tc.params, pk)
	decryptor := ckks.NewDecryptor(tc.params, tc.skIdeal)
	evaluator := ckks.NewEvaluator(tc.params, rlwe.NewMemEvaluationKeySet(rlk))

	ctLeft := encryptFloatValues(t, tc, encryptor, left)
	ctRight := encryptFloatValues(t, tc, encryptor, right)
	ctMul, err := evaluator.MulRelinNew(ctLeft, ctRight)
	if err != nil {
		t.Fatalf("mul-relin: %v", err)
	}
	if err := evaluator.Rescale(ctMul, ctMul); err != nil {
		t.Fatalf("rescale: %v", err)
	}

	got := decryptFloatValues(t, tc, decryptor, ctMul, len(want))

	t.Logf("protocol: AllocateShare -> SampleCRP -> GenShareRoundOne -> AggregateShares -> GenShareRoundTwo -> AggregateShares -> GenRelinearizationKey")
	t.Logf("left input:   %v", left)
	t.Logf("right input:  %v", right)
	t.Logf("want product: %v", want)
	t.Logf("got product:  %v", got)

	requireApproxSlice(t, want, got, 2e-2)
}

func TestLattigoMultiparty_GaloisKeyGen_Example(t *testing.T) {
	tc := newLattigoMultipartyExampleContext(t)
	input := []float64{10, 20, 30, 40}
	want := []float64{20, 30, 40, 10}

	ckg := multiparty.NewPublicKeyGenProtocol(tc.params)
	crpPk := ckg.SampleCRP(tc.crs)
	pkShare := make([]multiparty.PublicKeyGenShare, tc.parties)
	for i := range pkShare {
		pkShare[i] = ckg.AllocateShare()
		ckg.GenShare(tc.skShares[i], crpPk, &pkShare[i])
		if i > 0 {
			ckg.AggregateShares(pkShare[0], pkShare[i], &pkShare[0])
		}
	}
	pk := rlwe.NewPublicKey(tc.params)
	ckg.GenPublicKey(pkShare[0], crpPk, pk)

	gkg := make([]multiparty.GaloisKeyGenProtocol, tc.parties)
	for i := range gkg {
		if i == 0 {
			gkg[i] = multiparty.NewGaloisKeyGenProtocol(tc.params)
		} else {
			gkg[i] = gkg[0].ShallowCopy()
		}
	}

	galEl := tc.params.GaloisElementForRotation(1)
	shares := make([]multiparty.GaloisKeyGenShare, tc.parties)
	for i := range shares {
		shares[i] = gkg[i].AllocateShare()
	}

	crpGk := gkg[0].SampleCRP(tc.crs)
	for i := range shares {
		if err := gkg[i].GenShare(tc.skShares[i], galEl, crpGk, &shares[i]); err != nil {
			t.Fatalf("gen galois share: %v", err)
		}
	}
	for i := 1; i < tc.parties; i++ {
		if err := gkg[0].AggregateShares(shares[0], shares[i], &shares[0]); err != nil {
			t.Fatalf("aggregate galois share: %v", err)
		}
	}

	gk := rlwe.NewGaloisKey(tc.params)
	if err := gkg[0].GenGaloisKey(shares[0], crpGk, gk); err != nil {
		t.Fatalf("gen galois key: %v", err)
	}

	encryptor := ckks.NewEncryptor(tc.params, pk)
	decryptor := ckks.NewDecryptor(tc.params, tc.skIdeal)
	evaluator := ckks.NewEvaluator(tc.params, rlwe.NewMemEvaluationKeySet(nil, gk))

	ct := encryptFloatValues(t, tc, encryptor, input)
	rotated, err := evaluator.RotateNew(ct, 1)
	if err != nil {
		t.Fatalf("rotate: %v", err)
	}
	got := decryptFloatValues(t, tc, decryptor, rotated, len(want))

	t.Logf("protocol: AllocateShare -> SampleCRP -> GenShare -> AggregateShares -> GenGaloisKey")
	t.Logf("rotation: left by 1 slot")
	t.Logf("input values:  %v", input)
	t.Logf("output values: %v", got)

	requireApproxSlice(t, want, got, 1e-3)
}

func TestLattigoMultiparty_KeySwitchToZero_Example(t *testing.T) {
	tc := newLattigoMultipartyExampleContext(t)
	input := []float64{1.25, -2.5, 3.75, -4.5}

	ckg := multiparty.NewPublicKeyGenProtocol(tc.params)
	crpPk := ckg.SampleCRP(tc.crs)
	pkShare := make([]multiparty.PublicKeyGenShare, tc.parties)
	for i := range pkShare {
		pkShare[i] = ckg.AllocateShare()
		ckg.GenShare(tc.skShares[i], crpPk, &pkShare[i])
		if i > 0 {
			ckg.AggregateShares(pkShare[0], pkShare[i], &pkShare[0])
		}
	}
	pk := rlwe.NewPublicKey(tc.params)
	ckg.GenPublicKey(pkShare[0], crpPk, pk)

	ksp := make([]multiparty.KeySwitchProtocol, tc.parties)
	for i := range ksp {
		prot, err := multiparty.NewKeySwitchProtocol(tc.params, tc.params.Xe())
		if err != nil {
			t.Fatalf("new keyswitch protocol: %v", err)
		}
		if i == 0 {
			ksp[i] = prot
		} else {
			ksp[i] = ksp[0].ShallowCopy()
		}
	}

	encryptor := ckks.NewEncryptor(tc.params, pk)
	ct := encryptFloatValues(t, tc, encryptor, input)

	zero := rlwe.NewSecretKey(tc.params)
	shares := make([]multiparty.KeySwitchShare, tc.parties)
	for i := range shares {
		shares[i] = ksp[i].AllocateShare(ct.Level())
		ksp[i].GenShare(tc.skShares[i], zero, ct, &shares[i])
		if i > 0 {
			if err := ksp[0].AggregateShares(shares[0], shares[i], &shares[0]); err != nil {
				t.Fatalf("aggregate keyswitch share: %v", err)
			}
		}
	}

	zeroKeyCt := rlwe.NewCiphertext(tc.params, 1, ct.Level())
	ksp[0].KeySwitch(ct, shares[0], zeroKeyCt)
	got := decodeZeroKeyCiphertext(t, tc, zeroKeyCt, len(input))

	t.Logf("protocol: AllocateShare -> GenShare(sk_i -> 0) -> AggregateShares -> KeySwitch")
	t.Logf("input values:  %v", input)
	t.Logf("output values: %v", got)
	t.Logf("note: this is the exact pattern our code uses for collective decrypt")

	requireApproxSlice(t, input, got, 1e-3)
}

func TestLattigoMultiparty_Refresh_Example(t *testing.T) {
	tc := newLattigoMultipartyExampleContext(t)
	input := []float64{0.125, -0.75, 1.5, 2.25}

	ckg := multiparty.NewPublicKeyGenProtocol(tc.params)
	crpPk := ckg.SampleCRP(tc.crs)
	pkShare := make([]multiparty.PublicKeyGenShare, tc.parties)
	for i := range pkShare {
		pkShare[i] = ckg.AllocateShare()
		ckg.GenShare(tc.skShares[i], crpPk, &pkShare[i])
		if i > 0 {
			ckg.AggregateShares(pkShare[0], pkShare[i], &pkShare[0])
		}
	}
	pk := rlwe.NewPublicKey(tc.params)
	ckg.GenPublicKey(pkShare[0], crpPk, pk)

	minLevel, logBound, ok := mpckks.GetMinimumLevelForRefresh(128, tc.params.DefaultScale(), tc.parties, tc.params.Q())
	if !ok || minLevel+1 > tc.params.MaxLevel() {
		t.Fatalf("refresh parameters unsupported: ok=%v minLevel=%d maxLevel=%d", ok, minLevel, tc.params.MaxLevel())
	}

	encryptor := ckks.NewEncryptor(tc.params, pk)
	decryptor := ckks.NewDecryptor(tc.params, tc.skIdeal)
	evaluator := ckks.NewEvaluator(tc.params, nil)

	ct := encryptFloatValues(t, tc, encryptor, input)
	levelBeforeDrop := ct.Level()
	evaluator.DropLevel(ct, ct.Level()-(minLevel+1))
	levelBeforeRefresh := ct.Level()

	rfp := make([]mpckks.RefreshProtocol, tc.parties)
	for i := range rfp {
		prot, err := mpckks.NewRefreshProtocol(tc.params, 0, tc.params.Xe())
		if err != nil {
			t.Fatalf("new refresh protocol: %v", err)
		}
		if i == 0 {
			rfp[i] = prot
		} else {
			rfp[i] = rfp[0].ShallowCopy()
		}
	}

	crpRefresh := rfp[0].SampleCRP(tc.params.MaxLevel(), tc.crs)
	shares := make([]multiparty.RefreshShare, tc.parties)
	for i := range shares {
		shares[i] = rfp[i].AllocateShare(ct.Level(), tc.params.MaxLevel())
		if err := rfp[i].GenShare(tc.skShares[i], logBound, ct, crpRefresh, &shares[i]); err != nil {
			t.Fatalf("gen refresh share: %v", err)
		}
		if i > 0 {
			if err := rfp[0].AggregateShares(&shares[0], &shares[i], &shares[0]); err != nil {
				t.Fatalf("aggregate refresh share: %v", err)
			}
		}
	}

	refreshed := ckks.NewCiphertext(tc.params, 1, tc.params.MaxLevel())
	if err := rfp[0].Finalize(ct, crpRefresh, shares[0], refreshed); err != nil {
		t.Fatalf("finalize refresh: %v", err)
	}

	got := decryptFloatValues(t, tc, decryptor, refreshed, len(input))

	t.Logf("protocol: GetMinimumLevelForRefresh -> AllocateShare -> SampleCRP -> GenShare -> AggregateShares -> Finalize")
	t.Logf("input values:         %v", input)
	t.Logf("output values:        %v", got)
	t.Logf("ciphertext level:     fresh=%d dropped=%d refreshed=%d", levelBeforeDrop, levelBeforeRefresh, refreshed.Level())
	t.Logf("default scale reset?: %v", refreshed.Scale.Float64() == tc.params.DefaultScale().Float64())

	requireApproxSlice(t, input, got, 2e-2)
}
