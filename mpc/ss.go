package mpc

import (
	crand "crypto/rand"
	"math/big"

	mpc_core "github.com/hhcho/mpc-core"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

func plaintextFromPoly(cryptoParams *crypto.CryptoParams, poly ring.Poly, src *rlwe.Ciphertext, level int) *rlwe.Plaintext {
	pt := ckks.NewPlaintext(cryptoParams.Params, level)
	pt.Value = poly
	pt.Element.Value[0] = poly
	pt.Scale = src.Scale
	pt.IsNTT = src.IsNTT
	return pt
}

func (mpcObj *MPC) SSMultMat(a, b mpc_core.RMat) mpc_core.RMat {
	ar, am := mpcObj.BeaverPartitionMat(a)
	br, bm := mpcObj.BeaverPartitionMat(b)
	ab := mpcObj.BeaverMultMat(ar, am, br, bm)
	return mpcObj.BeaverReconstructMat(ab)
}

func (mpcObj *MPC) SSMultElemVecScalar(a mpc_core.RVec, b mpc_core.RElem) mpc_core.RVec {
	ar, am := mpcObj.BeaverPartitionVec(a)
	br, bm := mpcObj.BeaverPartition(b)
	x := mpc_core.InitRVec(mpcObj.rtype.Zero(), len(a))
	for i := range x {
		x[i] = mpcObj.BeaverMult(ar[i], am[i], br, bm)
	}
	return mpcObj.BeaverReconstructVec(x)
}

func (mpcObj *MPC) SSSquareElemVec(a mpc_core.RVec) mpc_core.RVec {
	ar, am := mpcObj.BeaverPartitionVec(a)
	x := mpcObj.BeaverMultElemVec(ar, am, ar, am)
	return mpcObj.BeaverReconstructVec(x)
}

func (mpcObj *MPC) SSMultElemVec(a, b mpc_core.RVec) mpc_core.RVec {
	ar, am := mpcObj.BeaverPartitionVec(a)
	br, bm := mpcObj.BeaverPartitionVec(b)
	x := mpcObj.BeaverMultElemVec(ar, am, br, bm)
	return mpcObj.BeaverReconstructVec(x)
}

func (mpcObj *MPC) SSMultElemMat(a, b mpc_core.RMat) mpc_core.RMat {
	ar, am := mpcObj.BeaverPartitionMat(a)
	br, bm := mpcObj.BeaverPartitionMat(b)
	x := mpcObj.BeaverMultElemMat(ar, am, br, bm)
	return mpcObj.BeaverReconstructMat(x)
}

// --- Secure share <-> ciphertext conversion (slot-domain, masked) ---
//
// Re-port of the v2 masked conversions (which used the patched Lattigo-v2 encoder methods
// EncodeRVecNew/DecodeRVec) to stock Lattigo v6, via full-precision []*big.Float slot
// encoding so field-sized masks survive. No party ever observes a cleartext value.

// elemToSignedBigFloat returns the centered (signed) fixed-point real value of e.
// Only the field types used by this codebase (LElem128/LElem256) are supported.
func elemToSignedBigFloat(e mpc_core.RElem, fracBits int) *big.Float {
	switch v := e.(type) {
	case mpc_core.LElem256:
		return v.ToSignedBigFloat(fracBits)
	case mpc_core.LElem128:
		return v.ToSignedBigFloat(fracBits)
	default:
		panic("mpc: SS<->cipher conversion only supports LElem128 or LElem256")
	}
}

// signedBigFloatToElem converts a signed fixed-point real value back into a field
// element: round(f * 2^fracBits) reduced modulo the field modulus. The explicit
// reduction is required because RElem.FromBigInt does not reduce on construction.
func signedBigFloatToElem(rtype mpc_core.RElem, f *big.Float, fracBits int) mpc_core.RElem {
	prec := f.Prec()
	if prec < 256 {
		prec = 256
	}
	scaled := new(big.Float).SetPrec(prec).Mul(f, bigFloatPow2(fracBits, prec))
	half := new(big.Float).SetPrec(prec).SetFloat64(0.5)
	if scaled.Sign() >= 0 {
		scaled.Add(scaled, half)
	} else {
		scaled.Sub(scaled, half)
	}
	r, _ := scaled.Int(nil)
	r.Mod(r, rtype.Modulus())
	return rtype.FromBigInt(r)
}

func bigFloatPow2(n int, prec uint) *big.Float {
	return new(big.Float).SetPrec(prec).SetInt(new(big.Int).Lsh(big.NewInt(1), uint(n)))
}

// sampleCenteredMask draws a per-element uniform mask in (-bound/2, bound/2], stored
// in the field representation. Each party samples its own mask locally.
func sampleCenteredMask(rtype mpc_core.RElem, nrow, ncol int, bound *big.Int) mpc_core.RMat {
	modulus := rtype.Modulus()
	boundHalf := new(big.Int).Rsh(bound, 1)
	mask := make(mpc_core.RMat, nrow)
	for i := range mask {
		mask[i] = make(mpc_core.RVec, ncol)
		for j := range mask[i] {
			t, err := crand.Int(crand.Reader, bound)
			if err != nil {
				panic(err)
			}
			if t.Cmp(boundHalf) >= 0 {
				t.Sub(t, bound)
			}
			t.Mod(t, modulus)
			mask[i][j] = rtype.FromBigInt(t)
		}
	}
	return mask
}

// encodeElemSlotsToCM encrypts row-major slot-domain field shares into a CipherMatrix
// under the collective public key, using full-precision big.Float encoding so that
// field-sized (masked) values are represented exactly. scale/level set the encoding of
// the resulting ciphertexts.
func (mpcObj *MPC) encodeElemSlotsToCM(cryptoParams *crypto.CryptoParams, share mpc_core.RMat, scale rlwe.Scale, level int) crypto.CipherMatrix {
	slots := cryptoParams.GetSlots()
	fracBits := mpcObj.GetFracBits()
	numCtxRow := len(share)
	nElemCol := len(share[0])
	numCtxCol := 1 + ((nElemCol - 1) / slots)

	pm := make(crypto.PlainMatrix, numCtxRow)
	for i := range pm {
		pm[i] = make(crypto.PlainVector, numCtxCol)
		for j := 0; j < numCtxCol; j++ {
			start := j * slots
			end := start + slots
			if end > nElemCol {
				end = nElemCol
			}
			buf := make([]*big.Float, slots)
			for k := range buf {
				buf[k] = new(big.Float)
			}
			for k := start; k < end; k++ {
				buf[k-start] = elemToSignedBigFloat(share[i][k], fracBits)
			}
			pt := ckks.NewPlaintext(cryptoParams.Params, level)
			pt.Scale = scale
			if err := cryptoParams.WithEncoder(func(encoder *ckks.Encoder) error {
				return encoder.Encode(buf, pt)
			}); err != nil {
				panic(err)
			}
			pm[i][j] = pt
		}
	}

	return crypto.EncryptPlaintextMatrix(cryptoParams, pm)
}

// decodePlainToElemSlots decodes plaintexts back to slot-domain field elements
// (nElemRow per row), reversing encodeElemSlotsToCM.
func (mpcObj *MPC) decodePlainToElemSlots(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, pm crypto.PlainMatrix, nElemRow int) mpc_core.RMat {
	slots := cryptoParams.GetSlots()
	fracBits := mpcObj.GetFracBits()
	numCtxRow := len(pm)

	// The masked plaintext carries a full-Q coefficient mask (~log2(Q) bits), so the
	// slot decode must keep far more precision than the pooled encoder (config prec):
	// otherwise the huge mask rounds away the signal and the field shares don't cancel.
	// v2's DecodeRVec decoded exact big.Int coefficients through a big-float FFT; we match
	// that by decoding with a high-precision encoder.
	decPrec := uint(cryptoParams.Params.LogQP() + 128)
	encoder := ckks.NewEncoder(cryptoParams.Params, decPrec)

	out := mpc_core.InitRMat(rtype.Zero(), numCtxRow, nElemRow)
	for i := range pm {
		for j := range pm[i] {
			buf := make([]*big.Float, slots)
			for k := range buf {
				buf[k] = new(big.Float).SetPrec(decPrec)
			}
			if err := encoder.Decode(pm[i][j], buf); err != nil {
				panic(err)
			}
			start := j * slots
			for k := 0; k < slots && start+k < nElemRow; k++ {
				out[i][start+k] = signedBigFloatToElem(rtype, buf[k], fracBits)
			}
		}
	}
	return out
}

// SSToCMat converts secret-shared values (row-major, slot-domain) into a CipherMatrix
// encrypting the reconstructed value. Faithful to the v2 masked-share protocol: each
// party masks its share, only x - Sum(mask_i) is revealed, and per-party encryptions of
// the masked share are summed to Enc(x). No party observes x.
func (mpcObj *MPC) SSToCMat(cryptoParams *crypto.CryptoParams, rm mpc_core.RMat) (cm crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		cm = make(crypto.CipherMatrix, 1)
		cm[0] = make(crypto.CipherVector, 1)
		return
	}

	rtype := rm.Type().Zero()
	if rtype.TypeID() != mpc_core.LElem256UniqueID && rtype.TypeID() != mpc_core.LElem128UniqueID {
		panic("SSToCMat only supported for LElem128 or LElem256")
	}

	numCtxRow := len(rm)
	nElemCol := len(rm[0])

	// Bound masks so that x - Sum(mask_i) and the hub share never wrap the field
	// (matches v2: modulus / (4 * #data-parties)).
	bound := new(big.Int).Set(rtype.Modulus())
	bound.Quo(bound, big.NewInt(int64(4*(mpcObj.GetNParty()-1))))

	mask := sampleCenteredMask(rtype, numCtxRow, nElemCol, bound)

	rmMask := rm.Copy()
	rmMask.Sub(mask)
	rmMask = mpcObj.RevealSymMat(rmMask) // reveals only x - Sum(mask_i)

	var share mpc_core.RMat
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		share = rmMask
		share.Add(mask)
	} else {
		share = mask
	}

	cm = mpcObj.encodeElemSlotsToCM(cryptoParams, share, cryptoParams.Params.DefaultScale(), cryptoParams.Params.MaxLevel())
	return mpcObj.Network.AggregateCMat(cryptoParams, cm)
}

func (mpcObj *MPC) SSToCVec(cryptoParams *crypto.CryptoParams, rv mpc_core.RVec) (cv crypto.CipherVector) {
	return mpcObj.SSToCMat(cryptoParams, mpc_core.RMat{rv})[0]
}

func (mpcObj *MPC) SStoCiphertext(cryptoParams *crypto.CryptoParams, rv mpc_core.RVec) *rlwe.Ciphertext {
	return mpcObj.SSToCVec(cryptoParams, rv)[0]
}

func (mpcObj *MPC) CMatToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, cm crypto.CipherMatrix, sourcePid, numCtxRow, numCtxCol, nElemRow int) (rm mpc_core.RMat) {
	rm = mpc_core.InitRMat(rtype.Zero(), numCtxRow, nElemRow)
	if mpcObj.GetPid() == 0 {
		return
	}
	if rtype.TypeID() != mpc_core.LElem256UniqueID && rtype.TypeID() != mpc_core.LElem128UniqueID {
		panic("CMatToSS only supported for LElem128 or LElem256")
	}

	if sourcePid > 0 {
		cm = mpcObj.Network.BroadcastCMat(cryptoParams, cm, sourcePid, numCtxRow, numCtxCol)
	}

	cm, level := crypto.FlattenLevels(cryptoParams, cm)

	ringQ := cryptoParams.Params.RingQ().AtLevel(level)
	nCoeffs := cryptoParams.Params.N()
	nParty := mpcObj.GetNParty()

	// v2 mask bound = (product of active Q moduli) / (2 * #data-parties).
	bound := new(big.Int).Set(cryptoParams.Params.RingQ().ModulusAtLevel[level])
	bound.Quo(bound, big.NewInt(int64(2*(nParty-1))))
	if bound.Sign() <= 0 {
		panic("CMatToSS: ciphertext level too low to mask securely")
	}
	boundHalf := new(big.Int).Rsh(bound, 1)

	// Private smudge sampler (independent PRNG, never the shared CRS).
	prng, err := sampling.NewPRNG()
	if err != nil {
		panic(err)
	}
	gauss := ring.NewGaussianSampler(prng, cryptoParams.Params.RingQ(), ring.DiscreteGaussian{Sigma: 3.19, Bound: 19}, false).AtLevel(level)

	// Per ciphertext, build the masked partial-decryption share h0 = sk*c1 + mask + e
	// (coefficient-domain mask, NTT'd), keeping each mask poly locally for the decode.
	maskPolys := make([][]ring.Poly, numCtxRow)
	shares := make([][]ring.Poly, numCtxRow)
	for i := range cm {
		maskPolys[i] = make([]ring.Poly, len(cm[i]))
		shares[i] = make([]ring.Poly, len(cm[i]))
		for j := range cm[i] {
			maskBig := make([]*big.Int, nCoeffs)
			for t := range maskBig {
				m, cerr := crand.Int(crand.Reader, bound)
				if cerr != nil {
					panic(cerr)
				}
				if m.Cmp(boundHalf) >= 0 {
					m.Sub(m, bound)
				}
				maskBig[t] = m
			}
			maskPoly := ringQ.NewPoly()
			ringQ.SetCoefficientsBigint(maskBig, maskPoly)
			ringQ.NTT(maskPoly, maskPoly)
			maskPolys[i][j] = maskPoly

			h0 := *maskPoly.CopyNew()
			ringQ.MulCoeffsMontgomeryThenAdd(cryptoParams.Sk.Value.Q, cm[i][j].Value[1], h0)
			eNoise := gauss.ReadNew()
			ringQ.NTT(eNoise, eNoise)
			ringQ.Add(h0, eNoise, h0)
			shares[i][j] = h0
		}
	}

	// Aggregate the decryption shares across parties: agg = sk*c1 + Sum mask + Sum e.
	agg := make([][]ring.Poly, numCtxRow)
	for i := range shares {
		agg[i] = mpcObj.Network.AggregateDecryptShareVec(cryptoParams, shares[i], level)
	}

	// c0 + agg encodes the masked plaintext (x + Sum mask); decode it and the local mask
	// to the field (linear big.Float slot decode = v2's patched DecodeRVec).
	slots := cryptoParams.GetSlots()
	maskWidth := len(cm[0]) * slots
	pmSum := make(crypto.PlainMatrix, numCtxRow)
	pmMask := make(crypto.PlainMatrix, numCtxRow)
	for i := range cm {
		pmSum[i] = make(crypto.PlainVector, len(cm[i]))
		pmMask[i] = make(crypto.PlainVector, len(cm[i]))
		for j := range cm[i] {
			sumPoly := *cm[i][j].Value[0].CopyNew()
			ringQ.Add(sumPoly, agg[i][j], sumPoly)
			pmSum[i][j] = plaintextFromPoly(cryptoParams, sumPoly, cm[i][j], level)
			pmMask[i][j] = plaintextFromPoly(cryptoParams, maskPolys[i][j], cm[i][j], level)
		}
	}

	decodedMask := mpcObj.decodePlainToElemSlots(cryptoParams, rtype, pmMask, maskWidth)

	// v2 sign convention: hub = decode(x + Sum mask) - mask_hub ; others = -mask_i.
	var shareFull mpc_core.RMat
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		shareFull = mpcObj.decodePlainToElemSlots(cryptoParams, rtype, pmSum, maskWidth)
		shareFull.Sub(decodedMask)
	} else {
		shareFull = mpc_core.InitRMat(rtype.Zero(), numCtxRow, maskWidth)
		shareFull.Sub(decodedMask)
	}

	// Output only the nElemRow data slots (element k lives at flat index k).
	for i := range rm {
		for j := 0; j < nElemRow; j++ {
			rm[i][j] = shareFull[i][j]
		}
	}
	return
}

func (mpcObj *MPC) CVecToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, cv crypto.CipherVector, sourcePid, numCtx, nElem int) (rm mpc_core.RVec) {
	return mpcObj.CMatToSS(cryptoParams, rtype, crypto.CipherMatrix{cv}, sourcePid, 1, numCtx, nElem)[0]
}

func (mpcObj *MPC) CiphertextToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, ct *rlwe.Ciphertext, sourcePid, n int) (rv mpc_core.RVec) {
	return mpcObj.CVecToSS(cryptoParams, rtype, crypto.CipherVector{ct}, sourcePid, 1, n)
}
