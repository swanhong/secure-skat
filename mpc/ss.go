package mpc

import (
	crand "crypto/rand"
	"math/big"

	mpc_core "github.com/hhcho/mpc-core"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

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
	scaled := new(big.Float).SetPrec(600).Mul(f, bigFloatPow2(fracBits))
	half := new(big.Float).SetPrec(600).SetFloat64(0.5)
	if scaled.Sign() >= 0 {
		scaled.Add(scaled, half)
	} else {
		scaled.Sub(scaled, half)
	}
	r, _ := scaled.Int(nil)
	r.Mod(r, rtype.Modulus())
	return rtype.FromBigInt(r)
}

func bigFloatPow2(n int) *big.Float {
	return new(big.Float).SetPrec(600).SetInt(new(big.Int).Lsh(big.NewInt(1), uint(n)))
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

	out := mpc_core.InitRMat(rtype.Zero(), numCtxRow, nElemRow)
	for i := range pm {
		for j := range pm[i] {
			buf := make([]*big.Float, slots)
			for k := range buf {
				buf[k] = new(big.Float)
			}
			if err := cryptoParams.WithEncoder(func(encoder *ckks.Encoder) error {
				return encoder.Decode(pm[i][j], buf)
			}); err != nil {
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

// CMatToSS converts a CipherMatrix encrypting x into secret-shared values (row-major,
// slot-domain). Faithful to the v2 masked-decryption protocol: parties homomorphically
// add their masks, only Enc(x - Sum(mask_i)) is collectively decrypted (so no party sees
// x), and the additive shares reconstruct x over the field.
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

	// Common level so the homomorphic subtraction and decode are well-defined.
	cm, level := crypto.FlattenLevels(cryptoParams, cm)
	scale := cm[0][0].Scale

	// Bound masks by the active CKKS modulus Q_level (as v2 did), not the field modulus: a mask
	// encodes to a coefficient ~mask_real*scale that must stay below Q_level/2, else Lattigo
	// wraps it mod Q_level and the masks stop cancelling. 4x margin; capped at the field bound.
	P := mpcObj.GetNParty()
	qLevel := cryptoParams.Params.RingQ().ModulusAtLevel[level]
	bf := new(big.Float).SetPrec(600).SetInt(qLevel)
	bf.Mul(bf, bigFloatPow2(mpcObj.GetFracBits()))
	bf.Quo(bf, new(big.Float).SetPrec(600).SetFloat64(scale.Float64()*float64(4*(P-1))))
	bound, _ := bf.Int(nil)
	fieldCap := new(big.Int).Quo(rtype.Modulus(), big.NewInt(int64(4*(P-1))))
	if bound.Cmp(fieldCap) > 0 {
		bound = fieldCap
	}
	if bound.Sign() <= 0 {
		panic("CMatToSS: ciphertext level too low to mask securely")
	}

	// Mask the FULL ciphertext width (as v2 did): CollectiveDecryptMat decrypts every slot, so
	// nonzero non-data slots (e.g. rotate-and-add residues) must be padded too, not just data.
	slots := cryptoParams.GetSlots()
	maskWidth := len(cm[0]) * slots
	mask := sampleCenteredMask(rtype, numCtxRow, maskWidth, bound)

	cmMask := mpcObj.encodeElemSlotsToCM(cryptoParams, mask, scale, level)
	cmMaskAggr := mpcObj.Network.AggregateCMat(cryptoParams, cmMask) // Enc(Sum mask_i) at all parties

	// Enc(x - Sum mask_i)
	masked := make(crypto.CipherMatrix, len(cm))
	for i := range cm {
		masked[i] = crypto.CSub(cryptoParams, cm[i], cmMaskAggr[i])
	}

	pm := mpcObj.Network.CollectiveDecryptMat(cryptoParams, masked, -1) // reveals only x - Sum mask_i (fully masked)
	revealed := mpcObj.decodePlainToElemSlots(cryptoParams, rtype, pm, maskWidth)

	// Full-width shares: hub = (x - Sum mask_i) + mask_hub ; others = mask_i. They sum to x.
	var shareFull mpc_core.RMat
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		revealed.Add(mask)
		shareFull = revealed
	} else {
		shareFull = mask
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
