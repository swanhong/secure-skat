package mpc

import (
	mpc_core "github.com/hhcho/mpc-core"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
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

// SSToCMat currently uses a simple reveal-then-encrypt fallback.
// This keeps the v6 migration moving while the original masked share conversion
// path is ported separately.
func (mpcObj *MPC) SSToCMat(cryptoParams *crypto.CryptoParams, rm mpc_core.RMat) (cm crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		cm = make(crypto.CipherMatrix, 1)
		cm[0] = make(crypto.CipherVector, 1)
		return
	}

	numCtxRow := len(rm)
	nElemCol := len(rm[0])
	numCtxCol := 1 + ((nElemCol - 1) / cryptoParams.GetSlots())
	hubPid := mpcObj.GetHubPid()

	revealed := mpcObj.RevealSymMat(rm)
	if mpcObj.GetPid() == hubPid {
		revealedFloat := revealed.ToFloat(mpcObj.GetFracBits())
		cm = make(crypto.CipherMatrix, numCtxRow)
		for i := range cm {
			cm[i], _ = crypto.EncryptFloatVector(cryptoParams, revealedFloat[i])
		}
	}

	return mpcObj.Network.BroadcastCMat(cryptoParams, cm, hubPid, numCtxRow, numCtxCol)
}

func (mpcObj *MPC) SSToCVec(cryptoParams *crypto.CryptoParams, rv mpc_core.RVec) (cv crypto.CipherVector) {
	return mpcObj.SSToCMat(cryptoParams, mpc_core.RMat{rv})[0]
}

func (mpcObj *MPC) SStoCiphertext(cryptoParams *crypto.CryptoParams, rv mpc_core.RVec) *rlwe.Ciphertext {
	return mpcObj.SSToCVec(cryptoParams, rv)[0]
}

// CMatToSS currently uses a simple decrypt-then-re-share fallback in which the
// hub holds the full additive share and the other parties hold zero.
func (mpcObj *MPC) CMatToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, cm crypto.CipherMatrix, sourcePid, numCtxRow, numCtxCol, nElemRow int) (rm mpc_core.RMat) {
	rm = mpc_core.InitRMat(rtype.Zero(), numCtxRow, nElemRow)
	if mpcObj.GetPid() == 0 {
		return
	}

	pm := mpcObj.Network.CollectiveDecryptMat(cryptoParams, cm, sourcePid)
	if mpcObj.GetPid() != mpcObj.GetHubPid() {
		return
	}

	decoded := make([][]float64, numCtxRow)
	for i := range decoded {
		row := crypto.DecodeFloatVector(cryptoParams, pm[i])
		if len(row) > nElemRow {
			row = row[:nElemRow]
		}
		decoded[i] = row
	}

	return mpc_core.FloatToRMat(rtype, decoded, mpcObj.GetFracBits())
}

func (mpcObj *MPC) CVecToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, cv crypto.CipherVector, sourcePid, numCtx, nElem int) (rm mpc_core.RVec) {
	return mpcObj.CMatToSS(cryptoParams, rtype, crypto.CipherMatrix{cv}, sourcePid, 1, numCtx, nElem)[0]
}

func (mpcObj *MPC) CiphertextToSS(cryptoParams *crypto.CryptoParams, rtype mpc_core.RElem, ct *rlwe.Ciphertext, sourcePid, n int) (rv mpc_core.RVec) {
	return mpcObj.CVecToSS(cryptoParams, rtype, crypto.CipherVector{ct}, sourcePid, 1, n)
}
