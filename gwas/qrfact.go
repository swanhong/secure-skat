package gwas

import (
	"math"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"

	"github.com/hhcho/sfgwas/mpc"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"go.dedis.ch/onet/v3/log"
)

func NetDQRenc(cryptoParams *crypto.CryptoParams, mpcObj *mpc.MPC, A crypto.CipherMatrix, nrowsAll []int) crypto.CipherMatrix {
	debug := false
	useBoolean := mpcObj.GetBooleanShareFlag()

	// A is column-encrypted matrix with rows split amongst parties
	slots := cryptoParams.GetSlots()

	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()
	fracBits := mpcObj.GetFracBits()
	dataBits := mpcObj.GetDataBits()
	nparty := len(nrowsAll)

	nrows := 0
	if pid < len(nrowsAll) {
		nrows = nrowsAll[pid]
	}
	ncols := len(A)
	totN := 0
	for i := 1; i < len(nrowsAll); i++ {
		totN += nrowsAll[i]
	}
	sqrtN := math.Sqrt(float64(totN))
	invN := 1.0 / float64(totN)
	invSqrtN := 1.0 / sqrtN

	// Householder vectors
	vList := make(crypto.CipherMatrix, ncols)

	/* Forward pass: Compute a list of Householder vectors */
	for rTest := range A {
		for cTest := range A[rTest] {
			if A[rTest][cTest] == nil {
				log.LLvl1("EXTREME CRITICAL ERROR: NetDQRenc receives A[", rTest, "][", cTest, "] as NULL inside PID", pid)
			}
		}
	}

	for col := 0; col < ncols; col++ {
		log.LLvl1(time.Now().Format(time.RFC3339), "DistriubtedQR, forward, column", col+1, "/", ncols)

		ncolCurr := ncols - col

		// Determine who has the current first row of A and their local indices
		upid, ctid, slotid := crypto.GlobalToPartyIndex(cryptoParams, nrowsAll, col, nparty)

		if true {
			log.LLvl1(time.Now().Format(time.RFC3339), "check location: upid, ctid, slotid", upid, ctid, slotid)
		}

		// Compute Householder vector from the current first column
		var z, zloc *rlwe.Ciphertext
		var uvec crypto.CipherVector

		// Compute zloc and aggregate to get z
		// Consumes one level
		if pid > 0 {

			zloc = crypto.SqSum(cryptoParams, A[0])
			uvec = crypto.CopyEncryptedVector(A[0])

			z = mpcObj.Network.AggregateSharesCT(cryptoParams, zloc)

		} else {
			z = crypto.Zero(cryptoParams)
		}

		// All parties compute sqrt (need party 0 for mpc)
		zSS := mpcObj.CiphertextToSS(cryptoParams, mpcObj.GetRType(), z, -1, 1)

		zSqrtSS, _ := mpcObj.SqrtAndSqrtInverse(zSS, useBoolean)

		var ssIn *rlwe.Ciphertext
		if pid > 0 && upid == pid {
			ssIn = uvec[ctid]
		}

		xSS := mpcObj.CiphertextToSS(cryptoParams, mpcObj.GetRType(), ssIn, upid, slots)[slotid]

		sgnSS := mpcObj.IsPositive(mpc_core.RVec{xSS}, useBoolean)[0]

		// Convert 0, 1 to -1, 1
		sgnSS = sgnSS.Mul(rtype.FromInt(2))
		if pid == mpcObj.GetHubPid() {
			sgnSS = sgnSS.Sub(rtype.FromInt(1))
		}

		alphaSS := mpcObj.SSMultElemVec(zSqrtSS, mpc_core.RVec{sgnSS})

		concatSS := mpc_core.RVec{alphaSS[0].Add(xSS), xSS}
		concatSqSS := mpcObj.SSMultElemVec(concatSS, concatSS)
		zUpdateSS := mpcObj.Trunc(concatSqSS[0].Sub(concatSqSS[1]), dataBits, fracBits)

		// Back to top level
		_, zNewSqrtInvSS := mpcObj.SqrtAndSqrtInverse(mpc_core.RVec{zSS[0].Add(zUpdateSS)}, false)

		scalar := rtype.FromFloat64(sqrtN, mpcObj.GetFracBits())
		zNewSqrtInvSS.MulScalar(scalar)
		zNewSqrtInvSS = mpcObj.TruncVec(zNewSqrtInvSS, dataBits, fracBits)

		alphaScaledSS := mpcObj.SSMultElemVec(alphaSS, zNewSqrtInvSS)
		alphaScaledSS = mpcObj.TruncVec(alphaScaledSS, dataBits, fracBits)

		alphaScaled := mpcObj.SStoCiphertext(cryptoParams, mpc_core.RVec{alphaScaledSS[0]})
		alphaScaled = crypto.Rebalance(cryptoParams, alphaScaled)
		alphaScaled = crypto.Mask(cryptoParams, alphaScaled, slotid, false)

		zNewSqrtInv := mpcObj.SStoCiphertext(cryptoParams, mpc_core.RVec{zNewSqrtInvSS[0]})
		zNewSqrtInv = crypto.Rebalance(cryptoParams, zNewSqrtInv)

		if debug || col <= 2 {
			log.LLvl1(time.Now().Format(time.RFC3339), "col", col,
				"zSS", mpcObj.RevealSym(zSS[0]).Float64(fracBits),
				"zSqrtSS", mpcObj.RevealSym(zSqrtSS[0]).Float64(fracBits),
				"xSS", mpcObj.RevealSym(xSS).Float64(fracBits),
				"sgnSS", mpcObj.RevealSym(sgnSS).Float64(0),
				"alphaSS", mpcObj.RevealSym(alphaSS[0]).Float64(fracBits),
				"zUpdateSS", mpcObj.RevealSym(zUpdateSS).Float64(fracBits),
				"zNewSqrtInvSS", mpcObj.RevealSym(zNewSqrtInvSS[0]).Float64(fracBits),
			)
		}

		if pid > 0 {
			// Consumes one level
			// TODO: postpone scaling until later, after bootstrapping
			// uvec * sqrtN
			uvec = crypto.CMultScalar(cryptoParams, uvec, zNewSqrtInv) //crypto.CMultScalar or uvec[p][c]

			if pid == upid {
				alphaScaled = crypto.Mask(cryptoParams, alphaScaled, slotid, false)
				uvec[ctid] = crypto.Add(cryptoParams, uvec[ctid], alphaScaled)
			}
		}

		// Save to list of all Householder vectors (all parties, including pid=0)
		vList[col] = uvec

		if pid > 0 {
			if debug {
				uvecDec := mpcObj.Network.CollectiveDecrypt(cryptoParams, uvec[0], 1)
				log.LLvl1(time.Now().Format(time.RFC3339), "col", col, "householder vector: ", crypto.DecodeFloatVector(cryptoParams, crypto.PlainVector{uvecDec})[:5])
			}

			// Compute 2 * v * v^T * A
			vMat := crypto.CipherMatrix{uvec}

			fn := func(cp *crypto.CryptoParams, a crypto.CipherVector,
				B crypto.CipherMatrix, j int) crypto.CipherVector {
				//prod := crypto.CMult(cp, a, B[j])
				//return crypto.CMultConst(cp, prod, 2)
				return crypto.CMult(cp, a, B[j])
			}

			// TODO: Bootstrap in the middle after computing v^T * A
			vvTA := DCMatMulAAtB(cryptoParams, mpcObj, vMat, A, nrowsAll, ncolCurr, fn)

			// fmt.Println("Finished matrix multiplication")
			for c := range A {
				scaled := crypto.CMultConstRescale(cryptoParams, vvTA[c], -2*invN, false)
				A[c] = crypto.CAdd(cryptoParams, A[c], scaled)
			}

			A = mpcObj.Network.BootstrapMatAll(cryptoParams, A)

			// Remove first column
			A = A[1:]

			// Remove first row
			if upid == pid {
				for c := range A {
					A[c][ctid] = crypto.Mask(cryptoParams, A[c][ctid], slotid, true)
				}
			}

			A, _ = crypto.FlattenLevels(cryptoParams, A)
		}

	}

	/* Backward pass: Construct Q matrix */

	// Initialize Q to identity matrix
	// Scaled by sqrt(N)
	Q := make(crypto.CipherMatrix, ncols)
	if pid > 0 {
		for c := range Q {
			col := make([]float64, nrows)

			upid, ctid, slotid := crypto.GlobalToPartyIndex(cryptoParams, nrowsAll, c, nparty)

			if upid == pid {
				//col[ctid * slots + slotid] = 1.0
				col[ctid*slots+slotid] = sqrtN
			}
			Q[c], _ = crypto.EncryptFloatVector(cryptoParams, col)
		}
	} else {
		// Party 0 needs valid (empty) CipherVectors so that CAdd in the backward pass doesn't panic
		for c := range Q {
			Q[c] = crypto.CZeros(cryptoParams, 1)
		}
	}

	// SYNCHRONIZE Q across parties so that they represent the same polynomial!
	// Because each party has only their locally valid values in the vector, we can just aggregate them
	// and bootstrap them to synchronize the polynomial.
	// WAIT: If Party 1 and Party 2 both use slots 0-999, aggregating them adds them maliciously!
	// We MUST NOT aggregate them because `vvTQ` aggregation will handle the distributed sum correctly!
	if pid > 0 {
		// Q = mpcObj.Network.AggregateCMat(cryptoParams, Q)
		Q = mpcObj.Network.BootstrapMatAll(cryptoParams, Q)

		// Iterate backwards through the list of Householder vectors to update Q
		for j := ncols - 1; j >= 0; j-- {
			log.LLvl1(time.Now().Format(time.RFC3339), "DistriubtedQR, backward, column", j+1, "/", ncols)

			upid, ctid, slotid := crypto.GlobalToPartyIndex(cryptoParams, nrowsAll, j, nparty)
			ncolCurr := ncols - j

			// Compute 2 * v * v^T * Q
			QSlice := make(crypto.CipherMatrix, 1) //check
			QSlice = Q[j:]

			vMat := make(crypto.CipherMatrix, 1) //check
			vMat = crypto.CipherMatrix{vList[j]}

			fn := func(cp *crypto.CryptoParams, a crypto.CipherVector,
				B crypto.CipherMatrix, j int) crypto.CipherVector {
				var cv crypto.CipherVector
				if j == 0 {
					// For the leftmost column, we are multiplying with a column in
					// the identity matrix
					cv = crypto.CZeros(cp, 1)
					if upid == pid { // only the chosen party computes
						// cv[0] = crypto.MaskWithScaling(cp, a[ctid], slotid, false, sqrtN)
						cv[0] = crypto.Mask(cp, a[ctid], slotid, false) // Accounted for lack of sqrtN scaling below
					}
				} else {
					cv = crypto.CMult(cp, a, B[j])
				}

				return cv
			}

			vvTQ := DCMatMulAAtB(cryptoParams, mpcObj, vMat, QSlice, nrowsAll, ncolCurr, fn)

			// Compute (I - 2 * v * v^T) * Q
			if pid > 0 {
				for c := 0; c < ncolCurr; c++ {
					var scalar float64
					if c == 0 {
						scalar = invSqrtN
					} else {
						scalar = invN
					}
					scaled := crypto.CMultConstRescale(cryptoParams, vvTQ[c], -2*scalar, false)
					Q[j+c] = crypto.CAdd(cryptoParams, Q[j+c], scaled)
				}
			}

			tmp := mpcObj.Network.BootstrapMatAll(cryptoParams, Q[j:j+ncolCurr])
			for c := range tmp {
				Q[j+c] = tmp[c]
			}
		}

		// Mask last ciphertext to zero out flanking elements
		// TODO is this necessary?
		for i := range Q {
			for j := range Q[i] {
				var N int
				if j == len(Q[i])-1 {
					N = ((nrowsAll[pid] - 1) % slots) + 1
				} else {
					N = slots
				}
				Q[i][j] = crypto.MaskTrunc(cryptoParams, Q[i][j], N)
			}
		}

		if debug {
			for i := range Q {
				dec := mpcObj.Network.CollectiveDecrypt(cryptoParams, Q[i][0], 1)
				log.LLvl1(time.Now().Format(time.RFC3339), "col", i, "Q[i]: ", crypto.DecodeFloatVector(cryptoParams, crypto.PlainVector{dec})[:5])
			}
		}

		return Q
	}
	return Q //nil for pid=0
}
