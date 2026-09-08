package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

const (
	auxiliaryPartyID = iota
	cohortAPartyID
	cohortBPartyID
)

func ShareSum(
	mpcObj *mpc.MPC,
	localValue *mat.Dense,
	rows int,
	columns int,
) mpc_core.RMat {
	rtype := mpcObj.GetRType()
	pid := mpcObj.GetPid()

	if pid == auxiliaryPartyID {
		return mpc_core.InitRMat(rtype.Zero(), rows, columns)
	}

	localRows := make([][]float64, rows)
	for row := range localRows {
		localRows[row] = localValue.RawRowView(row)
	}
	local := mpc_core.FloatToRMat(rtype, localRows, mpcObj.GetFracBits())

	pad := mpcObj.Network.Rand.RandMat(rtype, rows, columns)
	masked := local.Copy()
	masked.Sub(pad)

	otherPartyID := cohortAPartyID
	if pid == cohortAPartyID {
		otherPartyID = cohortBPartyID
	}

	var otherMasked mpc_core.RMat
	if pid < otherPartyID {
		otherMasked = mpcObj.Network.ReceiveRMat(
			rtype,
			rows,
			columns,
			otherPartyID,
		)
		mpcObj.Network.SendRData(masked, otherPartyID)
	} else {
		mpcObj.Network.SendRData(masked, otherPartyID)
		otherMasked = mpcObj.Network.ReceiveRMat(
			rtype,
			rows,
			columns,
			otherPartyID,
		)
	}

	share := pad.Copy()
	share.Add(otherMasked)
	return share
}

func sumShares(rtype mpc_core.RElem, values mpc_core.RVec) mpc_core.RElem {
	sum := rtype.Zero()
	for _, value := range values {
		sum = sum.Add(value)
	}
	return sum
}

func sharedDot(
	mpcObj *mpc.MPC,
	left, right mpc_core.RVec,
) mpc_core.RElem {
	if len(left) == 0 {
		return mpcObj.GetRType().Zero()
	}

	products := mpcObj.SSMultElemVec(left, right)
	products = mpcObj.TruncVec(
		products,
		mpcObj.GetDataBits(),
		mpcObj.GetFracBits(),
	)
	return sumShares(mpcObj.GetRType(), products)
}
