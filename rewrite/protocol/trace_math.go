package protocol

import (
	"math/rand"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

func TraceProbe(
	variantCount int,
	batchWidth int,
	probeCount int,
	seed int64,
) (
	traceProbe *mat.Dense,
	correctionProbe *mat.Dense,
	scale float64,
) {
	/*
		Exact:
		    probe = Identity(variantCount)
		    scale = 1

		Hutchinson:
		    probe = deterministic Rademacher matrix
		    scale = 1 / probeCount

		PublicTrace and PrivateCorrection use separate random streams.
	*/
	if batchWidth <= probeCount {
		if variantCount == 0 {
			return nil, nil, 1
		}

		identity := mat.NewDense(
			variantCount,
			variantCount,
			nil,
		)
		for variant := 0; variant < variantCount; variant++ {
			identity.Set(variant, variant, 1)
		}
		return identity, identity, 1
	}

	scale = 1 / float64(probeCount)
	if variantCount == 0 {
		return nil, nil, scale
	}

	newProbe := func(probSeed int64) *mat.Dense {
		rng := rand.New(rand.NewSource(probSeed))
		values := make([]float64, variantCount*probeCount)
		for index := range values {
			values[index] = -1
			if rng.Intn(2) == 1 {
				values[index] = 1
			}
		}

		return mat.NewDense(variantCount, probeCount, values)
	}

	return newProbe(seed), newProbe(seed + 1), scale

}

func scaleSharedRows(
	mpcObj *mpc.MPC,
	weight mpc_core.RVec,
	matrix mpc_core.RMat,
) mpc_core.RMat {
	/*
		Multiply each matrix row by the corresponding shared weight.

		scaled = Diag(weight) * matrix
	*/
	rows, columns := matrix.Dims()
	rowWeight := mpc_core.InitRMat(
		mpcObj.GetRType().Zero(),
		rows,
		columns,
	)

	for row := 0; row < rows; row++ {
		for column := 0; column < columns; column++ {
			rowWeight[row][column] = weight[row].Copy()
		}
	}

	scaled := mpcObj.SSMultElemMat(rowWeight, matrix)
	return mpcObj.TruncMat(
		scaled,
		mpcObj.GetDataBits(),
		mpcObj.GetFracBits(),
	)
}

func ComputeKppRight(
	mpcObj *mpc.MPC,
	weightedRight mpc_core.RMat,
	gtgRight mpc_core.RMat,
	pooledGtx mpc_core.RMat,
	theta mpc_core.RMat,
	weight mpc_core.RVec,
) mpc_core.RMat {
	/*
		Compute Kpp * right from inputs:
			weightedRight = Diag(weight) * right
		    gtgRight      = pooledGtG * weightedRight
		    pooledGtx     = Transpose(Gp) * X
		    theta         = pooledGtx * xtxInv

		1. Remove the covariate projection:
		    mppRight
		        = gtgRight
		        - theta * Transpose(pooledGtx) * weightedRight

		2. Apply the outer weight:
		    KppRight
		        = 0.5 * Diag(weight) * mppRight

		Output:
		    KppRight = Kpp * right
	*/
	if len(weight) == 0 {
		return nil
	}

	dataBits := mpcObj.GetDataBits()
	fracBits := mpcObj.GetFracBits()

	// 1. Remove the covariate projection.
	//   mppRight = gtgRight - theta * Transpose(pooledGtx) * weightedRight
	covariateRight := mpcObj.SSMultMat(pooledGtx.Transpose(), weightedRight)
	covariateRight = mpcObj.TruncMat(covariateRight, dataBits, fracBits)

	projection := mpcObj.SSMultMat(theta, covariateRight)
	projection = mpcObj.TruncMat(projection, dataBits, fracBits)

	mppRight := gtgRight.Copy()
	mppRight.Sub(projection)

	// 2. Apply the outer weight and the factor 0.5.
	//  kppRight = 0.5 * Diag(weight) * mppRight
	kppRight := scaleSharedRows(mpcObj, weight, mppRight)
	kppRight.MulScalar(mpcObj.GetRType().FromFloat64(0.5, fracBits))

	return mpcObj.TruncMat(kppRight, dataBits, fracBits)
}

func privateTraceContraction(
	power int,
	left *mat.Dense,
	right *mat.Dense,
	h *mat.Dense,
	privateWeightSquared []float64,
) *mat.Dense {
	/*
		C_k(A,B) = A * H^k * Dv^2 * Transpose(B)

		H = Dv^2 * Transpose(Gv) * Gv
	*/
	leftTimesPower := mat.DenseCopyOf(left)
	for exponent := 0; exponent < power; exponent++ {
		product := new(mat.Dense)
		product.Mul(leftTimesPower, h)
		leftTimesPower = product
	}

	weighted := scaleDenseColumns(
		leftTimesPower,
		privateWeightSquared,
	)
	result := new(mat.Dense)
	result.Mul(weighted, right.T())
	return result
}

func ComputeMixedThirdTrace(
	mpcObj *mpc.MPC,
	delta3Action mpc_core.RMat,
	correctionProbe *mat.Dense,
	correctionScale float64,
	theta mpc_core.RMat,
	c0XX mpc_core.RMat,
) mpc_core.RElem {
	/*
		Split delta3Action into:

		probeAction | gpXAction | thetaAction
			(m x r)    (m x c)     (m x c)
		probeAction -> correctionScale * Dot(correctionProbe, probeAction)
		gpXAction   -> 2 * Dot(theta, gpXAction)
		thetaAction -> Trace(Transpose(theta) * thetaAction * c0XX)

		mixedTrace
		    = correctionScale * Dot(correctionProbe, probeAction)
		    - 2 * Dot(theta, gpXAction)
		    + Trace(
		        (Transpose(theta) * thetaAction) * c0XX)
	*/
	rtype := mpcObj.GetRType()
	dataBits := mpcObj.GetDataBits()
	fracBits := mpcObj.GetFracBits()

	publicVariantCount, probeCount := correctionProbe.Dims()
	_, covariateCount := theta.Dims()

	probeAction := mpc_core.InitRMat(
		rtype.Zero(), publicVariantCount, probeCount,
	)
	gpXAction := mpc_core.InitRMat(
		rtype.Zero(), publicVariantCount, covariateCount,
	)
	thetaAction := mpc_core.InitRMat(
		rtype.Zero(), publicVariantCount, covariateCount,
	)

	for variant := 0; variant < publicVariantCount; variant++ {
		for probe := 0; probe < probeCount; probe++ {
			probeAction[variant][probe] =
				delta3Action[variant][probe].Copy()
		}
		for covariate := 0; covariate < covariateCount; covariate++ {
			gpXAction[variant][covariate] =
				delta3Action[variant][probeCount+covariate].Copy()
			thetaAction[variant][covariate] =
				delta3Action[variant][probeCount+covariateCount+covariate].Copy()
		}
	}

	probeProducts := mpc_core.InitRVec(
		rtype.Zero(),
		publicVariantCount*probeCount,
	)
	offset := 0
	for variant := 0; variant < publicVariantCount; variant++ {
		for probe := 0; probe < probeCount; probe++ {
			coefficient := rtype.FromFloat64(
				correctionScale*correctionProbe.At(variant, probe),
				fracBits,
			)
			probeProducts[offset] =
				probeAction[variant][probe].Mul(coefficient)
			offset++
		}
	}
	probeProducts = mpcObj.TruncVec(
		probeProducts,
		dataBits,
		fracBits,
	)
	probeTerm := sumShares(rtype, probeProducts)

	thetaValues := mpc_core.InitRVec(
		rtype.Zero(),
		publicVariantCount*covariateCount,
	)
	gpXValues := mpc_core.InitRVec(
		rtype.Zero(),
		publicVariantCount*covariateCount,
	)
	offset = 0
	for variant := 0; variant < publicVariantCount; variant++ {
		for covariate := 0; covariate < covariateCount; covariate++ {
			thetaValues[offset] = theta[variant][covariate].Copy()
			gpXValues[offset] = gpXAction[variant][covariate].Copy()
			offset++
		}
	}
	thetaGpXProducts := mpcObj.SSMultElemVec(
		thetaValues,
		gpXValues,
	)
	thetaGpXProducts = mpcObj.TruncVec(
		thetaGpXProducts,
		dataBits,
		fracBits,
	)
	thetaGpXTerm := sumShares(rtype, thetaGpXProducts)
	twiceThetaGpX := thetaGpXTerm.Mul(rtype.FromInt(2))

	thetaTransposeThetaAction := mpcObj.SSMultMat(
		theta.Transpose(),
		thetaAction,
	)
	thetaTransposeThetaAction = mpcObj.TruncMat(
		thetaTransposeThetaAction,
		dataBits,
		fracBits,
	)

	traceLeft := mpc_core.InitRVec(
		rtype.Zero(),
		covariateCount*covariateCount,
	)
	traceRight := mpc_core.InitRVec(
		rtype.Zero(),
		covariateCount*covariateCount,
	)
	offset = 0
	for row := 0; row < covariateCount; row++ {
		for column := 0; column < covariateCount; column++ {
			traceLeft[offset] =
				thetaTransposeThetaAction[row][column].Copy()
			traceRight[offset] = c0XX[column][row].Copy()
			offset++
		}
	}
	traceProducts := mpcObj.SSMultElemVec(traceLeft, traceRight)
	traceProducts = mpcObj.TruncVec(
		traceProducts,
		dataBits,
		fracBits,
	)
	traceTerm := sumShares(rtype, traceProducts)

	return probeTerm.Sub(twiceThetaGpX).Add(traceTerm)
}
