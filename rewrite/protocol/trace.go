package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

type privateTraceValues struct {
	c0GpGpDiag  []float64
	c1GpGpDiag  []float64
	c0GpGpProbe *mat.Dense
	c0GpX       *mat.Dense
	c1GpX       *mat.Dense
	c0XX        *mat.Dense
	c1XX        *mat.Dense
	c2XX        *mat.Dense
	e1          float64
	e2          float64
	e3          float64
}

func prepareLocalPrivateTraceTerms(
	gp *mat.Dense,
	gv *mat.Dense,
	x *mat.Dense,
	privateSignedWeight []float64,
	correctionProbe *mat.Dense,
) privateTraceValues {
	/*
		For one gene, define:

		GpGv = Transpose(Gp) * Gv
		XGv  = Transpose(X) * Gv
		Dv²  = Diag(privateSignedWeight²)
		H    = Dv² * Transpose(Gv) * Gv

		C_k(A,B) = A * H^k * Dv² * Transpose(B)

		The public-private diagonal terms are computed without
		materializing their m-by-m matrices:

		c0WeightedGpGv = GpGv * Dv²
		c1WeightedGpGv = GpGv * H * Dv²

		c0GpGpDiag[j]
		    = Dot(c0WeightedGpGv[j,:], GpGv[j,:])

		c1GpGpDiag[j]
		    = Dot(c1WeightedGpGv[j,:], GpGv[j,:])

		e1 = Trace(H)
		e2 = Trace(H²)
		e3 = Trace(H³)
	*/
	privateWeightSquared := make([]float64, len(privateSignedWeight))
	for variant, weight := range privateSignedWeight {
		privateWeightSquared[variant] = weight * weight
	}

	gpGv := new(mat.Dense)
	gpGv.Mul(gp.T(), gv)

	xGv := new(mat.Dense)
	xGv.Mul(x.T(), gv)

	gvGtG := new(mat.Dense)
	gvGtG.Mul(gv.T(), gv)
	h := scaleDenseRows(gvGtG, privateWeightSquared)

	gpGvTimesH := new(mat.Dense)
	gpGvTimesH.Mul(gpGv, h)

	c0WeightedGpGv := scaleDenseColumns(
		gpGv,
		privateWeightSquared,
	)
	c1WeightedGpGv := scaleDenseColumns(
		gpGvTimesH,
		privateWeightSquared,
	)

	publicVariantCount, _ := gpGv.Dims()
	c0GpGpDiag := make([]float64, publicVariantCount)
	c1GpGpDiag := make([]float64, publicVariantCount)
	for variant := 0; variant < publicVariantCount; variant++ {
		c0GpGpDiag[variant] = mat.Dot(
			c0WeightedGpGv.RowView(variant),
			gpGv.RowView(variant),
		)
		c1GpGpDiag[variant] = mat.Dot(
			c1WeightedGpGv.RowView(variant),
			gpGv.RowView(variant),
		)
	}

	gpGvTransposeProbe := new(mat.Dense)
	gpGvTransposeProbe.Mul(gpGv.T(), correctionProbe)
	weightedProbe := scaleDenseRows(
		gpGvTransposeProbe,
		privateWeightSquared,
	)
	c0GpGpProbe := new(mat.Dense)
	c0GpGpProbe.Mul(gpGv, weightedProbe)

	hSquared := new(mat.Dense)
	hSquared.Mul(h, h)
	hCubed := new(mat.Dense)
	hCubed.Mul(hSquared, h)

	e1 := 0.0
	e2 := 0.0
	e3 := 0.0
	privateVariantCount, _ := h.Dims()
	for variant := 0; variant < privateVariantCount; variant++ {
		e1 += h.At(variant, variant)
		e2 += hSquared.At(variant, variant)
		e3 += hCubed.At(variant, variant)
	}

	return privateTraceValues{
		c0GpGpDiag:  c0GpGpDiag,
		c1GpGpDiag:  c1GpGpDiag,
		c0GpGpProbe: c0GpGpProbe,
		c0GpX: privateTraceContraction(
			0, gpGv, xGv, h, privateWeightSquared,
		),
		c1GpX: privateTraceContraction(
			1, gpGv, xGv, h, privateWeightSquared,
		),
		c0XX: privateTraceContraction(
			0, xGv, xGv, h, privateWeightSquared,
		),
		c1XX: privateTraceContraction(
			1, xGv, xGv, h, privateWeightSquared,
		),
		c2XX: privateTraceContraction(
			2, xGv, xGv, h, privateWeightSquared,
		),
		e1: e1,
		e2: e2,
		e3: e3,
	}
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
		probeAction -> uScale * Dot(U, probeAction)
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
