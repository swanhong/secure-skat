package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

// SetupNull fits the pooled null model using moments divided by the total N.
// It returns beta, omega = (X^T X / N)^-1, and rssOverN = RSS / N.
func SetupNull(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	x *mat.Dense,
	y0 *mat.Dense,
	observe func(stage string) func(),
) (
	beta mpc_core.RMat,
	omega mpc_core.RMat,
	rssOverN mpc_core.RVec,
) {
	covariateCount := dataParams.C
	phenotypeCount := dataParams.PhenotypeCount

	var localXtx *mat.Dense
	var localXty *mat.Dense
	var localYty *mat.Dense

	done := observe("null_local_equations")
	if mpcObj.GetPid() != auxiliaryPartyID {
		invN := 1 / float64(dataParams.N)
		localXtx = new(mat.Dense)
		localXtx.Mul(x.T(), x)
		localXtx.Scale(invN, localXtx)

		localXty = new(mat.Dense)
		localXty.Mul(x.T(), y0)
		localXty.Scale(invN, localXty)

		localYty = mat.NewDense(1, phenotypeCount, nil)
		for pheno := 0; pheno < phenotypeCount; pheno++ {
			column := y0.ColView(pheno)
			localYty.Set(0, pheno, invN*mat.Dot(column, column))
		}
	}
	done()

	done = observe("null_aggregate_shares")
	xtx := ShareSum(mpcObj, localXtx, covariateCount, covariateCount)
	xty := ShareSum(mpcObj, localXty, covariateCount, phenotypeCount)
	yty := ShareSum(mpcObj, localYty, 1, phenotypeCount)[0]
	done()

	done = observe("null_factor_solve")
	lower, diagInv := FactorSPD(mpcObj, xtx)
	beta = SolveSPD(mpcObj, lower, diagInv, xty)

	rtype := mpcObj.GetRType()
	identity := mpc_core.InitRMat(rtype.Zero(), covariateCount, covariateCount)
	if mpcObj.GetPid() == mpcObj.GetHubPid() {
		one := rtype.FromFloat64(1.0, mpcObj.GetFracBits())
		for diag := 0; diag < covariateCount; diag++ {
			identity[diag][diag] = one.Copy()
		}
	}
	omega = SolveSPD(mpcObj, lower, diagInv, identity)
	done()

	done = observe("null_rss")
	dataBits := mpcObj.GetDataBits()
	fracBits := mpcObj.GetFracBits()

	xtxBeta := mpcObj.SSMultMat(xtx, beta)
	xtxBeta = mpcObj.TruncMat(xtxBeta, dataBits, fracBits)

	xtyBeta := mpcObj.SSMultElemMat(xty, beta)
	xtyBeta = mpcObj.TruncMat(xtyBeta, dataBits, fracBits)

	betaXtxBeta := mpcObj.SSMultElemMat(beta, xtxBeta)
	betaXtxBeta = mpcObj.TruncMat(betaXtxBeta, dataBits, fracBits)

	rssOverN = mpc_core.InitRVec(rtype.Zero(), phenotypeCount)
	for pheno := 0; pheno < phenotypeCount; pheno++ {
		value := yty[pheno].Copy()
		for cov := 0; cov < covariateCount; cov++ {
			value = value.Sub(xtyBeta[cov][pheno])
			value = value.Sub(xtyBeta[cov][pheno])
			value = value.Add(betaXtxBeta[cov][pheno])
		}
		rssOverN[pheno] = value
	}
	done()
	return beta, omega, rssOverN
}
