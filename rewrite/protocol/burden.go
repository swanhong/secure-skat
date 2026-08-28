package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

func ComputeLocalGtG(
	dataParams DataParams,
	batch GeneBatch,
	gp []*mat.Dense,
) []*mat.Dense {
	/*
		For each gene g:
		localGtG[g] = Transpose(Gp[g]) * Gp[g] / N
	*/
	invN := 1 / float64(dataParams.N)
	localGtG := make([]*mat.Dense, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}
		localGtG[position] = new(mat.Dense)
		localGtG[position].Mul(gp[geneIndex].T(), gp[geneIndex])
		localGtG[position].Scale(invN, localGtG[position])
	}
	return localGtG
}

func SharePooledGtx(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	localGtx []*mat.Dense,
) []mpc_core.RMat {
	/*
		For each gene g:
		pooledGtx[g]
		    = (Transpose(Gp[A,g]) * X[A]
		       + Transpose(Gp[B,g]) * X[B]) / N
		The pooled result remains secret-shared.
	*/
	shapes := make([][2]int, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		shapes[position] = [2]int{
			dataParams.Genes[geneIndex].VariantCount,
			dataParams.C,
		}
	}
	return shareDenseMatrices(
		mpcObj, localGtx, shapes, 1/float64(dataParams.N),
	)
}

func PublicGtGActionExact(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	localGtG []*mat.Dense,
	rightMatrix []mpc_core.RMat,
) []mpc_core.RMat {
	/*
		For each gene g and shared right-hand side R[g]:
		pooledGtG[g]
		    = (Transpose(Gp[A,g]) * Gp[A,g]
		       + Transpose(Gp[B,g]) * Gp[B,g]) / N
		gtgRight[g] = pooledGtG[g] * R[g]
	*/
	shapes := make([][2]int, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		shapes[position] = [2]int{variantCount, variantCount}
	}
	pooledGtG := shareDenseMatrices(mpcObj, localGtG, shapes, 1)

	gtgRightMatrix := make([]mpc_core.RMat, len(batch.GeneIndices))
	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}
		product := mpcObj.SSMultMat(pooledGtG[position], rightMatrix[position])
		gtgRightMatrix[position] = mpcObj.TruncMat(
			product, mpcObj.GetDataBits(), mpcObj.GetFracBits(),
		)
	}
	return gtgRightMatrix
}

func ComputeBurdenQuadratic(
	mpcObj *mpc.MPC,
	signedWeight mpc_core.RVec,
	gtgWeight mpc_core.RVec,
	gvGene gvGeneShares,
) mpc_core.RElem {
	/*
		For one gene g, with w = signedWeight:
		gtgWeight = (Transpose(Gp) * Gp / N) * w
		quadratic
		    = Dot(w, gtgWeight)
		    + 2 * Dot(w, BurdenCross)
		    + BurdenSquare
	*/
	gtgTerm := sharedDot(mpcObj, signedWeight, gtgWeight)
	crossTerm := sharedDot(mpcObj, signedWeight, gvGene.burdenCross)
	twiceCross := crossTerm.Mul(mpcObj.GetRType().FromInt(2))
	return gtgTerm.Add(twiceCross).Add(gvGene.burdenSquare)
}

func ComputeBurdenProjectionTerm(
	mpcObj *mpc.MPC,
	pooledGtx mpc_core.RMat,
	signedWeight mpc_core.RVec,
	privateBurdenXtb mpc_core.RVec,
	omega mpc_core.RMat,
) mpc_core.RElem {
	/*
		For one gene g, with w = signedWeight:

		xtb = Transpose(pooledGtx) * w + privateBurdenXtb
		projection = Dot(xtb, omega * xtb)

		With pooledGtx and privateBurdenXtb at 1/N and
		omega = N * (X^T X)^-1, this equals
		b^T X (X^T X)^-1 X^T b / N.
	*/
	xtb := privateBurdenXtb.Copy()
	if len(signedWeight) > 0 {
		publicXtb := mpcObj.SSMultMat(
			pooledGtx.Transpose(), burdenColumn(mpcObj, signedWeight),
		)
		publicXtb = mpcObj.TruncMat(
			publicXtb, mpcObj.GetDataBits(), mpcObj.GetFracBits(),
		)
		for covariate := range xtb {
			xtb[covariate] = xtb[covariate].Add(publicXtb[covariate][0])
		}
	}

	projected := mpcObj.SSMultMat(omega, burdenColumn(mpcObj, xtb))
	projected = mpcObj.TruncMat(
		projected, mpcObj.GetDataBits(), mpcObj.GetFracBits(),
	)
	projectedVector := mpc_core.InitRVec(mpcObj.GetRType().Zero(), len(xtb))
	for covariate := range projectedVector {
		projectedVector[covariate] = projected[covariate][0]
	}
	return sharedDot(mpcObj, xtb, projectedVector)
}

func ComputeBurdenVariance(
	mpcObj *mpc.MPC,
	batch GeneBatch,
	pooledGtx []mpc_core.RMat,
	signedWeight []mpc_core.RVec,
	gtgWeight []mpc_core.RVec,
	gvGene []gvGeneShares,
	omega mpc_core.RMat,
) mpc_core.RVec {
	/*
		For each gene g:
			V[g] = quadratic[g] - projection[g]
		V[g] remains secret-shared.
	*/
	geneBatchV := mpc_core.InitRVec(
		mpcObj.GetRType().Zero(), len(batch.GeneIndices),
	)
	for position, geneIndex := range batch.GeneIndices {
		quadratic := ComputeBurdenQuadratic(
			mpcObj, signedWeight[geneIndex], gtgWeight[position], gvGene[position],
		)
		projection := ComputeBurdenProjectionTerm(
			mpcObj, pooledGtx[position], signedWeight[geneIndex],
			gvGene[position].burdenXtb, omega,
		)
		geneBatchV[position] = quadratic.Sub(projection)
	}
	return geneBatchV
}

func burdenColumn(mpcObj *mpc.MPC, values mpc_core.RVec) mpc_core.RMat {
	column := mpc_core.InitRMat(mpcObj.GetRType().Zero(), len(values), 1)
	for row := range values {
		column[row][0] = values[row].Copy()
	}
	return column
}

func shareDenseMatrices(
	mpcObj *mpc.MPC,
	localMatrices []*mat.Dense,
	shapes [][2]int,
	scale float64,
) []mpc_core.RMat {
	/*
		For each batch position:
			sharedMatrix = localMatrix[A] + localMatrix[B]
		Matrices are serialized in batch-position, row, column order.
	*/
	total := 0
	for _, shape := range shapes {
		total += shape[0] * shape[1]
	}
	sharedMatrices := make([]mpc_core.RMat, len(shapes))
	if total == 0 {
		return sharedMatrices
	}

	var packed *mat.Dense
	if mpcObj.GetPid() != auxiliaryPartyID {
		values := make([]float64, total)
		offset := 0
		for position, shape := range shapes {
			for row := 0; row < shape[0]; row++ {
				for column := 0; column < shape[1]; column++ {
					values[offset] = scale * localMatrices[position].At(row, column)
					offset++
				}
			}
		}
		packed = mat.NewDense(1, total, values)
	}

	shared := ShareSum(mpcObj, packed, 1, total)[0]
	offset := 0
	for position, shape := range shapes {
		if shape[0] == 0 || shape[1] == 0 {
			continue
		}
		matrix := mpc_core.InitRMat(
			mpcObj.GetRType().Zero(), shape[0], shape[1],
		)
		for row := range matrix {
			for column := range matrix[row] {
				matrix[row][column] = shared[offset].Copy()
				offset++
			}
		}
		sharedMatrices[position] = matrix
	}
	return sharedMatrices
}
