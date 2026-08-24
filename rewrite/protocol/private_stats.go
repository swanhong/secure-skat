package protocol

import (
	"math"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

type gvLocalGene struct {
	signedWeight []float64
	wGtx         *mat.Dense
	cross        *mat.Dense
	wGvGtG       *mat.Dense
}

type gvGeneValues struct {
	wGtxSquare   *mat.Dense
	privateXtb   []float64
	burdenCross  []float64
	burdenSquare float64
}

type gvPhenoValues struct {
	wGtySquare float64
	wGtxGty    []float64
	wGtySum    float64
}

type gvGeneShares struct {
	wGtxSquare   mpc_core.RMat
	privateXtb   mpc_core.RVec
	burdenCross  mpc_core.RVec
	burdenSquare mpc_core.RElem
}

type gvPhenoShares struct {
	wGtySquare mpc_core.RVec
	wGtxGty    mpc_core.RMat
	wGtySum    mpc_core.RVec
}

func privateSignedWeights(counts []float64, sampleCount int) []float64 {
	/*
		For each B-private variant j with ALT count count[j]:

		base[j]         = max(count[j], 2*N-count[j]) / (2*N)
		weight[j]       = 25 * base[j]^24
		signedWeight[j] = -weight[j] if count[j] > N, otherwise weight[j]
	*/
	totalAlleles := float64(2 * sampleCount)
	signedWeight := make([]float64, len(counts))
	for variant, count := range counts {
		baseCount := totalAlleles - count
		sign := 1.0
		if count > float64(sampleCount) {
			baseCount = count
			sign = -1.0
		}
		signedWeight[variant] = sign * 25 * math.Pow(baseCount/totalAlleles, 24)
	}
	return signedWeight
}

func PrepareGvGeneTerms(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	gv []*mat.Dense,
	gp []*mat.Dense,
	x *mat.Dense,
) (gvLocal []gvLocalGene, gvGene []gvGeneShares) {
	/*
		For each gene g, let D[g] = Diag(signedWeight[g]):

		wGtx[g]              = D[g] * Transpose(Gv[g]) * X[B]
		wGtxSquare[g]        = Transpose(wGtx[g]) * wGtx[g]
		PrivateXtb[g]        = ColSum(wGtx[g])
		cross[g]             = Transpose(Gp[B][g]) * Gv[g] * D[g]
		wGvGtG[g]            = D[g] * Transpose(Gv[g]) * Gv[g] * D[g]
		BurdenCross[g]       = RowSum(cross[g])
		BurdenSquare[g]      = Sum(wGvGtG[g])
	*/
	geneCount := len(batch.GeneIndices)
	gvLocal = make([]gvLocalGene, geneCount)
	localGene := make([]gvGeneValues, geneCount)

	if mpcObj.GetPid() == cohortBPartyID {
		for position, geneIndex := range batch.GeneIndices {
			publicVariantCount := dataParams.Genes[geneIndex].VariantCount
			privateGenotypes := gv[geneIndex]
			if privateGenotypes == nil {
				continue
			}

			_, privateVariantCount := privateGenotypes.Dims()
			counts := make([]float64, privateVariantCount)
			for variant := range counts {
				counts[variant] = mat.Sum(privateGenotypes.ColView(variant))
			}
			signedWeight := privateSignedWeights(counts, dataParams.N)

			gtx := new(mat.Dense)
			gtx.Mul(privateGenotypes.T(), x)
			wGtx := scaleDenseRows(gtx, signedWeight)

			wGtxSquare := new(mat.Dense)
			wGtxSquare.Mul(wGtx.T(), wGtx)
			privateXtb := denseColumnSums(wGtx)

			var cross *mat.Dense
			burdenCross := make([]float64, publicVariantCount)
			if publicVariantCount > 0 {
				cross = new(mat.Dense)
				cross.Mul(gp[geneIndex].T(), privateGenotypes)
				cross = scaleDenseColumns(cross, signedWeight)
				burdenCross = denseRowSums(cross)
			}

			gvGtG := new(mat.Dense)
			gvGtG.Mul(privateGenotypes.T(), privateGenotypes)
			wGvGtG := scaleDenseColumns(gvGtG, signedWeight)
			wGvGtG = scaleDenseRows(wGvGtG, signedWeight)

			gvLocal[position] = gvLocalGene{
				signedWeight: signedWeight,
				wGtx:         wGtx,
				cross:        cross,
				wGvGtG:       wGvGtG,
			}
			localGene[position] = gvGeneValues{
				wGtxSquare:   wGtxSquare,
				privateXtb:   privateXtb,
				burdenCross:  burdenCross,
				burdenSquare: mat.Sum(wGvGtG),
			}
		}
	}

	packedLocalGene := packGvGeneValues(dataParams, batch, localGene)
	_, sharedLength := packedLocalGene.Dims()
	shared := ShareSum(mpcObj, packedLocalGene, 1, sharedLength)[0]
	gvGene = unpackGvGeneShares(mpcObj, dataParams, batch, shared)
	return gvLocal, gvGene
}

func ComputeGvPhenoTerms(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	gv []*mat.Dense,
	y0 mat.Vector,
	gvLocal []gvLocalGene,
) gvPhenoShares {
	/*
		For each gene g and current phenotype t:

		wGty[g,t]       = signedWeight[g] .* (Transpose(Gv[g]) * y0[B,t])
		wGtySquare[g,t] = Dot(wGty[g,t], wGty[g,t])
		wGtxGty[g,t]    = Transpose(wGtx[g]) * wGty[g,t]
		wGtySum[g,t]    = Sum(wGty[g,t])
	*/
	geneCount := len(batch.GeneIndices)
	covariateCount := dataParams.C
	localPheno := make([]gvPhenoValues, geneCount)

	if mpcObj.GetPid() == cohortBPartyID {
		for position, geneIndex := range batch.GeneIndices {
			privateGenotypes := gv[geneIndex]
			if privateGenotypes == nil {
				continue
			}

			_, privateVariantCount := privateGenotypes.Dims()
			wGty := make([]float64, privateVariantCount)
			for variant := range wGty {
				gty := mat.Dot(privateGenotypes.ColView(variant), y0)
				wGty[variant] = gvLocal[position].signedWeight[variant] * gty
			}

			wGtySquare := 0.0
			wGtySum := 0.0
			for _, value := range wGty {
				wGtySquare += value * value
				wGtySum += value
			}

			wGtyVector := mat.NewVecDense(privateVariantCount, wGty)
			wGtxGty := new(mat.VecDense)
			wGtxGty.MulVec(
				gvLocal[position].wGtx.T(), wGtyVector,
			)
			wGtxGtyValues := make([]float64, covariateCount)
			for covariate := range wGtxGtyValues {
				wGtxGtyValues[covariate] = wGtxGty.AtVec(covariate)
			}
			localPheno[position] = gvPhenoValues{
				wGtySquare: wGtySquare,
				wGtxGty:    wGtxGtyValues,
				wGtySum:    wGtySum,
			}
		}
	}

	packedLocalPheno := packGvPhenoValues(localPheno, covariateCount)
	_, sharedLength := packedLocalPheno.Dims()
	shared := ShareSum(mpcObj, packedLocalPheno, 1, sharedLength)[0]
	return unpackGvPhenoShares(mpcObj, geneCount, covariateCount, shared)
}

func AssembleGvQL(
	mpcObj *mpc.MPC,
	gvPheno gvPhenoShares,
	gvGene []gvGeneShares,
	beta mpc_core.RVec,
) (geneBatchQ, geneBatchL mpc_core.RVec) {
	/*
		For each gene g and current phenotype t:

		Q[g,t] = wGtySquare[g,t]
		         - 2 * Dot(wGtxGty[g,t], beta[:,t])
		         + Dot(beta[:,t], wGtxSquare[g] * beta[:,t])
		L[g,t] = wGtySum[g,t] - Dot(PrivateXtb[g], beta[:,t])
	*/
	geneCount := len(gvPheno.wGtySquare)
	covariateCount := len(beta)
	rtype := mpcObj.GetRType()
	dataBits := mpcObj.GetDataBits()
	fracBits := mpcObj.GetFracBits()

	geneBatchQ = mpc_core.InitRVec(rtype.Zero(), geneCount)
	geneBatchL = mpc_core.InitRVec(rtype.Zero(), geneCount)
	betaColumn := mpc_core.InitRMat(rtype.Zero(), covariateCount, 1)
	for covariate := range beta {
		betaColumn[covariate][0] = beta[covariate].Copy()
	}

	for position := 0; position < geneCount; position++ {
		wGtxGtyBetaTerms := mpcObj.SSMultElemVec(
			gvPheno.wGtxGty[position], beta,
		)
		wGtxGtyBetaTerms = mpcObj.TruncVec(
			wGtxGtyBetaTerms, dataBits, fracBits,
		)
		wGtxGtyBeta := sumShares(rtype, wGtxGtyBetaTerms)

		privateXtbBetaTerms := mpcObj.SSMultElemVec(
			gvGene[position].privateXtb, beta,
		)
		privateXtbBetaTerms = mpcObj.TruncVec(
			privateXtbBetaTerms, dataBits, fracBits,
		)
		privateXtbBeta := sumShares(rtype, privateXtbBetaTerms)

		qMatrixBeta := mpcObj.SSMultMat(
			gvGene[position].wGtxSquare, betaColumn,
		)
		qMatrixBeta = mpcObj.TruncMat(qMatrixBeta, dataBits, fracBits)
		qMatrixBetaColumn := make(mpc_core.RVec, covariateCount)
		for covariate := range qMatrixBetaColumn {
			qMatrixBetaColumn[covariate] = qMatrixBeta[covariate][0]
		}

		qCorrectionTerms := mpcObj.SSMultElemVec(beta, qMatrixBetaColumn)
		qCorrectionTerms = mpcObj.TruncVec(qCorrectionTerms, dataBits, fracBits)
		qCorrection := sumShares(rtype, qCorrectionTerms)

		twicewGtxGtyBeta := wGtxGtyBeta.Mul(rtype.FromInt(2))
		geneBatchQ[position] = gvPheno.wGtySquare[position].Sub(
			twicewGtxGtyBeta,
		).Add(qCorrection)
		geneBatchL[position] = gvPheno.wGtySum[position].Sub(
			privateXtbBeta,
		)
	}
	return geneBatchQ, geneBatchL
}

func ComputeGvQL(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	gv []*mat.Dense,
	y0 mat.Vector,
	gvLocal []gvLocalGene,
	gvGene []gvGeneShares,
	beta mpc_core.RVec,
) (geneBatchQ, geneBatchL mpc_core.RVec) {
	/*
		For the current phenotype t:

		gvPheno[t]       = ComputeGvPhenoTerms(Gv[B], y0[B][:,t], gvLocal)
		(gvQ[t], gvL[t]) = AssembleGvQL(gvPheno[t], gvGene, beta[:,t])
	*/
	gvPheno := ComputeGvPhenoTerms(
		mpcObj, dataParams, batch, gv, y0, gvLocal,
	)
	return AssembleGvQL(mpcObj, gvPheno, gvGene, beta)
}

func packGvGeneValues(
	dataParams DataParams,
	batch GeneBatch,
	localGene []gvGeneValues,
) *mat.Dense {
	/*
		For each gene g, serialize the public-shape tuple:

		(wGtxSquare[g], PrivateXtb[g], BurdenCross[g], BurdenSquare[g])
	*/
	covariateCount := dataParams.C
	sharedLength := 0
	for _, geneIndex := range batch.GeneIndices {
		sharedLength += covariateCount*covariateCount + covariateCount
		sharedLength += dataParams.Genes[geneIndex].VariantCount + 1
	}
	packed := make([]float64, sharedLength)
	offset := 0
	for position, geneIndex := range batch.GeneIndices {
		values := localGene[position]
		if values.wGtxSquare != nil {
			for row := 0; row < covariateCount; row++ {
				for column := 0; column < covariateCount; column++ {
					packed[offset+row*covariateCount+column] = values.wGtxSquare.At(row, column)
				}
			}
		}
		offset += covariateCount * covariateCount

		copy(packed[offset:], values.privateXtb)
		offset += covariateCount
		copy(packed[offset:], values.burdenCross)
		offset += dataParams.Genes[geneIndex].VariantCount
		packed[offset] = values.burdenSquare
		offset++
	}
	return mat.NewDense(1, sharedLength, packed)
}

func packGvPhenoValues(
	localPheno []gvPhenoValues,
	covariateCount int,
) *mat.Dense {
	/*
		For each gene g and current phenotype t, serialize:

		(wGtySquare[g,t], wGtxGty[g,t], wGtySum[g,t])
	*/
	fieldsPerGene := covariateCount + 2
	packed := make([]float64, len(localPheno)*fieldsPerGene)
	for position, values := range localPheno {
		offset := position * fieldsPerGene
		packed[offset] = values.wGtySquare
		copy(packed[offset+1:], values.wGtxGty)
		packed[offset+covariateCount+1] = values.wGtySum
	}
	return mat.NewDense(1, len(packed), packed)
}

func unpackGvGeneShares(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	shared mpc_core.RVec,
) []gvGeneShares {
	covariateCount := dataParams.C
	gvGene := make([]gvGeneShares, len(batch.GeneIndices))
	offset := 0
	for position, geneIndex := range batch.GeneIndices {
		publicVariantCount := dataParams.Genes[geneIndex].VariantCount
		wGtxSquare := mpc_core.InitRMat(
			mpcObj.GetRType().Zero(), covariateCount, covariateCount,
		)
		for row := range wGtxSquare {
			for column := range wGtxSquare[row] {
				wGtxSquare[row][column] = shared[offset].Copy()
				offset++
			}
		}
		privateXtb := shared[offset : offset+covariateCount].Copy()
		offset += covariateCount
		burdenCross := shared[offset : offset+publicVariantCount].Copy()
		offset += publicVariantCount
		burdenSquare := shared[offset].Copy()
		offset++

		gvGene[position] = gvGeneShares{
			wGtxSquare:   wGtxSquare,
			privateXtb:   privateXtb,
			burdenCross:  burdenCross,
			burdenSquare: burdenSquare,
		}
	}
	return gvGene
}

func unpackGvPhenoShares(
	mpcObj *mpc.MPC,
	geneCount int,
	covariateCount int,
	shared mpc_core.RVec,
) gvPhenoShares {
	gvPheno := gvPhenoShares{
		wGtySquare: mpc_core.InitRVec(mpcObj.GetRType().Zero(), geneCount),
		wGtxGty: mpc_core.InitRMat(
			mpcObj.GetRType().Zero(), geneCount, covariateCount,
		),
		wGtySum: mpc_core.InitRVec(mpcObj.GetRType().Zero(), geneCount),
	}
	fieldsPerGene := covariateCount + 2
	for position := 0; position < geneCount; position++ {
		offset := position * fieldsPerGene
		gvPheno.wGtySquare[position] = shared[offset].Copy()
		for covariate := 0; covariate < covariateCount; covariate++ {
			gvPheno.wGtxGty[position][covariate] = shared[offset+1+covariate].Copy()
		}
		gvPheno.wGtySum[position] = shared[offset+covariateCount+1].Copy()
	}
	return gvPheno
}

func scaleDenseRows(matrix *mat.Dense, scale []float64) *mat.Dense {
	rows, columns := matrix.Dims()
	result := mat.NewDense(rows, columns, nil)
	for row := 0; row < rows; row++ {
		for column := 0; column < columns; column++ {
			result.Set(row, column, scale[row]*matrix.At(row, column))
		}
	}
	return result
}

func scaleDenseColumns(matrix *mat.Dense, scale []float64) *mat.Dense {
	rows, columns := matrix.Dims()
	result := mat.NewDense(rows, columns, nil)
	for row := 0; row < rows; row++ {
		for column := 0; column < columns; column++ {
			result.Set(row, column, scale[column]*matrix.At(row, column))
		}
	}
	return result
}

func denseColumnSums(matrix *mat.Dense) []float64 {
	_, columns := matrix.Dims()
	sums := make([]float64, columns)
	for column := range sums {
		sums[column] = mat.Sum(matrix.ColView(column))
	}
	return sums
}

func denseRowSums(matrix *mat.Dense) []float64 {
	rows, _ := matrix.Dims()
	sums := make([]float64, rows)
	for row := range sums {
		sums[row] = mat.Sum(matrix.RowView(row))
	}
	return sums
}

func sumShares(rtype mpc_core.RElem, values mpc_core.RVec) mpc_core.RElem {
	sum := rtype.Zero()
	for _, value := range values {
		sum = sum.Add(value)
	}
	return sum
}
