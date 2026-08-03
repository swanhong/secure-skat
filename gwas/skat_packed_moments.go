package gwas

import (
	"math"

	mpc_core "github.com/hhcho/mpc-core"
)

type sharedPrivateMomentTables struct {
	psi1Diag, xi1Diag               mpc_core.RVec
	psi2, xi2, pi0, pi1, pi2, psi1Z mpc_core.RMat
	s1, s2, s3                      mpc_core.RElem
}

type packedMomentGene struct {
	weight, dp2            mpc_core.RVec
	u, theta               mpc_core.RMat
	diagGamma              mpc_core.RVec
	tauProbe, psiProbe     [][]float64
	tauScale, psiScale     float64
	private                sharedPrivateMomentTables
	u1, acted              mpc_core.RMat
	tau1, tau2, tau3       mpc_core.RElem
	delta1, delta2, delta3 mpc_core.RElem
}

func floatVectorShares(rtype mpc_core.RElem, fracBits int, values []float64, n int) mpc_core.RVec {
	out := mpc_core.InitRVec(rtype.Zero(), n)
	for i := range values {
		out[i] = rtype.FromFloat64(values[i], fracBits)
	}
	return out
}

func floatMatrixShares(rtype mpc_core.RElem, fracBits int, values [][]float64, rows, columns int) mpc_core.RMat {
	out := mpc_core.InitRMat(rtype.Zero(), rows, columns)
	for i := range values {
		for j := range values[i] {
			out[i][j] = rtype.FromFloat64(values[i][j], fracBits)
		}
	}
	return out
}

func sharePrivateMomentTables(rtype mpc_core.RElem, fracBits, m, c, probes int, values privateMomentTables) sharedPrivateMomentTables {
	return sharedPrivateMomentTables{
		psi1Diag: floatVectorShares(rtype, fracBits, values.psi1Diag, m),
		xi1Diag:  floatVectorShares(rtype, fracBits, values.xi1Diag, m),
		psi2:     floatMatrixShares(rtype, fracBits, values.psi2, m, c),
		xi2:      floatMatrixShares(rtype, fracBits, values.xi2, m, c),
		pi0:      floatMatrixShares(rtype, fracBits, values.pi0, c, c),
		pi1:      floatMatrixShares(rtype, fracBits, values.pi1, c, c),
		pi2:      floatMatrixShares(rtype, fracBits, values.pi2, c, c),
		psi1Z:    floatMatrixShares(rtype, fracBits, values.psi1Z, m, probes),
		s1:       rtype.FromFloat64(values.s1, fracBits),
		s2:       rtype.FromFloat64(values.s2, fracBits),
		s3:       rtype.FromFloat64(values.s3, fracBits),
	}
}

func (ast *AssocTest) localGammaShares(window GeneBatchWindow, local []windowLocalContraction) []mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	out := make([]mpc_core.RMat, len(window.Tiles))
	for tileIndex, tile := range window.Tiles {
		out[tileIndex] = mpc_core.InitRMat(rtype.Zero(), tile.Variants, tile.Variants)
		if local[tileIndex].Gamma == nil {
			continue
		}
		for i := 0; i < tile.Variants; i++ {
			for j := 0; j < tile.Variants; j++ {
				out[tileIndex][i][j] = rtype.FromFloat64(local[tileIndex].Gamma.At(i, j), mpcObj.GetFracBits())
			}
		}
	}
	return out
}

func (ast *AssocTest) rawGammaWeights(window GeneBatchWindow, local []windowLocalContraction, weights []mpc_core.RVec) []mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	out := make([]mpc_core.RVec, len(window.Tiles))
	for tileIndex, tile := range window.Tiles {
		out[tileIndex] = mpc_core.InitRVec(rtype.Zero(), tile.Variants)
		if tile.Variants == 0 {
			continue
		}
		weight := asCol(weights[tileIndex])
		for start := 0; start < tile.Variants; start += gtgChunkRows {
			end := min(start+gtgChunkRows, tile.Variants)
			chunk := mpc_core.InitRMat(rtype.Zero(), end-start, tile.Variants)
			if local[tileIndex].Gamma != nil {
				for i := start; i < end; i++ {
					for j := 0; j < tile.Variants; j++ {
						chunk[i-start][j] = rtype.FromFloat64(local[tileIndex].Gamma.At(i, j), mpcObj.GetFracBits())
					}
				}
			}
			product := mpcObj.SSMultMat(chunk, weight)
			for i := range product {
				out[tileIndex][start+i] = product[i][0]
			}
		}
		out[tileIndex] = mpcObj.TruncVec(out[tileIndex], mpcObj.GetDataBits(), mpcObj.GetFracBits())
	}
	return out
}

func (ast *AssocTest) preparePackedMomentGenes(window GeneBatchWindow, local []windowLocalContraction, weights []mpc_core.RVec, null skatNull, probes int) []packedMomentGene {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	u := ast.packedWindowU(window, local, null.c)
	var theta []mpc_core.RMat
	if probes > 0 {
		theta = ast.packedWindowTheta(u, null)
	}

	states := make([]packedMomentGene, len(window.Tiles))
	var dp2Columns []mpc_core.RMat
	if probes > 0 {
		weightColumns := make([]mpc_core.RMat, len(weights))
		for tile := range weights {
			weightColumns[tile] = asCol(weights[tile])
		}
		dp2Columns = ast.batchElemMul(weightColumns, weightColumns)
	}
	N := float64(ast.skatTotalNumInds())
	for tileIndex, tile := range window.Tiles {
		state := &states[tileIndex]
		state.weight = weights[tileIndex]
		state.u = u[tileIndex]
		if probes == 0 {
			continue
		}
		state.dp2 = col0(dp2Columns[tileIndex])
		state.theta = theta[tileIndex]
		state.diagGamma = mpc_core.InitRVec(rtype.Zero(), tile.Variants)
		if local[tileIndex].Gamma != nil {
			for j := 0; j < tile.Variants; j++ {
				state.diagGamma[j] = rtype.FromFloat64(local[tileIndex].Gamma.At(j, j), mpcObj.GetFracBits())
			}
		}
		state.tauProbe, state.tauScale, _ = skatTraceProbes(tile.Variants, probes, int64(tile.Gene)*1000003+1)
		state.psiProbe, state.psiScale = state.tauProbe, state.tauScale
		if tile.Variants > probes {
			state.psiProbe, state.psiScale, _ = skatTraceProbes(tile.Variants, probes, int64(tile.Gene)*1000003+7)
		}
		tables := makePrivateMomentTables(local[tileIndex].Private, tile.Variants, null.c, N, state.psiProbe)
		probeCount := 0
		if tile.Variants > 0 {
			probeCount = len(state.psiProbe[0])
		}
		state.private = sharePrivateMomentTables(rtype, mpcObj.GetFracBits(), tile.Variants, null.c, probeCount, tables)
	}
	return states
}

func matrixColumns(rtype mpc_core.RElem, matrix mpc_core.RMat, start, end int) mpc_core.RMat {
	out := mpc_core.InitRMat(rtype.Zero(), len(matrix), end-start)
	for i := range matrix {
		copy(out[i], matrix[i][start:end])
	}
	return out
}

func joinMatrixColumns(rtype mpc_core.RElem, matrices ...mpc_core.RMat) mpc_core.RMat {
	rows, columns := 0, 0
	for _, matrix := range matrices {
		if len(matrix) > 0 {
			if rows == 0 {
				rows = len(matrix)
			} else if len(matrix) != rows {
				panic("matrix row count mismatch")
			}
			columns += len(matrix[0])
		}
	}
	out := mpc_core.InitRMat(rtype.Zero(), rows, columns)
	for i := 0; i < rows; i++ {
		offset := 0
		for _, matrix := range matrices {
			if len(matrix) > 0 {
				copy(out[i][offset:], matrix[i])
				offset += len(matrix[i])
			}
		}
	}
	return out
}

func probeWeights(rtype mpc_core.RElem, weight mpc_core.RVec, probes [][]float64) mpc_core.RMat {
	columns := 0
	if len(probes) > 0 {
		columns = len(probes[0])
	}
	out := mpc_core.InitRMat(rtype.Zero(), len(weight), columns)
	for i := range weight {
		for j, sign := range probes[i] {
			out[i][j] = weight[i].Copy()
			if sign < 0 {
				out[i][j] = out[i][j].Neg()
			}
		}
	}
	return out
}

func subMatrices(left, right []mpc_core.RMat) []mpc_core.RMat {
	if len(left) != len(right) {
		panic("matrix count mismatch")
	}
	out := make([]mpc_core.RMat, len(left))
	for gene := range left {
		out[gene] = subMat(left[gene], right[gene])
	}
	return out
}

func (ast *AssocTest) packedMAction(states []packedMomentGene, input, gammaInput []mpc_core.RMat) []mpc_core.RMat {
	rtype := ast.general.mpcObj[0].GetRType()
	uT := make([]mpc_core.RMat, len(states))
	theta := make([]mpc_core.RMat, len(states))
	for gene := range states {
		uT[gene] = transposeShares(rtype, states[gene].u)
		theta[gene] = states[gene].theta
	}
	uTY := ast.batchMatMul(uT, input)
	return subMatrices(gammaInput, ast.batchMatMul(theta, uTY))
}

func (ast *AssocTest) scalePublicVector(values mpc_core.RVec, factors []float64) mpc_core.RVec {
	if len(values) != len(factors) {
		panic("public scale count mismatch")
	}
	mpcObj := ast.general.mpcObj[0]
	out := make(mpc_core.RVec, len(values))
	for i := range values {
		out[i] = values[i].Mul(mpcObj.GetRType().FromFloat64(factors[i], mpcObj.GetFracBits()))
	}
	return mpcObj.TruncVec(out, mpcObj.GetDataBits(), mpcObj.GetFracBits())
}

func (ast *AssocTest) packedHutchinsonWave1(bucket GeneBatchBucket, window GeneBatchWindow, local []windowLocalContraction, states []packedMomentGene, c int) []mpc_core.RVec {
	rtype := ast.general.mpcObj[0].GetRType()
	R := len(states[0].tauProbe[0])
	basis := make([]mpc_core.RMat, len(states))
	dp2 := make([]mpc_core.RVec, len(states))
	for gene := range states {
		basis[gene] = joinMatrixColumns(rtype, states[gene].private.psi1Z, states[gene].private.psi2, states[gene].theta)
		dp2[gene] = states[gene].dp2
	}
	weightedBasis := ast.batchScaleRows(dp2, basis)
	input := make([]mpc_core.RMat, len(states))
	for gene := range states {
		input[gene] = joinMatrixColumns(rtype,
			probeWeights(rtype, states[gene].weight, states[gene].tauProbe),
			weightedBasis[gene], asCol(states[gene].weight))
	}
	gammaInput := ast.packedGramAction(bucket, window, local, input, 2*R+2*c+1)
	mInput := ast.packedMAction(states, input, gammaInput)

	first := make([]mpc_core.RMat, len(states))
	correction := make([]mpc_core.RMat, len(states))
	for gene := range states {
		first[gene] = matrixColumns(rtype, mInput[gene], 0, R)
		correction[gene] = matrixColumns(rtype, mInput[gene], R, 2*R+2*c)
	}
	u1 := ast.batchScaleMatrices(ast.batchScaleRows(stateWeights(states), first), 0.5)
	acted := ast.batchScaleRows(dp2, correction)

	diagonal := make([]mpc_core.RVec, len(states))
	d2Theta := make([]mpc_core.RMat, len(states))
	u := make([]mpc_core.RMat, len(states))
	for gene := range states {
		diagonal[gene] = states[gene].diagGamma
		d2Theta[gene] = matrixColumns(rtype, input[gene], 2*R+c, 2*R+2*c)
		u[gene] = states[gene].u
	}
	tau1Left := ast.batchVectorDots(dp2, diagonal)
	tau1Right := ast.batchMatrixDots(d2Theta, u)
	tau1 := make(mpc_core.RVec, len(states))
	for gene := range states {
		tau1[gene] = tau1Left[gene].Sub(tau1Right[gene])
	}
	tau1 = ast.ssPMul(tau1, 0.5)
	tau2 := ast.batchMatrixDots(u1, u1)
	factors := make([]float64, len(states))
	for gene := range states {
		factors[gene] = states[gene].tauScale
	}
	tau2 = ast.scalePublicVector(tau2, factors)

	gammaW := make([]mpc_core.RVec, len(states))
	for gene := range states {
		states[gene].u1 = u1[gene]
		states[gene].acted = acted[gene]
		states[gene].tau1 = tau1[gene]
		states[gene].tau2 = tau2[gene]
		gammaW[gene] = col0(matrixColumns(rtype, gammaInput[gene], 2*R+2*c, 2*R+2*c+1))
	}
	return gammaW
}

func stateWeights(states []packedMomentGene) []mpc_core.RVec {
	weights := make([]mpc_core.RVec, len(states))
	for gene := range states {
		weights[gene] = states[gene].weight
	}
	return weights
}

func (ast *AssocTest) packedBurdenVariance(states []packedMomentGene, gammaW []mpc_core.RVec, local []windowLocalContraction, null skatNull) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	c := null.c
	N := float64(ast.skatTotalNumInds())

	d := make([]mpc_core.RVec, len(states))
	qvv := mpc_core.InitRVec(rtype.Zero(), len(states))
	xv := make([]mpc_core.RMat, len(states))
	uT := make([]mpc_core.RMat, len(states))
	w := make([]mpc_core.RMat, len(states))
	for gene := range states {
		m := len(states[gene].weight)
		d[gene] = mpc_core.InitRVec(rtype.Zero(), m)
		xv[gene] = mpc_core.InitRMat(rtype.Zero(), c, 1)
		priv := local[gene].Private
		if priv != nil && len(priv.Weight) > 0 {
			if len(priv.BurdenCross) != m || len(priv.BurdenXtz) != c {
				panic("incomplete private burden cache")
			}
			d[gene] = floatVectorShares(rtype, mpcObj.GetFracBits(), priv.BurdenCross, m)
			qvv[gene] = rtype.FromFloat64(priv.BurdenZZ, mpcObj.GetFracBits())
			for j := 0; j < c; j++ {
				xv[gene][j][0] = rtype.FromFloat64(priv.BurdenXtz[j]/math.Sqrt(N), mpcObj.GetFracBits())
			}
		}
		uT[gene] = transposeShares(rtype, states[gene].u)
		w[gene] = asCol(states[gene].weight)
	}

	qpp := ast.batchVectorDots(stateWeights(states), gammaW)
	qpv := ast.batchVectorDots(stateWeights(states), d)
	x := ast.batchMatMul(uT, w)
	for gene := range x {
		if len(x[gene]) == 0 {
			x[gene] = mpc_core.InitRMat(rtype.Zero(), c, 1)
		}
		for j := 0; j < c; j++ {
			x[gene][j][0] = x[gene][j][0].Add(xv[gene][j][0])
		}
	}
	omp := make([]mpc_core.RMat, len(states))
	for gene := range omp {
		omp[gene] = null.omp
	}
	corr := ast.batchMatrixDots(x, ast.batchMatMul(omp, x))
	out := mpc_core.InitRVec(rtype.Zero(), len(states))
	for gene := range out {
		out[gene] = qpp[gene].Add(qpv[gene]).Add(qpv[gene]).Add(qvv[gene]).Sub(corr[gene])
	}
	return out
}

func (ast *AssocTest) packedExactMoments(states []packedMomentGene, gamma []mpc_core.RMat) {
	rtype := ast.general.mpcObj[0].GetRType()
	uT := make([]mpc_core.RMat, len(states))
	theta := make([]mpc_core.RMat, len(states))
	for gene := range states {
		uT[gene] = transposeShares(rtype, states[gene].u)
		theta[gene] = states[gene].theta
	}
	m := subMatrices(gamma, ast.batchMatMul(theta, uT))
	dm := ast.batchScaleRows(stateWeights(states), m)
	dmdT := make([]mpc_core.RMat, len(states))
	for gene := range states {
		dmdT[gene] = transposeShares(rtype, dm[gene])
	}
	dmdT = ast.batchScaleRows(stateWeights(states), dmdT)
	k := make([]mpc_core.RMat, len(states))
	for gene := range states {
		k[gene] = transposeShares(rtype, dmdT[gene])
	}
	k = ast.batchScaleMatrices(k, 0.5)
	k2 := ast.batchMatMul(k, k)
	tau1 := traceShares(rtype, k)
	tau2 := traceShares(rtype, k2)
	tau3 := ast.batchTraceProducts(k2, k)

	basis := make([]mpc_core.RMat, len(states))
	for gene := range states {
		basis[gene] = joinMatrixColumns(rtype, states[gene].private.psi1Z, states[gene].private.psi2, states[gene].theta)
	}
	acted := ast.batchScaleRows(stateWeights(states), basis)
	acted = ast.batchMatMul(k, acted)
	acted = ast.batchScaleRows(stateWeights(states), acted)
	for gene := range acted {
		for i := range acted[gene] {
			acted[gene][i].MulScalar(rtype.FromInt(2))
		}
		states[gene].acted = acted[gene]
		states[gene].tau1 = tau1[gene]
		states[gene].tau2 = tau2[gene]
		states[gene].tau3 = tau3[gene]
	}
}

func (ast *AssocTest) packedHutchinsonWave2(bucket GeneBatchBucket, window GeneBatchWindow, local []windowLocalContraction, states []packedMomentGene, u []mpc_core.RMat) {
	u1 := make([]mpc_core.RMat, len(states))
	for gene := range states {
		u1[gene] = states[gene].u1
	}
	input := ast.batchScaleRows(stateWeights(states), u1)
	rhs := len(input[0][0])
	gammaInput := ast.packedGramAction(bucket, window, local, input, rhs)
	wave2States := make([]packedMomentGene, len(states))
	for gene := range states {
		wave2States[gene].u = u[gene]
		wave2States[gene].theta = states[gene].theta
	}
	mInput := ast.packedMAction(wave2States, input, gammaInput)
	u2 := ast.batchScaleMatrices(ast.batchScaleRows(stateWeights(states), mInput), 0.5)
	factors := make([]float64, len(states))
	for gene := range states {
		factors[gene] = states[gene].tauScale
	}
	tau3 := ast.scalePublicVector(ast.batchMatrixDots(u1, u2), factors)
	for gene := range states {
		states[gene].tau3 = tau3[gene]
		states[gene].weight = nil
		states[gene].u = nil
		states[gene].theta = nil
		states[gene].u1 = nil
		states[gene].tauScale = 0
	}
}

func (ast *AssocTest) packedMomentCorrections(states []packedMomentGene, null skatNull) {
	rtype := ast.general.mpcObj[0].GetRType()
	n := len(states)
	theta := make([]mpc_core.RMat, n)
	dp2 := make([]mpc_core.RVec, n)
	psi1Diag := make([]mpc_core.RVec, n)
	xi1Diag := make([]mpc_core.RVec, n)
	psi2 := make([]mpc_core.RMat, n)
	xi2 := make([]mpc_core.RMat, n)
	pi0 := make([]mpc_core.RMat, n)
	pi1 := make([]mpc_core.RMat, n)
	pi2 := make([]mpc_core.RMat, n)
	omp := make([]mpc_core.RMat, n)
	for gene := range states {
		theta[gene] = states[gene].theta
		dp2[gene] = states[gene].dp2
		psi1Diag[gene] = states[gene].private.psi1Diag
		xi1Diag[gene] = states[gene].private.xi1Diag
		psi2[gene] = states[gene].private.psi2
		xi2[gene] = states[gene].private.xi2
		pi0[gene] = states[gene].private.pi0
		pi1[gene] = states[gene].private.pi1
		pi2[gene] = states[gene].private.pi2
		omp[gene] = null.omp
	}

	omPi0 := ast.batchMatMul(omp, pi0)
	rd1 := ast.batchRowDots(psi2, theta)
	rd2 := ast.batchRowDots(ast.batchMatMul(theta, pi0), theta)
	diagC := make([]mpc_core.RVec, n)
	for gene := range states {
		diagC[gene] = mpc_core.InitRVec(rtype.Zero(), len(dp2[gene]))
		for j := range diagC[gene] {
			diagC[gene][j] = psi1Diag[gene][j].Sub(rd1[gene][j]).Sub(rd1[gene][j]).Add(rd2[gene][j])
		}
	}
	trDp2C := ast.batchVectorDots(dp2, diagC)

	psi2Omp := ast.batchMatMul(psi2, omp)
	termB := subMatrices(xi2, ast.batchMatMul(psi2, omPi0))
	termD := subMatrices(pi1, ast.batchMatMul(pi0, omPi0))
	rdA := ast.batchRowDots(psi2Omp, psi2)
	rdB := ast.batchRowDots(termB, theta)
	rdD := ast.batchRowDots(ast.batchMatMul(theta, termD), theta)
	diagC2 := make([]mpc_core.RVec, n)
	for gene := range states {
		diagC2[gene] = mpc_core.InitRVec(rtype.Zero(), len(dp2[gene]))
		for j := range diagC2[gene] {
			diagC2[gene][j] = xi1Diag[gene][j].Sub(rdA[gene][j]).Sub(rdB[gene][j]).Sub(rdB[gene][j]).Add(rdD[gene][j])
		}
	}
	trDp2C2 := ast.batchVectorDots(dp2, diagC2)

	psiAcc := mpc_core.InitRVec(rtype.Zero(), n)
	dmv2 := make([]mpc_core.RMat, n)
	dmvt := make([]mpc_core.RMat, n)
	psiFactors := make([]float64, n)
	for gene := range states {
		probes := 0
		if len(states[gene].psiProbe) > 0 {
			probes = len(states[gene].psiProbe[0])
		}
		for j := range states[gene].psiProbe {
			for p, sign := range states[gene].psiProbe[j] {
				if sign > 0 {
					psiAcc[gene] = psiAcc[gene].Add(states[gene].acted[j][p])
				} else if sign < 0 {
					psiAcc[gene] = psiAcc[gene].Sub(states[gene].acted[j][p])
				}
			}
		}
		dmv2[gene] = matrixColumns(rtype, states[gene].acted, probes, probes+null.c)
		dmvt[gene] = matrixColumns(rtype, states[gene].acted, probes+null.c, probes+2*null.c)
		psiFactors[gene] = states[gene].psiScale
	}
	trPsi1 := ast.scalePublicVector(psiAcc, psiFactors)
	quadPsi2 := ast.batchMatrixDots(theta, dmv2)
	thetaT := make([]mpc_core.RMat, n)
	for gene := range states {
		thetaT[gene] = transposeShares(rtype, theta[gene])
	}
	quadTheta := ast.batchMatMul(thetaT, dmvt)
	for gene := range quadTheta {
		if len(quadTheta[gene]) == 0 {
			quadTheta[gene] = mpc_core.InitRMat(rtype.Zero(), null.c, null.c)
		}
	}
	quadThetaPi0 := ast.batchTraceProducts(quadTheta, pi0)
	trDp2MppDp2C := mpc_core.InitRVec(rtype.Zero(), n)
	for gene := range states {
		trDp2MppDp2C[gene] = trPsi1[gene].Sub(quadPsi2[gene].Mul(rtype.FromInt(2))).Add(quadThetaPi0[gene])
	}

	omPi1 := ast.batchMatMul(omp, pi1)
	m0sq := ast.batchMatMul(omPi0, omPi0)
	trOmPi0 := traceShares(rtype, omPi0)
	trOmPi1 := traceShares(rtype, omPi1)
	trOmPi2 := ast.batchTraceProducts(omp, pi2)
	trM0sq := traceShares(rtype, m0sq)
	trM0cu := ast.batchTraceProducts(m0sq, omPi0)
	trP1OP0O := ast.batchTraceProducts(omPi1, omPi0)

	d1Base := mpc_core.InitRVec(rtype.Zero(), n)
	d2Private := mpc_core.InitRVec(rtype.Zero(), n)
	d3Private := mpc_core.InitRVec(rtype.Zero(), n)
	for gene := range states {
		d1Base[gene] = states[gene].private.s1.Sub(trOmPi0[gene])
		d2Private[gene] = states[gene].private.s2.Sub(trOmPi1[gene]).Sub(trOmPi1[gene]).Add(trM0sq[gene])
		d3Private[gene] = states[gene].private.s3.
			Sub(trOmPi2[gene].Mul(rtype.FromInt(3))).
			Add(trP1OP0O[gene].Mul(rtype.FromInt(3))).
			Sub(trM0cu[gene])
	}
	d1 := ast.ssPMul(d1Base, 0.5)
	d2a := ast.ssPMul(trDp2C, 0.5)
	d2b := ast.ssPMul(d2Private, 0.25)
	d3a := ast.ssPMul(trDp2MppDp2C, 3.0/8.0)
	d3b := ast.ssPMul(trDp2C2, 3.0/8.0)
	d3c := ast.ssPMul(d3Private, 1.0/8.0)
	for gene := range states {
		states[gene].delta1 = d1[gene]
		states[gene].delta2 = d2a[gene].Add(d2b[gene])
		states[gene].delta3 = d3a[gene].Add(d3b[gene]).Add(d3c[gene])
		states[gene].dp2 = nil
		states[gene].diagGamma = nil
		states[gene].tauProbe = nil
		states[gene].psiProbe = nil
		states[gene].psiScale = 0
		states[gene].private = sharedPrivateMomentTables{}
		states[gene].acted = nil
		states[gene].u = nil
		if len(states[gene].u1) == 0 {
			states[gene].weight = nil
			states[gene].theta = nil
			states[gene].tauScale = 0
		}
	}
}
