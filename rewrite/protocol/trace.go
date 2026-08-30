package protocol

import (
	"sort"

	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
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

type privateTraceShares struct {
	c0GpGpDiag  mpc_core.RVec
	c1GpGpDiag  mpc_core.RVec
	c0GpGpProbe mpc_core.RMat
	c0GpX       mpc_core.RMat
	c1GpX       mpc_core.RMat
	c0XX        mpc_core.RMat
	c1XX        mpc_core.RMat
	c2XX        mpc_core.RMat
	e1          mpc_core.RElem
	e2          mpc_core.RElem
	e3          mpc_core.RElem
}

func prepareLocalPrivateTraceTerms(
	gp *mat.Dense,
	gv *mat.Dense,
	x *mat.Dense,
	privateSignedWeight []float64,
	correctionProbe *mat.Dense,
	invN float64,
) privateTraceValues {
	/*
		For one gene, define:

			GpGv = Transpose(Gp) * Gv / N
			XGv  = Transpose(X) * Gv / N
			Dv²  = Diag(privateSignedWeight²)
			H    = Dv² * Transpose(Gv) * Gv / N

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
	xGv := new(mat.Dense)
	xGv.Mul(x.T(), gv)
	xGv.Scale(invN, xGv)

	gvGtG := new(mat.Dense)
	gvGtG.Mul(gv.T(), gv)
	gvGtG.Scale(invN, gvGtG)
	h := scaleDenseRows(gvGtG, privateWeightSquared)

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

	privateValues := privateTraceValues{
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

	// If there is no public data, then return private trace terms only
	if gp == nil {
		return privateValues
	}

	gpGv := new(mat.Dense)
	gpGv.Mul(gp.T(), gv)
	gpGv.Scale(invN, gpGv)

	gpGvH := new(mat.Dense)
	gpGvH.Mul(gpGv, h)

	c0WGpGv := scaleDenseColumns(gpGv, privateWeightSquared)
	c1WGpGv := scaleDenseColumns(gpGvH, privateWeightSquared)

	publicVariantCount, _ := gpGv.Dims()
	privateValues.c0GpGpDiag = make([]float64, publicVariantCount)
	privateValues.c1GpGpDiag = make([]float64, publicVariantCount)
	for variant := 0; variant < publicVariantCount; variant++ {
		privateValues.c0GpGpDiag[variant] =
			mat.Dot(c0WGpGv.RowView(variant), gpGv.RowView(variant))
		privateValues.c1GpGpDiag[variant] =
			mat.Dot(c1WGpGv.RowView(variant), gpGv.RowView(variant))
	}

	gpGvTransposeProbe := new(mat.Dense)
	gpGvTransposeProbe.Mul(gpGv.T(), correctionProbe)
	wProbe := scaleDenseRows(gpGvTransposeProbe, privateWeightSquared)
	privateValues.c0GpGpProbe = new(mat.Dense)
	privateValues.c0GpGpProbe.Mul(gpGv, wProbe)
	privateValues.c0GpX = privateTraceContraction(
		0, gpGv, xGv, h, privateWeightSquared,
	)
	privateValues.c1GpX = privateTraceContraction(
		1, gpGv, xGv, h, privateWeightSquared,
	)

	return privateValues
}

func PreparePrivateTraceTerms(
	mpcObj *mpc.MPC,
	dataParams DataParams,
	batch GeneBatch,
	gp []*mat.Dense,
	gv []*mat.Dense,
	x *mat.Dense,
	gvLocal []gvLocalGene,
	correctionProbe []*mat.Dense,
) []privateTraceShares {
	geneCount := len(batch.GeneIndices)
	covariateCount := dataParams.C
	probeCounts := make([]int, geneCount)
	localValues := make([]privateTraceValues, geneCount)
	packedLength := 0

	// Compute the fixed public packed length from public dimensions
	for position, geneIndex := range batch.GeneIndices {
		publicVariantCount :=
			dataParams.Genes[geneIndex].VariantCount

		if publicVariantCount > 0 {
			_, probeCounts[position] =
				correctionProbe[position].Dims()
		}

		packedLength +=
			2*publicVariantCount +
				publicVariantCount*probeCounts[position]

		packedLength +=
			2*publicVariantCount*covariateCount +
				3*covariateCount*covariateCount +
				3
	}

	// B: Compute the local private trace terms.
	if mpcObj.GetPid() == cohortBPartyID {
		for position, geneIndex := range batch.GeneIndices {
			if gv[geneIndex] == nil {
				continue
			}

			localValues[position] =
				prepareLocalPrivateTraceTerms(
					gp[geneIndex],
					gv[geneIndex],
					x,
					gvLocal[position].signedWeight,
					correctionProbe[position],
					1/float64(dataParams.N),
				)
		}
	}

	// A and B: Allocate the same public-length row; A leaves it zero.
	var packedLocal *mat.Dense
	if mpcObj.GetPid() != auxiliaryPartyID {
		packedLocal = mat.NewDense(
			1,
			packedLength,
			nil,
		)
	}

	// B: Pack the local private trace terms into its row.
	if mpcObj.GetPid() == cohortBPartyID {
		packed := packedLocal.RawRowView(0)
		offset := 0

		packMatrix := func(
			matrix *mat.Dense,
			rows int,
			columns int,
		) {
			if matrix != nil {
				for row := 0; row < rows; row++ {
					for column := 0; column < columns; column++ {
						packed[offset+row*columns+column] =
							matrix.At(row, column)
					}
				}
			}
			offset += rows * columns
		}

		for position, geneIndex := range batch.GeneIndices {
			publicVariantCount :=
				dataParams.Genes[geneIndex].VariantCount
			probeCount := probeCounts[position]
			values := localValues[position]

			copy(
				packed[offset:],
				values.c0GpGpDiag,
			)
			offset += publicVariantCount

			copy(
				packed[offset:],
				values.c1GpGpDiag,
			)
			offset += publicVariantCount

			packMatrix(
				values.c0GpGpProbe,
				publicVariantCount,
				probeCount,
			)
			packMatrix(
				values.c0GpX,
				publicVariantCount,
				covariateCount,
			)
			packMatrix(
				values.c1GpX,
				publicVariantCount,
				covariateCount,
			)
			packMatrix(
				values.c0XX,
				covariateCount,
				covariateCount,
			)
			packMatrix(
				values.c1XX,
				covariateCount,
				covariateCount,
			)
			packMatrix(
				values.c2XX,
				covariateCount,
				covariateCount,
			)

			packed[offset] = values.e1
			packed[offset+1] = values.e2
			packed[offset+2] = values.e3
			offset += 3
		}
	}

	// Share the packed B-local record without revealing it.
	shared := ShareSum(
		mpcObj,
		packedLocal,
		1,
		packedLength,
	)[0]

	// Unpack each party's share into per-gene trace terms.
	terms := make([]privateTraceShares, geneCount)
	offset := 0

	unpackVector := func(length int) mpc_core.RVec {
		vector := shared[offset : offset+length].Copy()
		offset += length
		return vector
	}

	unpackMatrix := func(
		rows int,
		columns int,
	) mpc_core.RMat {
		if rows == 0 {
			return nil
		}

		matrix := mpc_core.InitRMat(
			mpcObj.GetRType().Zero(),
			rows,
			columns,
		)
		for row := 0; row < rows; row++ {
			for column := 0; column < columns; column++ {
				matrix[row][column] =
					shared[offset].Copy()
				offset++
			}
		}
		return matrix
	}

	for position, geneIndex := range batch.GeneIndices {
		publicVariantCount :=
			dataParams.Genes[geneIndex].VariantCount
		probeCount := probeCounts[position]

		terms[position] = privateTraceShares{
			c0GpGpDiag: unpackVector(
				publicVariantCount,
			),
			c1GpGpDiag: unpackVector(
				publicVariantCount,
			),
			c0GpGpProbe: unpackMatrix(
				publicVariantCount,
				probeCount,
			),
			c0GpX: unpackMatrix(
				publicVariantCount,
				covariateCount,
			),
			c1GpX: unpackMatrix(
				publicVariantCount,
				covariateCount,
			),
			c0XX: unpackMatrix(
				covariateCount,
				covariateCount,
			),
			c1XX: unpackMatrix(
				covariateCount,
				covariateCount,
			),
			c2XX: unpackMatrix(
				covariateCount,
				covariateCount,
			),
			e1: shared[offset].Copy(),
			e2: shared[offset+1].Copy(),
			e3: shared[offset+2].Copy(),
		}
		offset += 3
	}

	return terms
}

func ComputeDelta(
	mpcObj *mpc.MPC,
	privateTerms privateTraceShares,
	weightSquared mpc_core.RVec,
	theta mpc_core.RMat,
	omega mpc_core.RMat,
	delta3Action mpc_core.RMat,
	correctionProbe *mat.Dense,
	correctionScale float64,
) (
	delta1 mpc_core.RElem,
	delta2 mpc_core.RElem,
	delta3 mpc_core.RElem,
) {
	/*
		Compute delta1, delta2, and delta3 from the shared terms
	*/
	rtype := mpcObj.GetRType()
	dataBits := mpcObj.GetDataBits()
	fracBits := mpcObj.GetFracBits()
	covariateCount, _ := omega.Dims()

	// *************************
	// 1. Compute privateTrace1, privateTrace2, and privateTrace3.
	// *************************
	// JQk = omega * c_kXX
	jQ0 := mpcObj.SSMultMat(omega, privateTerms.c0XX)
	jQ0 = mpcObj.TruncMat(jQ0, dataBits, fracBits)

	jQ1 := mpcObj.SSMultMat(omega, privateTerms.c1XX)
	jQ1 = mpcObj.TruncMat(jQ1, dataBits, fracBits)

	jQ2 := mpcObj.SSMultMat(omega, privateTerms.c2XX)
	jQ2 = mpcObj.TruncMat(jQ2, dataBits, fracBits)

	// Compute the traces of JQ0, JQ1, JQ2, JQ0^2, JQ1JQ0, and JQ0^3.
	jQ0Squared := mpcObj.SSMultMat(jQ0, jQ0)
	jQ0Squared = mpcObj.TruncMat(jQ0Squared, dataBits, fracBits)

	jQ1JQ0 := mpcObj.SSMultMat(jQ1, jQ0)
	jQ1JQ0 = mpcObj.TruncMat(jQ1JQ0, dataBits, fracBits)

	jQ0Cubed := mpcObj.SSMultMat(jQ0Squared, jQ0)
	jQ0Cubed = mpcObj.TruncMat(jQ0Cubed, dataBits, fracBits)

	// compute traces of the matrices
	traceJQ0 := rtype.Zero()
	traceJQ1 := rtype.Zero()
	traceJQ2 := rtype.Zero()
	traceJQ0Squared := rtype.Zero()
	traceJQ1JQ0 := rtype.Zero()
	traceJQ0Cubed := rtype.Zero()
	for diagonal := 0; diagonal < covariateCount; diagonal++ {
		traceJQ0 = traceJQ0.Add(jQ0[diagonal][diagonal])
		traceJQ1 = traceJQ1.Add(jQ1[diagonal][diagonal])
		traceJQ2 = traceJQ2.Add(jQ2[diagonal][diagonal])
		traceJQ0Squared = traceJQ0Squared.Add(
			jQ0Squared[diagonal][diagonal],
		)
		traceJQ1JQ0 = traceJQ1JQ0.Add(
			jQ1JQ0[diagonal][diagonal],
		)
		traceJQ0Cubed = traceJQ0Cubed.Add(
			jQ0Cubed[diagonal][diagonal],
		)
	}

	// privateTrace1 = e1 - Trace(JQ0)
	// privateTrace2 = e2 - 2*Trace(JQ1) + Trace(JQ0^2)
	// privateTrace3 = e3 - 3*Trace(JQ2) + 3*Trace(JQ1JQ0) - Trace(JQ0^3)
	privateTrace1 := privateTerms.e1.Sub(traceJQ0)
	privateTrace2 := privateTerms.e2.
		Sub(traceJQ1.Mul(rtype.FromInt(2))).
		Add(traceJQ0Squared)
	privateTrace3 := privateTerms.e3.
		Sub(traceJQ2.Mul(rtype.FromInt(3))).
		Add(traceJQ1JQ0.Mul(rtype.FromInt(3))).
		Sub(traceJQ0Cubed)

	// *************************
	// 2. Handle m=0 using private-only terms.
	// *************************
	// if m=0, then just return 0.5*privateTrace1, 0.25*privateTrace2, 0.125*privateTrace3
	if len(weightSquared) == 0 {
		scaled := mpc_core.RVec{
			privateTrace1.Mul(
				rtype.FromFloat64(0.5, fracBits),
			),
			privateTrace2.Mul(
				rtype.FromFloat64(0.25, fracBits),
			),
			privateTrace3.Mul(
				rtype.FromFloat64(0.125, fracBits),
			),
		}
		scaled = mpcObj.TruncVec(scaled, dataBits, fracBits)
		return scaled[0], scaled[1], scaled[2]
	}

	// *************************
	// 3. Compute diag0.
	// *************************
	// diag0 = c0GpGpDiag - 2*Diag(c0GpX * Transpose(theta)) + Diag(theta * c0XX * Transpose(theta))
	p0ThetaTranspose := mpcObj.SSMultMat(
		privateTerms.c0GpX,
		theta.Transpose(),
	)
	p0ThetaTranspose = mpcObj.TruncMat(
		p0ThetaTranspose,
		dataBits,
		fracBits,
	)

	thetaQ0 := mpcObj.SSMultMat(
		theta,
		privateTerms.c0XX,
	)
	thetaQ0 = mpcObj.TruncMat(
		thetaQ0,
		dataBits,
		fracBits,
	)

	thetaQ0ThetaTranspose := mpcObj.SSMultMat(
		thetaQ0,
		theta.Transpose(),
	)
	thetaQ0ThetaTranspose = mpcObj.TruncMat(
		thetaQ0ThetaTranspose,
		dataBits,
		fracBits,
	)

	publicVariantCount := len(weightSquared)
	diag0 := mpc_core.InitRVec(
		rtype.Zero(),
		publicVariantCount,
	)
	for variant := 0; variant < publicVariantCount; variant++ {
		diag0[variant] = privateTerms.c0GpGpDiag[variant].
			Sub(
				p0ThetaTranspose[variant][variant].
					Mul(rtype.FromInt(2)),
			).
			Add(thetaQ0ThetaTranspose[variant][variant])
	}

	// *************************
	// 4. Compute diag1.
	// *************************
	// diag1 = c1GpGpDiag
	//       - Diag(c0GpX * omega * Transpose(c0GpX))
	//       - 2*Diag((c1GpX - c0GpX * omega * c0XX) * Transpose(theta))
	//       + Diag(theta * (c1XX - c0XX * omega * c0XX) * Transpose(theta))
	p0J := mpcObj.SSMultMat(
		privateTerms.c0GpX,
		omega,
	)
	p0J = mpcObj.TruncMat(
		p0J,
		dataBits,
		fracBits,
	)

	p0JP0Transpose := mpcObj.SSMultMat(
		p0J,
		privateTerms.c0GpX.Transpose(),
	)
	p0JP0Transpose = mpcObj.TruncMat(
		p0JP0Transpose,
		dataBits,
		fracBits,
	)

	p0JQ0 := mpcObj.SSMultMat(
		p0J,
		privateTerms.c0XX,
	)
	p0JQ0 = mpcObj.TruncMat(
		p0JQ0,
		dataBits,
		fracBits,
	)

	p1Correction := privateTerms.c1GpX.Copy()
	p1Correction.Sub(p0JQ0)

	p1CorrectionThetaTranspose := mpcObj.SSMultMat(
		p1Correction,
		theta.Transpose(),
	)
	p1CorrectionThetaTranspose = mpcObj.TruncMat(
		p1CorrectionThetaTranspose,
		dataBits,
		fracBits,
	)

	q0JQ0 := mpcObj.SSMultMat(
		privateTerms.c0XX,
		jQ0,
	)
	q0JQ0 = mpcObj.TruncMat(
		q0JQ0,
		dataBits,
		fracBits,
	)

	q1Correction := privateTerms.c1XX.Copy()
	q1Correction.Sub(q0JQ0)

	thetaQ1Correction := mpcObj.SSMultMat(
		theta,
		q1Correction,
	)
	thetaQ1Correction = mpcObj.TruncMat(
		thetaQ1Correction,
		dataBits,
		fracBits,
	)

	thetaQ1CorrectionThetaTranspose := mpcObj.SSMultMat(
		thetaQ1Correction,
		theta.Transpose(),
	)
	thetaQ1CorrectionThetaTranspose = mpcObj.TruncMat(
		thetaQ1CorrectionThetaTranspose,
		dataBits,
		fracBits,
	)

	diag1 := mpc_core.InitRVec(
		rtype.Zero(),
		publicVariantCount,
	)
	for variant := 0; variant < publicVariantCount; variant++ {
		diag1[variant] = privateTerms.c1GpGpDiag[variant].
			Sub(p0JP0Transpose[variant][variant]).
			Sub(
				p1CorrectionThetaTranspose[variant][variant].
					Mul(rtype.FromInt(2)),
			).
			Add(
				thetaQ1CorrectionThetaTranspose[variant][variant],
			)
	}
	// *************************
	// 5. Compute mixedTrace.
	// *************************
	// Compute the mixed public-private third-order trace correction.
	mixedTrace := ComputeMixedThirdTrace(
		mpcObj,
		delta3Action,
		correctionProbe,
		correctionScale,
		theta,
		privateTerms.c0XX,
	)

	// *************************
	// 6. Assemble delta1, delta2, and delta3.
	// *************************
	// delta1 = 0.5*privateTrace1
	// delta2 = 0.5*Dot(weightSquared, diag0) + 0.25*privateTrace2
	// delta3 = 0.375*mixedTrace + 0.375*Dot(weightSquared, diag1) + 0.125*privateTrace3
	diag0Term := sharedDot(
		mpcObj,
		weightSquared,
		diag0,
	)
	diag1Term := sharedDot(
		mpcObj,
		weightSquared,
		diag1,
	)

	scaled := mpc_core.RVec{
		privateTrace1.Mul(
			rtype.FromFloat64(0.5, fracBits),
		),
		diag0Term.Mul(
			rtype.FromFloat64(0.5, fracBits),
		),
		privateTrace2.Mul(
			rtype.FromFloat64(0.25, fracBits),
		),
		mixedTrace.Mul(
			rtype.FromFloat64(0.375, fracBits),
		),
		diag1Term.Mul(
			rtype.FromFloat64(0.375, fracBits),
		),
		privateTrace3.Mul(
			rtype.FromFloat64(0.125, fracBits),
		),
	}
	scaled = mpcObj.TruncVec(
		scaled,
		dataBits,
		fracBits,
	)

	delta1 = scaled[0]
	delta2 = scaled[1].Add(scaled[2])
	delta3 = scaled[3].Add(scaled[4]).Add(scaled[5])
	return delta1, delta2, delta3
}

func GtGTransformGaloisElements(
	heParams ckks.Parameters,
	cryptoParams CryptoParams,
) []uint64 {
	/*
		Compute the rotation keys required by Hutchinson GtG transforms.

		Input:
		    cryptoParams.Batches = packed gene batches

		1. Select Hutchinson batches:
		    batch.W > cryptoParams.R

		2. Add the transform diagonals:
		    diagonal[d] = d * nu
		    nu          = Slots / batch.W

		Output:
		    Galois elements required by the GtG transforms

		The caller must combine these elements with rotation keys
		required by other protocol operations.
	*/
	uniqueElements := make(map[uint64]struct{})

	for _, batch := range cryptoParams.Batches {
		if batch.W <= cryptoParams.R {
			continue
		}

		nu := cryptoParams.Slots / batch.W
		diagonalIndices := make([]int, batch.W)
		for diagonal := range diagonalIndices {
			diagonalIndices[diagonal] = diagonal * nu
		}

		transformParameters := lintrans.Parameters{
			DiagonalsIndexList:        diagonalIndices,
			LevelQ:                    heParams.MaxLevel(),
			LevelP:                    heParams.MaxLevelP(),
			Scale:                     heParams.DefaultScale(),
			LogDimensions:             heParams.LogMaxDimensions(),
			LogBabyStepGiantStepRatio: 0,
		}

		for _, element := range lintrans.GaloisElements(
			heParams,
			transformParameters,
		) {
			uniqueElements[element] = struct{}{}
		}
	}

	elements := make([]uint64, 0, len(uniqueElements))
	for element := range uniqueElements {
		elements = append(elements, element)
	}
	sort.Slice(elements, func(left, right int) bool {
		return elements[left] < elements[right]
	})

	return elements
}

func ApplyGtG(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	localGtG []*mat.Dense,
	rightMatrix []mpc_core.RMat,
	rhsCount int,
) []mpc_core.RMat {
	/*
		Compute normalized pooledGtG * rightMatrix using packed HE.

		Inputs:
		    localGtG[position]   = local Transpose(Gp) * Gp / N
		    rightMatrix[position] = shared right-hand-side matrix
		    rhsCount             = number of right-hand-side columns

		1. Pack the shared right-hand-side matrices:
		    group = column / H
		    slot  = group*Slots
		          + variant*nu
		          + LaneBase
		          + column%H

		2. Convert the packed shares to ciphertexts:
		    ShareToCipher(packedRight)

		3. Apply each cohort's local GtG transform:
		    localResult[i] = localGtG[i] * rightMatrix

		4. Rerandomize each local result:
		    localResult[i] += Enc(0)

		5. Aggregate and convert back to shares:
		    pooledResult = localResult[A] + localResult[B]
		    packedResult = CipherToShare(pooledResult)

		6. Unpack the shared result matrices.

		Output:
		    gtgRightMatrix[position]
		        = (pooledGtG[position] / N) * rightMatrix[position]
	*/
	slots := cryptoParams.Slots
	nu := slots / batch.W
	h := nu / len(batch.GeneIndices)
	groupCount := (rhsCount + h - 1) / h
	rtype := mpcObj.GetRType()

	// 1. Pack the shared right-hand-side matrices.
	packedRight := mpc_core.InitRVec(
		rtype.Zero(),
		groupCount*slots,
	)

	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		laneBase := position * h

		for variant := 0; variant < variantCount; variant++ {
			for column := 0; column < rhsCount; column++ {
				group := column / h
				lane := laneBase + column%h
				slot := group*slots + variant*nu + lane

				packedRight[slot] =
					rightMatrix[position][variant][column].Copy()
			}
		}
	}

	// 2. Convert the packed shares to ciphertexts.
	packedRightCipher := mpcObj.SSToCVec(
		heParams,
		packedRight,
	)
	if mpcObj.GetPid() == auxiliaryPartyID {
		packedRightCipher = make(
			securecrypto.CipherVector,
			groupCount,
		)
	}

	pooledResult := make(
		securecrypto.CipherVector,
		groupCount,
	)

	if mpcObj.GetPid() != auxiliaryPartyID {
		// 3. Encode and apply the cohort-local GtG transform.
		inputLevel := packedRightCipher[0].Level()

		diagonalIndices := make([]int, batch.W)
		diagonals := make(
			lintrans.Diagonals[float64],
			batch.W,
		)

		for diagonal := 0; diagonal < batch.W; diagonal++ {
			diagonalIndex := diagonal * nu
			diagonalIndices[diagonal] = diagonalIndex

			values := make([]float64, slots)

			for position, geneIndex := range batch.GeneIndices {
				variantCount :=
					dataParams.Genes[geneIndex].VariantCount
				laneBase := position * h
				gamma := localGtG[position]

				if gamma == nil {
					continue
				}

				for row := 0; row < variantCount; row++ {
					column := (row + diagonal) % batch.W
					if column >= variantCount {
						continue
					}

					for lane := 0; lane < h; lane++ {
						slot := row*nu + laneBase + lane
						values[slot] = gamma.At(row, column)
					}
				}
			}

			diagonals[diagonalIndex] = values
		}

		transformParameters := lintrans.Parameters{
			DiagonalsIndexList: diagonalIndices,
			LevelQ:             inputLevel,
			LevelP:             heParams.Params.MaxLevelP(),
			Scale: heParams.Params.GetOptimalScalingFactor(
				heParams.Params.DefaultScale(),
				heParams.Params.DefaultScale(),
				inputLevel,
			),
			LogDimensions:             heParams.Params.LogMaxDimensions(),
			LogBabyStepGiantStepRatio: 0,
		}

		transform := lintrans.NewTransformation(
			heParams.Params,
			transformParameters,
		)
		if err := heParams.WithEncoder(
			func(encoder *ckks.Encoder) error {
				return lintrans.Encode(
					encoder,
					diagonals,
					transform,
				)
			},
		); err != nil {
			panic(err)
		}

		localResult := make(
			securecrypto.CipherVector,
			groupCount,
		)
		if err := heParams.WithEvaluator(
			func(evaluator *ckks.Evaluator) error {
				linearEvaluator :=
					lintrans.NewEvaluator(evaluator)

				for group, ciphertext := range packedRightCipher {
					result, err :=
						linearEvaluator.EvaluateNew(
							ciphertext,
							transform,
						)
					if err != nil {
						return err
					}
					if err := evaluator.Rescale(
						result,
						result,
					); err != nil {
						return err
					}

					localResult[group] = result
				}
				return nil
			},
		); err != nil {
			panic(err)
		}

		// 4. Rerandomize the cohort-local result.
		freshZero := securecrypto.CZeros(
			heParams,
			groupCount,
		)
		freshZero = securecrypto.DropLevel(
			heParams,
			securecrypto.CipherMatrix{freshZero},
			localResult[0].Level(),
		)[0]
		localResult = securecrypto.CAdd(
			heParams,
			localResult,
			freshZero,
		)

		// 5. Aggregate the cohort-local results.
		pooledResult = mpcObj.Network.AggregateCVec(
			heParams,
			localResult,
		)
	}

	packedResult := mpcObj.CVecToSS(
		heParams,
		rtype,
		pooledResult,
		-1,
		groupCount,
		groupCount*slots,
	)

	// 6. Unpack the shared result matrices.
	gtgRightMatrix := make(
		[]mpc_core.RMat,
		len(batch.GeneIndices),
	)

	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		if variantCount == 0 {
			continue
		}

		laneBase := position * h
		result := mpc_core.InitRMat(
			rtype.Zero(),
			variantCount,
			rhsCount,
		)

		for variant := 0; variant < variantCount; variant++ {
			for column := 0; column < rhsCount; column++ {
				group := column / h
				lane := laneBase + column%h
				slot := group*slots + variant*nu + lane

				result[variant][column] =
					packedResult[slot].Copy()
			}
		}

		gtgRightMatrix[position] = result
	}

	return gtgRightMatrix
}

func PublicGtGAction(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	localGtG []*mat.Dense,
	rightMatrix []mpc_core.RMat,
	rhsCount int,
) []mpc_core.RMat {
	/*
		Compute normalized pooledGtG * rightMatrix.

		Inputs:
		    localGtG   = cohort-local Transpose(Gp) * Gp / N
		    rightMatrix = shared right-hand-side matrices
		    rhsCount    = number of right-hand-side columns

		1. Select the public execution mode:
		    Exact       if batch.W <= cryptoParams.R
		    Hutchinson  if batch.W > cryptoParams.R

		2. Compute the same semantic action:
		    gtgRightMatrix
		        = (pooledGtG / N) * rightMatrix

		Output:
		    shared gtgRightMatrix
	*/
	if batch.W <= cryptoParams.R {
		return PublicGtGActionExact(
			mpcObj,
			dataParams,
			batch,
			localGtG,
			rightMatrix,
		)
	}

	return ApplyGtG(
		mpcObj,
		heParams,
		dataParams,
		cryptoParams,
		batch,
		localGtG,
		rightMatrix,
		rhsCount,
	)
}

func ComputeGeneBatchKernelStatistics(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	gp []*mat.Dense,
	gv []*mat.Dense,
	x *mat.Dense,
	localGtx []*mat.Dense,
	localGtG []*mat.Dense,
	gvLocal []gvLocalGene,
	gvGene []gvGeneShares,
	omega mpc_core.RMat,
	weight []mpc_core.RVec,
	signedWeight []mpc_core.RVec,
	seed int64,
	observe func(stage string) func(),
) (geneBatchV, geneBatchS1, geneBatchS2, geneBatchS3 mpc_core.RVec) {
	/*
		Compute the phenotype-independent kernel statistics for one gene batch.

		Inputs in global gene:
		    gp, gv, weight, signedWeight

		Inputs in batch:
		    localGtx, localGtG, gvLocal, gvGene

		For each gene:
		    pooledGtx = (Transpose(Gp[A]) * X[A] + Transpose(Gp[B]) * X[B]) / N
		    pooledGtG = (Transpose(Gp[A]) * Gp[A] + Transpose(Gp[B]) * Gp[B]) / N

		    S1 = Trace(K/N)
		    S2 = Trace((K/N)²)
		    S3 = Trace((K/N)³)
		    V  = Burden variance / N

		Outputs V, S1, S2, S3 remain secret-shared and use batch order.
	*/
	geneCount := len(batch.GeneIndices)

	// 1. Share pooledGtx.
	// pooledGtx = (Transpose(Gp[A]) * X[A] + Transpose(Gp[B]) * X[B]) / N
	done := observe("kernel_inputs")
	pooledGtx := SharePooledGtx(mpcObj, dataParams, batch, localGtx)

	// 2. Share diagGtG.
	// diagGtG = diag(localGtG[A] + localGtG[B])
	packedLength := 0
	for _, geneIndex := range batch.GeneIndices {
		packedLength += dataParams.Genes[geneIndex].VariantCount
	}

	diagGtG := make([]mpc_core.RVec, geneCount)
	if packedLength > 0 {
		var packedLocal *mat.Dense

		if mpcObj.GetPid() != auxiliaryPartyID {
			packedLocal = mat.NewDense(1, packedLength, nil)
			packed := packedLocal.RawRowView(0)
			offset := 0

			for position, geneIndex := range batch.GeneIndices {
				variantCount := dataParams.Genes[geneIndex].VariantCount
				for variant := 0; variant < variantCount; variant++ {
					packed[offset+variant] = localGtG[position].At(variant, variant)
				}
				offset += variantCount
			}
		}

		sharedDiag := ShareSum(mpcObj, packedLocal, 1, packedLength)[0]
		offset := 0

		for position, geneIndex := range batch.GeneIndices {
			variantCount := dataParams.Genes[geneIndex].VariantCount
			if variantCount > 0 {
				diagGtG[position] = sharedDiag[offset : offset+variantCount].Copy()
			}
			offset += variantCount
		}
	}
	done()

	// 3. Compute S1, S2, S3 and reuse the first GtG action for gtgWeight.
	geneBatchS1, geneBatchS2, geneBatchS3, gtgWeight := ComputeKernelTraces(
		mpcObj, heParams, dataParams, cryptoParams, batch, gp, gv, x, gvLocal,
		localGtG, pooledGtx, diagGtG, omega, weight, signedWeight, seed, observe,
	)

	// 4. Compute the normalized Burden variance using gtgWeight = pooledGtG * signedWeight.
	done = observe("burden_variance")
	geneBatchV = ComputeBurdenVariance(
		mpcObj, batch, pooledGtx, signedWeight, gtgWeight, gvGene, omega,
	)
	done()

	return geneBatchV, geneBatchS1, geneBatchS2, geneBatchS3
}

func ComputeKernelTraces(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	gp []*mat.Dense,
	gv []*mat.Dense,
	x *mat.Dense,
	gvLocal []gvLocalGene,
	localGtG []*mat.Dense,
	pooledGtx []mpc_core.RMat,
	diagGtG []mpc_core.RVec,
	omega mpc_core.RMat,
	weight []mpc_core.RVec,
	signedWeight []mpc_core.RVec,
	seed int64,
	observe func(stage string) func(),
) (geneBatchS1, geneBatchS2, geneBatchS3 mpc_core.RVec, gtgWeight []mpc_core.RVec) {
	/*
		Compute the three kernel traces and the public GtG burden action.

		Inputs in global gene order:
		    gp, gv, weight, signedWeight

		Inputs in batch order:
		    gvLocal, localGtG, pooledGtx, diagGtG

		Outputs in batch order:
		    S1, S2, S3, gtgWeight
	*/
	geneCount := len(batch.GeneIndices)
	covariateCount := dataParams.C
	rtype := mpcObj.GetRType()
	dataBits := mpcObj.GetDataBits()
	fracBits := mpcObj.GetFracBits()

	// 1. Generate the public trace and private-correction probes.
	traceProbe := make([]*mat.Dense, geneCount)
	correctionProbe := make([]*mat.Dense, geneCount)
	probeScale := make([]float64, geneCount)

	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount

		traceProbe[position], correctionProbe[position], probeScale[position] =
			TraceProbe(variantCount, batch.W, cryptoParams.R, seed)
	}

	// 2. Prepare the fixed-public-shape private trace terms.
	privateTerms := PreparePrivateTraceTerms(
		mpcObj, dataParams, batch, gp, gv, x, gvLocal, correctionProbe,
	)

	// 3. Compute theta and weightSquared for each gene.
	//
	// theta         = pooledGtx * omega
	// weightSquared = weight .* weight
	batchWeight := make([]mpc_core.RVec, geneCount)
	batchSignedWeight := make([]mpc_core.RVec, geneCount)
	theta := make([]mpc_core.RMat, geneCount)
	weightSquared := make([]mpc_core.RVec, geneCount)

	for position, geneIndex := range batch.GeneIndices {
		batchWeight[position] = weight[geneIndex]
		batchSignedWeight[position] = signedWeight[geneIndex]

		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}

		theta[position] = mpcObj.SSMultMat(pooledGtx[position], omega)
		theta[position] = mpcObj.TruncMat(theta[position], dataBits, fracBits)

		weightSquared[position] = mpcObj.SSMultElemVec(batchWeight[position], batchWeight[position])
		weightSquared[position] = mpcObj.TruncVec(weightSquared[position], dataBits, fracBits)
	}

	// 3a. Compute tau1 before the public GtG actions.
	flattenMatrix := func(matrix mpc_core.RMat) mpc_core.RVec {
		rows, columns := matrix.Dims()
		values := mpc_core.InitRVec(rtype.Zero(), rows*columns)

		offset := 0
		for row := 0; row < rows; row++ {
			for column := 0; column < columns; column++ {
				values[offset] = matrix[row][column].Copy()
				offset++
			}
		}
		return values
	}

	tau1 := mpc_core.InitRVec(rtype.Zero(), geneCount)
	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}

		weightedTheta := scaleSharedRows(mpcObj, weightSquared[position], theta[position])
		diagTerm := sharedDot(mpcObj, weightSquared[position], diagGtG[position])
		projectionTerm := sharedDot(
			mpcObj, flattenMatrix(weightedTheta), flattenMatrix(pooledGtx[position]),
		)

		scaled := mpc_core.RVec{
			diagTerm.Sub(projectionTerm).Mul(rtype.FromFloat64(0.5, fracBits)),
		}
		tau1[position] = mpcObj.TruncVec(scaled, dataBits, fracBits)[0]
	}

	// 4. Build the first combined GtG right-hand side.
	// delta3Basis = ConcatColumns(c0GpGpProbe, c0GpX, theta)
	// weightedProbe = Diag(weight) * traceProbe
	// weightedBasisRight = Diag(weight)^2 * delta3Basis
	// firstRight = ConcatColumns(weightedProbe, weightedBasisRight, signedWeight)
	firstRight := make([]mpc_core.RMat, geneCount)
	weightedProbe := make([]mpc_core.RMat, geneCount)
	weightedBasisRight := make([]mpc_core.RMat, geneCount)

	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		if variantCount == 0 {
			continue
		}

		_, probeCount := traceProbe[position].Dims()
		basisColumnCount := probeCount + 2*covariateCount

		sharedProbe := mpc_core.InitRMat(rtype.Zero(), variantCount, probeCount)
		if mpcObj.GetPid() == mpcObj.GetHubPid() {
			for variant := 0; variant < variantCount; variant++ {
				for column := 0; column < probeCount; column++ {
					sharedProbe[variant][column] =
						rtype.FromFloat64(traceProbe[position].At(variant, column), fracBits)
				}
			}
		}

		weightedProbe[position] = scaleSharedRows(mpcObj, batchWeight[position], sharedProbe)

		delta3Basis := mpc_core.InitRMat(rtype.Zero(), variantCount, basisColumnCount)
		for variant := 0; variant < variantCount; variant++ {
			offset := 0

			for column := 0; column < probeCount; column++ {
				delta3Basis[variant][offset] = privateTerms[position].c0GpGpProbe[variant][column].Copy()
				offset++
			}

			for column := 0; column < covariateCount; column++ {
				delta3Basis[variant][offset] = privateTerms[position].c0GpX[variant][column].Copy()
				offset++
			}

			for column := 0; column < covariateCount; column++ {
				delta3Basis[variant][offset] = theta[position][variant][column].Copy()
				offset++
			}
		}

		basisRight := scaleSharedRows(mpcObj, batchWeight[position], delta3Basis)
		weightedBasisRight[position] = scaleSharedRows(mpcObj, batchWeight[position], basisRight)

		firstRight[position] =
			mpc_core.InitRMat(rtype.Zero(), variantCount, probeCount+basisColumnCount+1)
		for variant := 0; variant < variantCount; variant++ {
			offset := 0

			for column := 0; column < probeCount; column++ {
				firstRight[position][variant][offset] = weightedProbe[position][variant][column].Copy()
				offset++
			}

			for column := 0; column < basisColumnCount; column++ {
				firstRight[position][variant][offset] =
					weightedBasisRight[position][variant][column].Copy()
				offset++
			}

			firstRight[position][variant][offset] = batchSignedWeight[position][variant].Copy()
		}
	}

	// 5. Compute the first public GtG action.
	// firstGtgAction = pooledGtG * firstRight
	// Hutchinson rhsCount = 2*R + 2*C + 1.
	// Exact mode ignores rhsCount because its widths are gene-specific.
	done := observe("first_gtg_action")
	firstGtgAction := PublicGtGAction(
		mpcObj, heParams, dataParams, cryptoParams, batch, localGtG, firstRight,
		2*cryptoParams.R+2*covariateCount+1,
	)
	done()

	// 6. Split the first action and compute kProbe and basisAction.
	// kProbe = Kpp * traceProbe
	// basisAction = Kpp * Diag(weight) * delta3Basis
	// gtgWeight = (pooled GtG / N) * signedWeight
	kProbe := make([]mpc_core.RMat, geneCount)
	basisAction := make([]mpc_core.RMat, geneCount)
	gtgWeight = make([]mpc_core.RVec, geneCount)

	for position, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		if variantCount == 0 {
			continue
		}

		_, probeCount := traceProbe[position].Dims()
		basisColumnCount := probeCount + 2*covariateCount

		gtgProbe := mpc_core.InitRMat(rtype.Zero(), variantCount, probeCount)
		gtgBasis := mpc_core.InitRMat(rtype.Zero(), variantCount, basisColumnCount)
		gtgWeight[position] = mpc_core.InitRVec(rtype.Zero(), variantCount)

		for variant := 0; variant < variantCount; variant++ {
			for column := 0; column < probeCount; column++ {
				gtgProbe[variant][column] = firstGtgAction[position][variant][column].Copy()
			}

			for column := 0; column < basisColumnCount; column++ {
				gtgBasis[variant][column] = firstGtgAction[position][variant][probeCount+column].Copy()
			}

			gtgWeight[position][variant] =
				firstGtgAction[position][variant][probeCount+basisColumnCount].Copy()
		}

		kProbe[position] = ComputeKppRight(
			mpcObj, weightedProbe[position], gtgProbe,
			pooledGtx[position], theta[position], batchWeight[position],
		)
		basisAction[position] = ComputeKppRight(
			mpcObj, weightedBasisRight[position], gtgBasis,
			pooledGtx[position], theta[position], batchWeight[position],
		)
	}

	// 7. Compute the dependent second public GtG action.
	// weightedKProbe = Diag(weight) * kProbe
	// secondGtgAction = pooledGtG * weightedKProbe
	weightedKProbe := make([]mpc_core.RMat, geneCount)
	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}

		weightedKProbe[position] = scaleSharedRows(mpcObj, batchWeight[position], kProbe[position])
	}

	done = observe("second_gtg_action")
	secondGtgAction := PublicGtGAction(
		mpcObj, heParams, dataParams, cryptoParams, batch,
		localGtG, weightedKProbe, cryptoParams.R,
	)
	done()

	kSquaredProbe := make([]mpc_core.RMat, geneCount)
	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}

		kSquaredProbe[position] = ComputeKppRight(
			mpcObj, weightedKProbe[position], secondGtgAction[position],
			pooledGtx[position], theta[position], batchWeight[position],
		)
	}

	// 8. Compute tau2, tau3, and delta3Action.
	// tau2 = probeScale * Dot(kProbe, kProbe)
	// tau3 = probeScale * Dot(kProbe, kSquaredProbe)
	// delta3Action = 2 * Diag(weight) * basisAction
	tau2 := mpc_core.InitRVec(rtype.Zero(), geneCount)
	tau3 := mpc_core.InitRVec(rtype.Zero(), geneCount)
	delta3Action := make([]mpc_core.RMat, geneCount)

	for position, geneIndex := range batch.GeneIndices {
		if dataParams.Genes[geneIndex].VariantCount == 0 {
			continue
		}

		flatKProbe := flattenMatrix(kProbe[position])
		tau2Term := sharedDot(mpcObj, flatKProbe, flatKProbe)
		tau3Term := sharedDot(mpcObj, flatKProbe, flattenMatrix(kSquaredProbe[position]))

		scaledTau := mpc_core.RVec{
			tau2Term.Mul(rtype.FromFloat64(probeScale[position], fracBits)),
			tau3Term.Mul(rtype.FromFloat64(probeScale[position], fracBits)),
		}
		scaledTau = mpcObj.TruncVec(scaledTau, dataBits, fracBits)

		tau2[position] = scaledTau[0]
		tau3[position] = scaledTau[1]

		delta3Action[position] = scaleSharedRows(mpcObj, batchWeight[position], basisAction[position])
		delta3Action[position].MulScalar(rtype.FromInt(2))
	}

	// 9. Compute the three private trace corrections.
	delta1 := mpc_core.InitRVec(rtype.Zero(), geneCount)
	delta2 := mpc_core.InitRVec(rtype.Zero(), geneCount)
	delta3 := mpc_core.InitRVec(rtype.Zero(), geneCount)

	done = observe("private_trace_correction")
	for position := range batch.GeneIndices {
		delta1[position], delta2[position], delta3[position] = ComputeDelta(
			mpcObj, privateTerms[position], weightSquared[position], theta[position],
			omega, delta3Action[position], correctionProbe[position], probeScale[position],
		)
	}
	done()

	// 10. Assemble and return the final traces.
	//
	// S1 = tau1 + delta1
	// S2 = tau2 + delta2
	// S3 = tau3 + delta3
	geneBatchS1 = mpc_core.InitRVec(rtype.Zero(), geneCount)
	geneBatchS2 = mpc_core.InitRVec(rtype.Zero(), geneCount)
	geneBatchS3 = mpc_core.InitRVec(rtype.Zero(), geneCount)

	for position := range batch.GeneIndices {
		geneBatchS1[position] = tau1[position].Add(delta1[position])
		geneBatchS2[position] = tau2[position].Add(delta2[position])
		geneBatchS3[position] = tau3[position].Add(delta3[position])
	}

	return geneBatchS1, geneBatchS2, geneBatchS3, gtgWeight
}
