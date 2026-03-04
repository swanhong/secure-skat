package gwas

import (
	"fmt"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
)

// covAllOnes is a global flag to ensure the all-ones covariate is added only once.
// This is copied from assoc.go
var covAllOnes bool

// computeResidual performs QR and projects covariates out to get the null model residual
func (ast *AssocTest) computeResidual() crypto.CipherMatrix {
	cryptoParams := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()

	nrowsAll := make([]int, ast.general.config.NumMainParties+1)
	nrowsTotal := 0
	for p := 1; p <= ast.general.config.NumMainParties; p++ {
		nrowsAll[p] = ast.general.config.NumInds[p]
		nrowsTotal += nrowsAll[p]
	}
	nrowsTotalInv := 1.0 / float64(nrowsTotal)

	C := ast.inputCov
	if !covAllOnes {
		log.LLvl1("Adding an all-ones covariate for SKAT Null Model")
		if pid > 0 { // Party 0 doesn't encode data
			arr := make([]float64, nrowsAll[pid])
			for i := range arr {
				arr[i] = 1.0
			}
			pv, _ := crypto.EncodeFloatVector(cryptoParams, arr)
			C = append([]crypto.PlainVector{pv}, C...)
		} else {
			// Party 0 dummy vector
			C = append([]crypto.PlainVector{make(crypto.PlainVector, 0)}, C...)
		}
		covAllOnes = true
	} else {
		log.LLvl1("SKAT Warning: assumes the first covariate is all ones")
	}

	// Joint QR (NetDQRenc inside requires Party 0)
	// SKAT runs without PCA covariates, passing nil

	// Joint QR (NetDQRenc inside requires Party 0)
	// SKAT runs without PCA covariates, passing nil

	Q := ast.computeCombinedQV2(C, nil)

	SaveMatrixToFile(cryptoParams, mpcObj, Q, nrowsAll[pid], -1, ast.general.OutPath("Qcomb.txt"))

	// Project covariates out of y: ynew = (I - Q*Q')*y
	ymat := make(crypto.PlainMatrix, 1)
	if pid > 0 {
		ymat[0] = ast.pheno
	} else {
		ymat[0] = make(crypto.PlainVector, 0)
	}

	mmplainfn := func(cp *crypto.CryptoParams, a crypto.CipherVector, B crypto.PlainMatrix, j int) crypto.CipherVector {
		return crypto.CRescale(cp, crypto.CPMult(cp, a, B[j]))
	}

	var ynew crypto.CipherMatrix
	// DCMatMulAAtBPlain only aggregates across data parties. Party 0 can skip it or just get nil
	if pid > 0 {
		ynew = DCMatMulAAtBPlain(cryptoParams, mpcObj, Q, ymat, nrowsAll, 1, mmplainfn) // Level -2
		ynew[0] = crypto.CMultConstRescale(cryptoParams, ynew[0], nrowsTotalInv, true)

	} else {
		ynew = make(crypto.CipherMatrix, 1)
		ynew[0] = nil
	}

	// BootstrapVecAll requires Party 0
	ynew[0] = mpcObj.Network.BootstrapVecAll(cryptoParams, ynew[0])

	if pid > 0 {
		ynew[0] = crypto.CMultConst(cryptoParams, ynew[0], -1.0, true)
		ynew[0] = crypto.CPAdd(cryptoParams, ynew[0], ast.pheno)
	}

	SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[pid], -1, ast.general.OutPath("ynew.txt"))

	return ynew
}

// weightsCalculation computes SKAT weights (1/25)*(1-MAF)^24 via MPC
func (ast *AssocTest) weightsCalculation(dosageSum []float64, nsnps_block int) (mpc_core.RVec, mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()

	xSumRVec := mpc_core.InitRVec(rtype.Zero(), nsnps_block)

	if pid > 0 {
		for j := 0; j < nsnps_block; j++ {
			xSumRVec[j] = rtype.FromFloat64(dosageSum[j], 0) // exact integer sums
		}
	}

	// Calculate MAF securely
	// Total number of individuals N is perfectly public across all parties, so 2N is public.
	// We can compute p_j = xSum / (2N) securely by multiplying the secret-shared xSumRVec by the plaintext scalar 1/(2N).
	totalIndivs := ast.general.config.NumInds[1] + ast.general.config.NumInds[2]
	inv2N := rtype.FromFloat64(1.0/float64(2*totalIndivs), mpcObj.GetFracBits())

	p_j := xSumRVec.Copy()
	p_j.MulScalar(inv2N)
	// p_j inherently possesses mpcObj.GetFracBits() precision now (0 bits * frac bits = frac bits)

	one := rtype.FromFloat64(1.0, mpcObj.GetFracBits())
	var onesRVec mpc_core.RVec
	if pid == 1 {
		onesRVec = mpc_core.InitRVec(one, len(p_j))
	} else {
		onesRVec = mpc_core.InitRVec(rtype.Zero(), len(p_j))
	}
	onesRVec.Sub(p_j)
	w_term := onesRVec

	// Compute (1 - p_j)^24 via squaring
	w2 := mpcObj.SSMultElemVec(w_term, w_term)
	w2 = mpcObj.TruncVec(w2, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w4 := mpcObj.SSMultElemVec(w2, w2)
	w4 = mpcObj.TruncVec(w4, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w8 := mpcObj.SSMultElemVec(w4, w4)
	w8 = mpcObj.TruncVec(w8, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w16 := mpcObj.SSMultElemVec(w8, w8)
	w16 = mpcObj.TruncVec(w16, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w24 := mpcObj.SSMultElemVec(w16, w8)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	inv25 := rtype.FromFloat64(1.0/25.0, mpcObj.GetFracBits())
	w24.MulScalar(inv25)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	return p_j, w24
}

// ScoreCalculation calculates the final SKAT Score statistic iteratively
func (ast *AssocTest) ScoreCalculation(S_vec crypto.CipherVector, w_enc crypto.CipherVector) crypto.CipherVector {
	cryptoParams := ast.general.cps

	// Compute [s_j^2]
	S2 := crypto.CMult(cryptoParams, S_vec, S_vec)

	// Compute [w_j^2] and then [w_j^2 * s_j^2]
	w2 := crypto.CMult(cryptoParams, w_enc, w_enc)
	w2S2 := crypto.CMult(cryptoParams, w2, S2)

	// Sum across all SNPs in this block
	qBlock := crypto.InnerSumAll(cryptoParams, w2S2)

	return crypto.CipherVector{qBlock}
}

// Main SKAT computation function calling separated steps
func (ast *AssocTest) ComputeSKATStatistics() (qStat crypto.CipherVector, S_all crypto.CipherMatrix, outFilter []bool) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	cryptoParams := ast.general.cps

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting SKAT Phase 1: Null Model & Score Vector")

	// Step 1: Null Model Residuals
	ynew := ast.computeResidual()


	Y := crypto.CipherMatrix{ynew[0]}
	if pid == 0 {
		Y = crypto.CipherMatrix{nil}
	}

	var finalQStat crypto.CipherVector
	numBlocks := ast.general.config.GenoNumBlocks

	for block := 0; block < numBlocks; block++ {
		// Determine number of SNPs in this block
		var nsnps_block int
		if pid > 0 {
			isPgen := ast.general.IsPgen()
			blockSize := ast.general.genoBlockSizes[block]
			shift := uint64(0)
			for i := 0; i < block; i++ {
				shift += uint64(ast.general.genoBlockSizes[i])
			}

			if isPgen {
				if ast.general.gwasParams.snpFilt == nil {
					nsnps_block = blockSize
				} else {
					nsnps_block = SumBool(ast.general.gwasParams.snpFilt[shift : shift+uint64(blockSize)])
				}
			} else {
				nsnps_block = int(ast.general.genoBlocks[block].NumColsToKeep())
			}
		}

		if pid == 1 {
			mpcObj.Network.SendInt(nsnps_block, 0)
		} else if pid == 0 {
			nsnps_block = mpcObj.Network.ReceiveInt(1)
		}

		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Processing block %d/%d: Loading %d SNPs", block, numBlocks, nsnps_block))

		// Step 2: loadData & multi-party GenoBlockMult implicitly
		var S_block crypto.CipherMatrix
		var dosageSum []float64
		var filterBlock []bool

		if pid > 0 {
			S_block, dosageSum, _, filterBlock = ast.GenoBlockMult(block, Y, false)
			outFilter = append(outFilter, filterBlock...)
		} else {
			S_block = nil
			dosageSum = make([]float64, nsnps_block)
			outFilter = append(outFilter, make([]bool, nsnps_block)...)
		}

		// Step 3: weightsCalculation
		_, w_block_RVec := ast.weightsCalculation(dosageSum, nsnps_block)

		// Aggregate across parties
		S_block_aggr := mpcObj.Network.AggregateCMat(cryptoParams, S_block)
		S_block_aggr = mpcObj.Network.CollectiveBootstrapMat(cryptoParams, S_block_aggr, -1)

		if pid > 0 {
			S_vec := S_block_aggr[0]
			S_all = append(S_all, S_vec)

			w_enc := mpcObj.SSToCVec(cryptoParams, w_block_RVec)

			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{S_vec}, nsnps_block, -1, ast.general.OutPath(fmt.Sprintf("S_vec_block%d.txt", block)))
			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{w_enc}, nsnps_block, -1, ast.general.OutPath(fmt.Sprintf("w_enc_block%d.txt", block)))

			// Step 4: ScoreCalculation
			qBlockRes := ast.ScoreCalculation(S_vec, w_enc)

			SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{qBlockRes}, 1, -1, ast.general.OutPath(fmt.Sprintf("qBlock_block%d.txt", block)))

			if finalQStat == nil {
				finalQStat = qBlockRes
			} else {
				finalQStat = crypto.CAdd(cryptoParams, finalQStat, qBlockRes)
			}
		}
	}

	return finalQStat, S_all, outFilter
}
