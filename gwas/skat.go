package gwas

import (
	"fmt"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
)

func (ast *AssocTest) skatNumInds() []int {
	filtNumInds := ast.general.gwasParams.FiltNumInds()
	if len(filtNumInds) == ast.general.config.NumMainParties+1 {
		return filtNumInds
	}

	return ast.general.config.NumInds
}

func (ast *AssocTest) skatTotalNumInds() int {
	nrows := ast.skatNumInds()
	total := 0
	for p := 1; p <= ast.general.config.NumMainParties; p++ {
		total += nrows[p]
	}
	return total
}

func (ast *AssocTest) zeroPlainVectorForNonDataParty() crypto.PlainVector {
	zeros := make([]float64, ast.general.cps.GetSlots())
	pv, _ := crypto.EncodeFloatVector(ast.general.cps, zeros)
	return pv
}

// computeResidual performs QR and projects covariates out to get the null model residual
func (ast *AssocTest) computeResidual() crypto.CipherMatrix {
	cryptoParams := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()

	debug := ast.general.config.Debug

	nrowsAll := ast.skatNumInds()
	nrowsTotal := ast.skatTotalNumInds()
	nrowsTotalInv := 1.0 / float64(nrowsTotal)

	// covAllOnes specifies whether the input covariates already contain an all-ones column.
	// SKAT needs an intercept; only prepend one when the input does not already carry it.
	covAllOnes := ast.general.config.CovAllOnes
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
			// Keep one zero-filled ciphertext so distributed QR sees aligned local shapes.
			C = append([]crypto.PlainVector{ast.zeroPlainVectorForNonDataParty()}, C...)
		}
	} else {
		log.LLvl1("CovAllOnes=true: assuming first covariate is the all-ones intercept")
	}

	// Joint QR (NetDQRenc inside requires Party 0)
	// SKAT runs without PCA covariates, passing nil
	Q := ast.computeCombinedQV2(C, nil)

	// In SKAT, the first covariate is ALWAYS the identically precise intercept.
	// Override Q[0] with 1.0 to eliminate NetDQRenc Goldschmidt approximation noise.
	if pid > 0 {
		slots := cryptoParams.GetSlots()
		ct := crypto.CZeros(cryptoParams, 1)[0]
		ct = crypto.AddConst(cryptoParams, ct, 1.0)

		QFirst := make(crypto.CipherVector, ((nrowsAll[pid]-1)/slots)+1)
		for i := range QFirst {
			nElem := slots
			if i == len(QFirst)-1 {
				nElem = nrowsAll[pid] - (len(QFirst)-1)*slots
			}
			QFirst[i] = crypto.MaskTrunc(cryptoParams, ct, nElem)
		}
		Q[0] = QFirst
		Q, _ = crypto.FlattenLevels(cryptoParams, Q)
	}

	// Save vertically partitioned arrays by iterating individually to avoid
	// MPC extracting ciphertexts composed of independent data halves (-1).
	if debug {
		for p := 1; p <= ast.general.GetConfig().NumMainParties; p++ {
			SaveMatrixToFile(cryptoParams, mpcObj, Q, nrowsAll[p], p, ast.general.OutPath("Qcomb.txt"))
		}
	}

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

	// WARNING: BootstrapVecAll also resolves to DummyBootstrapping which corrupts
	// ciphertexts if run locally without AggregateSk. We bypass it as PN14QP438
	// has ample depth.
	// ynew[0] = mpcObj.Network.BootstrapVecAll(cryptoParams, ynew[0])

	if pid > 0 {
		ynew[0] = crypto.CMultConst(cryptoParams, ynew[0], -1.0, true)
		ynew[0] = crypto.CPAdd(cryptoParams, ynew[0], ast.pheno)
	}

	if debug {
		for p := 1; p <= ast.general.GetConfig().NumMainParties; p++ {
			SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[p], p, ast.general.OutPath("ynew.txt"))
		}
	}

	return ynew
}

// weightsCalculation computes the shared beta-density term used by the SKAT score.
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
	totalIndivs := ast.skatTotalNumInds()
	inv2N := rtype.FromFloat64(1.0/float64(2*totalIndivs), mpcObj.GetFracBits())

	p_j := xSumRVec.Copy()
	p_j.MulScalar(inv2N)
	// p_j inherently possesses mpcObj.GetFracBits() precision now (0 bits * frac bits = frac bits)

	one := rtype.FromFloat64(1.0, mpcObj.GetFracBits())
	var onesRVec mpc_core.RVec
	if pid == mpcObj.GetHubPid() {
		onesRVec = mpc_core.InitRVec(one, len(p_j))
	} else {
		onesRVec = mpc_core.InitRVec(rtype.Zero(), len(p_j))
	}
	onesRVec.Sub(p_j)
	pVec := onesRVec // p = 1 - p_j

	// SKAT weight base is 1-MAF = max(p_j, 1-p_j), invariant to allele orientation.
	// betaBase = p_j + [p_j < p]*(p - p_j).
	useBoolean := mpcObj.GetBooleanShareFlag()
	majorSelect := mpcObj.LessThan(p_j, pVec, useBoolean)
	majorSelect.MulScalar(one)

	betaBase := p_j.Copy()
	majorDelta := pVec.Copy()
	majorDelta.Sub(p_j)
	majorDelta = mpcObj.SSMultElemVec(majorDelta, majorSelect)
	majorDelta = mpcObj.TruncVec(majorDelta, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	betaBase.Add(majorDelta)

	// Compute betaBase^24 = (1 - MAF)^24 via squaring
	w2 := mpcObj.SSMultElemVec(betaBase, betaBase)
	w2 = mpcObj.TruncVec(w2, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w4 := mpcObj.SSMultElemVec(w2, w2)
	w4 = mpcObj.TruncVec(w4, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w8 := mpcObj.SSMultElemVec(w4, w4)
	w8 = mpcObj.TruncVec(w8, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w16 := mpcObj.SSMultElemVec(w8, w8)
	w16 = mpcObj.TruncVec(w16, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	w24 := mpcObj.SSMultElemVec(w16, w8)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	// In standard SKAT, the weight is dbeta(MAF, 1, 25)
	// The beta density is f(x) = x^(a-1)*(1-x)^(b-1) / B(a,b)
	// B(1, 25) = 1/25. So f(x) = 25 * (1-x)^24
	betaConst := rtype.FromFloat64(25.0, mpcObj.GetFracBits())
	w24.MulScalar(betaConst)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	return p_j, w24
}

// ScoreCalculation calculates the final SKAT Score statistic iteratively
func (ast *AssocTest) ScoreCalculation(S_vec crypto.CipherVector, w_enc crypto.CipherVector) (crypto.CipherVector, crypto.CipherVector) {
	cryptoParams := ast.general.cps

	// Compute [s_j^2]
	S2 := crypto.CMult(cryptoParams, S_vec, S_vec)

	// Compute [w_j^2] and then [w_j^2 * s_j^2]
	w2 := crypto.CMult(cryptoParams, w_enc, w_enc)
	w2S2 := crypto.CMult(cryptoParams, w2, S2)

	// Sum across all SNPs in this block
	qSkatBlock := crypto.InnerSumAll(cryptoParams, w2S2)

	// Compute [w_j * s_j] for Burden
	wS := crypto.CMult(cryptoParams, w_enc, S_vec)
	qBurdenBlock := crypto.InnerSumAll(cryptoParams, wS)

	return crypto.CipherVector{qSkatBlock}, crypto.CipherVector{qBurdenBlock}
}

// Main SKAT computation function calling separated steps
func (ast *AssocTest) ComputeSKATStatistics() (qStat crypto.CipherVector, qBurden crypto.CipherVector, S_all crypto.CipherMatrix, outFilter []bool) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	cryptoParams := ast.general.cps
	debug := ast.general.config.Debug

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting SKAT Phase 1: Null Model & Score Vector")

	// Step 1: Null Model Residuals
	ynew := ast.computeResidual()

	Y := crypto.CipherMatrix{ynew[0]}
	if pid == 0 {
		Y = crypto.CipherMatrix{nil}
	}

	var finalQStat crypto.CipherVector
	var finalBurdenStat crypto.CipherVector
	numBlocks := ast.general.config.GenoNumBlocks
	hubPid := mpcObj.GetHubPid()

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

		if pid == hubPid {
			mpcObj.Network.SendInt(nsnps_block, 0)
		} else if pid == 0 {
			nsnps_block = mpcObj.Network.ReceiveInt(hubPid)
		}

		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Processing block %d/%d: Loading %d SNPs", block, numBlocks, nsnps_block))

		if nsnps_block == 0 {
			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Block %d/%d skipped (empty)", block, numBlocks))
			continue
		}

		// Step 2: loadData & multi-party GenoBlockMult implicitly
		var S_block crypto.CipherMatrix
		var dosageSum []float64
		var filterBlock []bool

		if pid > 0 {
			S_block, dosageSum, _, filterBlock = ast.GenoBlockMultSKAT(block, Y, false)
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

			if debug {
				SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{S_vec}, nsnps_block, -1, ast.general.OutPath(fmt.Sprintf("S_vec_block%d.txt", block)))
				SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{w_enc}, nsnps_block, -1, ast.general.OutPath(fmt.Sprintf("w_enc_block%d.txt", block)))
			}

			// Step 4: ScoreCalculation
			qBlockRes, qBurdenBlockRes := ast.ScoreCalculation(S_vec, w_enc)

			if debug {
				SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{qBlockRes}, 1, -1, ast.general.OutPath(fmt.Sprintf("qBlock_block%d.txt", block)))
				SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{qBurdenBlockRes}, 1, -1, ast.general.OutPath(fmt.Sprintf("qBurdenBlock_block%d.txt", block)))
			}

			if finalQStat == nil {
				finalQStat = qBlockRes
				finalBurdenStat = qBurdenBlockRes
			} else {
				finalQStat = crypto.CAdd(cryptoParams, finalQStat, qBlockRes)
				finalBurdenStat = crypto.CAdd(cryptoParams, finalBurdenStat, qBurdenBlockRes)
			}
		}
	}

	if pid > 0 && finalBurdenStat != nil {
		// Secure SKAT computes Burden as the square of the sum of block scores
		finalBurdenStat = crypto.CMult(cryptoParams, finalBurdenStat, finalBurdenStat)
	}

	return finalQStat, finalBurdenStat, S_all, outFilter
}
