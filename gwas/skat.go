package gwas

import (
	"fmt"
	"os"
	"strings"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"go.dedis.ch/onet/v3/log"
)

type SKATBlockData struct {
	NumSnps   int
	ScoreVec  crypto.CipherVector
	DosageSum []float64
	PEnc      crypto.CipherVector
	PBarEnc   crypto.CipherVector
	WeightEnc crypto.CipherVector
}

func (ast *AssocTest) skatNumInds() []int {
	filtNumInds := ast.general.gwasParams.FiltNumInds()
	if len(filtNumInds) == ast.general.config.NumMainParties+1 {
		return filtNumInds
	}

	return ast.general.config.NumInds
}

func (ast *AssocTest) skatTotalNumInds() int {
	total := 0
	for p := 1; p <= ast.general.config.NumMainParties; p++ {
		total += ast.skatNumInds()[p]
	}
	return total
}

func (ast *AssocTest) zeroPlainVectorForNonDataParty() crypto.PlainVector {
	zeros := make([]float64, ast.general.cps.GetSlots())
	pv, _ := crypto.EncodeFloatVector(ast.general.cps, zeros)
	return pv
}

func (ast *AssocTest) computeEncryptedResidualRSS(ynew crypto.CipherMatrix) crypto.CipherVector {
	cryptoParams := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()

	if pid == 0 {
		return nil
	}

	var rssLocal *rlwe.Ciphertext
	if pid > 0 && len(ynew) > 0 && ynew[0] != nil {
		ynewSq := crypto.CMult(cryptoParams, ynew[0], ynew[0])
		rssLocal = crypto.InnerSumAll(cryptoParams, ynewSq)
	} else {
		rssLocal = crypto.CZeros(cryptoParams, 1)[0]
	}

	return mpcObj.Network.AggregateCVec(cryptoParams, crypto.CipherVector{rssLocal})
}

// computeResidual performs QR and projects covariates out to get the null model residual
func (ast *AssocTest) computeResidual() crypto.CipherMatrix {
	cryptoParams := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()

	nrowsAll := ast.skatNumInds()
	nrowsTotal := 0
	for p := 1; p <= ast.general.config.NumMainParties; p++ {
		nrowsTotal += nrowsAll[p]
	}
	nrowsTotalInv := 1.0 / float64(nrowsTotal)

	// covAllOnes specifies whether the input covariates already contain an all-ones column.
	// For SKAT, we assume it does NOT locally, so we explicitly prepend it.
	C := ast.inputCov

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

	if ast.general.config.Debug {
		// Save vertically partitioned arrays by iterating individually to avoid
		// MPC extracting ciphertexts composed of independent data halves (-1).
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
		return crypto.CPMult(cp, a, B[j])
	}

	var ynew crypto.CipherMatrix
	var QTY crypto.CipherMatrix
	// DCMatMulAAtBPlain only aggregates across data parties. Party 0 can skip it or just get nil
	if pid > 0 {
		var err error
		ynew, QTY, err = DCMatMulAAtBPlainWithIntmd(cryptoParams, mpcObj, Q, ymat, nrowsAll, 1, mmplainfn) // Level -2
		if err != nil {
			log.Lvl1("Error in DCMatMulAAtBPlainWithIntmd: ", err)
		}

		if ast.general.config.Debug {
			for p := 1; p <= ast.general.GetConfig().NumMainParties; p++ {
				if QTY != nil {
					// Manually decrypt the 0-th slot of each ciphertext in QTY[0]
					ptList := mpcObj.Network.CollectiveDecryptVec(cryptoParams, QTY[0], p)

					if pid == p {
						f, err := os.Create(ast.general.OutPath("qty.txt"))
						if err == nil {
							var vals []string
							for _, pt := range ptList {
								v := crypto.DecodeFloatVector(cryptoParams, crypto.PlainVector{pt})[0]
								vals = append(vals, fmt.Sprintf("%.6e", v))
							}
							f.WriteString(strings.Join(vals, ",") + "\n")
							f.Close()
						}
					}
				}
			}
			// Save y_proj before rescaling/subtraction
			for p := 1; p <= ast.general.GetConfig().NumMainParties; p++ {
				SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[p], p, ast.general.OutPath("y_proj_raw.txt"))
			}
		}

		ynew[0] = crypto.CMultConstRescale(cryptoParams, ynew[0], nrowsTotalInv, true)

		if ast.general.config.Debug {
			// Save y_proj AFTER rescaling but before subtracting from pheno
			for p := 1; p <= ast.general.GetConfig().NumMainParties; p++ {
				SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[p], p, ast.general.OutPath("y_proj_rescaled.txt"))
			}
		}

	} else {
		ynew = make(crypto.CipherMatrix, 1)
		ynew[0] = nil
	}

	// WARNING: BootstrapVecAll also resolves to DummyBootstrapping which corrupts
	// ciphertexts if run locally without AggregateSk. We bypass it as PN14QP438
	// has ample depth.
	// ynew[0] = mpcObj.Network.BootstrapVecAll(cryptoParams, ynew[0])

	if pid > 0 {
		ynew[0] = crypto.CPSubOther(cryptoParams, ast.pheno, ynew[0])
	}

	if ast.general.config.Debug {
		for p := 1; p <= ast.general.GetConfig().NumMainParties; p++ {
			SaveMatrixToFile(cryptoParams, mpcObj, ynew, nrowsAll[p], p, ast.general.OutPath("ynew.txt"))
		}
	}

	return ynew
}

// weightsCalculation computes per-variant p_bar = dosageSum/(2N), p = 1-p_bar,
// then uses maf = min(p, p_bar) to form beta weights Beta(maf; 1, 25)
// = 25 * (1-maf)^24.
//
// This matches the current SKAT package notation, where the weighted linear
// kernel is K = G W W G' and the user-supplied "weights" correspond to the
// diagonal entries of W rather than the older paper's W^2 notation.
func (ast *AssocTest) weightsCalculation(dosageSum []float64, nsnps_block int) (mpc_core.RVec, mpc_core.RVec, mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()
	useBoolean := mpcObj.GetBooleanShareFlag()

	dosageSumRVec := mpc_core.InitRVec(rtype.Zero(), nsnps_block)

	if pid > 0 {
		for j := 0; j < nsnps_block; j++ {
			dosageSumRVec[j] = rtype.FromFloat64(dosageSum[j], 0) // exact integer sums
		}
	}

	// The secure genotype orientation already yields dosage_bar = 2N - dosage.
	// Therefore dosageSum/(2N) is p_bar relative to the plain reference, and
	// p = 1-p_bar must be derived explicitly in the shared domain.
	totalIndivs := ast.skatTotalNumInds()
	inv2N := rtype.FromFloat64(1.0/float64(2*totalIndivs), mpcObj.GetFracBits())

	pBarVec := dosageSumRVec.Copy()
	pBarVec.MulScalar(inv2N)

	oneFrac := rtype.FromFloat64(1.0, mpcObj.GetFracBits())
	var oneShared mpc_core.RVec
	if pid == mpcObj.GetHubPid() {
		oneShared = mpc_core.InitRVec(oneFrac, nsnps_block)
	} else {
		oneShared = mpc_core.InitRVec(rtype.Zero(), nsnps_block)
	}
	oneShared.Sub(pBarVec)
	pVec := oneShared

	// Use MAF = min(p, 1-p) so the weighting is invariant to allele orientation.
	// Beta(MAF; 1, 25) = 25 * (1-MAF)^24, so the exponent base is max(p, 1-p).
	majorSelect := mpcObj.LessThan(pBarVec, pVec, useBoolean)
	majorSelect.MulScalar(oneFrac)

	betaBase := pBarVec.Copy()
	majorDelta := pVec.Copy()
	majorDelta.Sub(pBarVec)
	majorDelta = mpcObj.SSMultElemVec(majorDelta, majorSelect)
	majorDelta = mpcObj.TruncVec(majorDelta, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	betaBase.Add(majorDelta)

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

	// In the current SKAT package notation, beta(MAF; 1, 25) = 25 * (1-MAF)^24
	// is the weight vector supplied to SKAT(..., weights = ...).
	// Since maf = min(p, 1-p), this equals 25 * max(p, 1-p)^24.
	betaConst := rtype.FromFloat64(25.0, mpcObj.GetFracBits())
	w24.MulScalar(betaConst)
	w24 = mpcObj.TruncVec(w24, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	return pVec, pBarVec, w24
}

// ScoreCalculation calculates the final SKAT Score statistic iteratively.
// It also returns key intermediate vectors to simplify debugging of score aggregation.
func (ast *AssocTest) ScoreCalculation(scoreVec crypto.CipherVector, weightEnc crypto.CipherVector) (
	crypto.CipherVector,
	crypto.CipherVector,
	crypto.CipherVector,
	crypto.CipherVector,
	crypto.CipherVector,
	crypto.CipherVector,
) {
	cryptoParams := ast.general.cps

	// Compute [s_j^2]
	S2 := crypto.CMult(cryptoParams, scoreVec, scoreVec)

	// Compute [w_j^2] and then [w_j^2 * s_j^2].
	// This matches SKAT::SKAT()'s current weighted linear kernel notation.
	w2 := crypto.CMult(cryptoParams, weightEnc, weightEnc)
	w2S2 := crypto.CMult(cryptoParams, w2, S2)

	// Sum across all SNPs in this block
	qSkatBlock := crypto.InnerSumAll(cryptoParams, w2S2)

	// Compute [w_j * s_j] for Burden
	wS := crypto.CMult(cryptoParams, weightEnc, scoreVec)
	qBurdenBlock := crypto.InnerSumAll(cryptoParams, wS)

	return crypto.CipherVector{qSkatBlock}, crypto.CipherVector{qBurdenBlock}, S2, w2, w2S2, wS
}

func (ast *AssocTest) ComputeSKATStep1Residuals() crypto.CipherMatrix {
	log.LLvl1(time.Now().Format(time.RFC3339), "SKAT Step 1/4: Null model residuals")
	return ast.computeResidual()
}

func (ast *AssocTest) skatBlockNumSnps(block int) int {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	hubPid := mpcObj.GetHubPid()
	if pid == 0 {
		return mpcObj.Network.ReceiveInt(hubPid)
	}

	var nsnpsBlock int
	isPgen := ast.general.IsPgen()
	blockSize := ast.general.genoBlockSizes[block]
	shift := uint64(0)
	for i := 0; i < block; i++ {
		shift += uint64(ast.general.genoBlockSizes[i])
	}

	if isPgen {
		if ast.general.gwasParams.snpFilt == nil {
			nsnpsBlock = blockSize
		} else {
			nsnpsBlock = SumBool(ast.general.gwasParams.snpFilt[shift : shift+uint64(blockSize)])
		}
	} else {
		nsnpsBlock = int(ast.general.genoBlocks[block].NumColsToKeep())
	}

	if pid == hubPid {
		mpcObj.Network.SendInt(nsnpsBlock, 0)
	}
	return nsnpsBlock
}

func (ast *AssocTest) ComputeSKATStep2LoadBlockScore(block int, Y crypto.CipherMatrix) SKATBlockData {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	cryptoParams := ast.general.cps
	numBlocks := ast.general.config.GenoNumBlocks

	nsnpsBlock := ast.skatBlockNumSnps(block)
	log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("SKAT Step 2/4: block %d/%d score loading (%d SNPs)", block+1, numBlocks, nsnpsBlock))

	var SBlock crypto.CipherMatrix
	var dosageSum []float64
	if pid > 0 {
		SBlock, dosageSum, _, _ = ast.GenoBlockMultSKAT(block, Y, false)
	} else {
		SBlock = nil
		dosageSum = make([]float64, nsnpsBlock)
	}

	SBlockAggr := mpcObj.Network.AggregateCMat(cryptoParams, SBlock)
	if pid > 0 && len(SBlockAggr) > 0 && len(SBlockAggr[0]) > 0 {
		if mpcObj.Network.CanCollectiveBootstrap(cryptoParams, SBlockAggr[0][0].Level()) {
			SBlockAggr = mpcObj.Network.CollectiveBootstrapMat(cryptoParams, SBlockAggr, -1)
		} else {
			log.LLvl1(time.Now().Format(time.RFC3339), "SKAT Step 2/4: skipping score bootstrap at level", SBlockAggr[0][0].Level())
		}
	}

	blockData := SKATBlockData{
		NumSnps:   nsnpsBlock,
		DosageSum: dosageSum,
	}
	if pid > 0 {
		blockData.ScoreVec = SBlockAggr[0]
	}
	return blockData
}

func (ast *AssocTest) ComputeSKATStep3BlockWeights(blockData *SKATBlockData) {
	mpcObj := ast.general.mpcObj[0]
	log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("SKAT Step 3/4: block weight calculation (%d SNPs)", blockData.NumSnps))
	pBlockRVec, pBarBlockRVec, weightBlockRVec := ast.weightsCalculation(blockData.DosageSum, blockData.NumSnps)

	if mpcObj.GetPid() == 0 {
		return
	}

	if ast.general.config.Debug {
		blockData.PEnc = mpcObj.SSToCVec(ast.general.cps, pBlockRVec)
		blockData.PBarEnc = mpcObj.SSToCVec(ast.general.cps, pBarBlockRVec)
	}
	blockData.WeightEnc = mpcObj.SSToCVec(ast.general.cps, weightBlockRVec)
}

func (ast *AssocTest) saveSKATStep3Outputs(block int, blockData SKATBlockData) {
	cryptoParams := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	if mpcObj.GetPid() == 0 || !ast.general.config.Debug {
		return
	}

	SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{blockData.PEnc}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("p_block%d.txt", block)))
	SaveMatrixComplexPartsToFile(
		cryptoParams,
		mpcObj,
		crypto.CipherMatrix{blockData.PEnc},
		blockData.NumSnps,
		-1,
		ast.general.OutPath(fmt.Sprintf("p_block%d_real.txt", block)),
		ast.general.OutPath(fmt.Sprintf("p_block%d_imag.txt", block)),
	)
	SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{blockData.PBarEnc}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("p_bar_block%d.txt", block)))
	SaveMatrixComplexPartsToFile(
		cryptoParams,
		mpcObj,
		crypto.CipherMatrix{blockData.PBarEnc},
		blockData.NumSnps,
		-1,
		ast.general.OutPath(fmt.Sprintf("p_bar_block%d_real.txt", block)),
		ast.general.OutPath(fmt.Sprintf("p_bar_block%d_imag.txt", block)),
	)
	SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{blockData.ScoreVec}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("S_vec_block%d.txt", block)))
	SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{blockData.WeightEnc}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("w_enc_block%d.txt", block)))
	SaveMatrixComplexPartsToFile(
		cryptoParams,
		mpcObj,
		crypto.CipherMatrix{blockData.WeightEnc},
		blockData.NumSnps,
		-1,
		ast.general.OutPath(fmt.Sprintf("w_enc_block%d_real.txt", block)),
		ast.general.OutPath(fmt.Sprintf("w_enc_block%d_imag.txt", block)),
	)
}

func (ast *AssocTest) ComputeSKATStep4BlockStatistics(block int, blockData SKATBlockData) (crypto.CipherVector, crypto.CipherVector) {
	cryptoParams := ast.general.cps
	mpcObj := ast.general.mpcObj[0]
	if mpcObj.GetPid() == 0 {
		return nil, nil
	}

	log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("SKAT Step 4/4: block statistic aggregation (%d SNPs)", blockData.NumSnps))
	qBlockRes, qBurdenBlockRes, S2, w2, w2S2, wS := ast.ScoreCalculation(blockData.ScoreVec, blockData.WeightEnc)

	SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{qBlockRes}, 1, -1, ast.general.OutPath(fmt.Sprintf("qBlock_block%d.txt", block)))
	SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{qBurdenBlockRes}, 1, -1, ast.general.OutPath(fmt.Sprintf("qBurdenBlock_block%d.txt", block)))

	if ast.general.config.Debug {
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{S2}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("S2_block%d.txt", block)))
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{w2}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("w2_block%d.txt", block)))
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{w2S2}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("w2S2_block%d.txt", block)))
		SaveMatrixToFile(cryptoParams, mpcObj, crypto.CipherMatrix{wS}, blockData.NumSnps, -1, ast.general.OutPath(fmt.Sprintf("wS_block%d.txt", block)))
		SaveMatrixComplexPartsToFile(
			cryptoParams,
			mpcObj,
			crypto.CipherMatrix{w2},
			blockData.NumSnps,
			-1,
			ast.general.OutPath(fmt.Sprintf("w2_block%d_real.txt", block)),
			ast.general.OutPath(fmt.Sprintf("w2_block%d_imag.txt", block)),
		)
		SaveMatrixComplexPartsToFile(
			cryptoParams,
			mpcObj,
			crypto.CipherMatrix{w2S2},
			blockData.NumSnps,
			-1,
			ast.general.OutPath(fmt.Sprintf("w2S2_block%d_real.txt", block)),
			ast.general.OutPath(fmt.Sprintf("w2S2_block%d_imag.txt", block)),
		)
	}

	return qBlockRes, qBurdenBlockRes
}

func (ast *AssocTest) scalarCiphertextToShares(stat crypto.CipherVector) mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	if stat == nil || len(stat) == 0 || stat[0] == nil {
		return mpc_core.InitRVec(rtype.Zero(), 1)
	}
	return mpcObj.CiphertextToSS(ast.general.cps, rtype, stat[0], -1, 1)
}

// ComputeSKATStatistics returns the final secure SKAT and Burden statistics.
func (ast *AssocTest) ComputeSKATStatistics() (qStat crypto.CipherVector, qBurden crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	cryptoParams := ast.general.cps
	rtype := mpcObj.GetRType()

	ynew := ast.ComputeSKATStep1Residuals()
	nullRSS := ast.computeEncryptedResidualRSS(ynew)
	Y := crypto.CipherMatrix{ynew[0]}
	if pid == 0 {
		Y = crypto.CipherMatrix{nil}
	}

	finalQSS := mpc_core.InitRVec(rtype.Zero(), 1)
	finalBurdenSS := mpc_core.InitRVec(rtype.Zero(), 1)
	numBlocks := ast.general.config.GenoNumBlocks

	for block := 0; block < numBlocks; block++ {
		blockData := ast.ComputeSKATStep2LoadBlockScore(block, Y)

		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("SKAT Progress: block %d/%d (%.1f%%)", block+1, numBlocks, 100.0*float64(block+1)/float64(numBlocks)))

		ast.ComputeSKATStep3BlockWeights(&blockData)
		ast.saveSKATStep3Outputs(block, blockData)

		qBlockRes, qBurdenBlockRes := ast.ComputeSKATStep4BlockStatistics(block, blockData)
		finalQSS.Add(ast.scalarCiphertextToShares(qBlockRes))
		finalBurdenSS.Add(ast.scalarCiphertextToShares(qBurdenBlockRes))
	}

	finalBurdenSS = mpcObj.SSSquareElemVec(finalBurdenSS)
	finalBurdenSS = mpcObj.TruncVec(finalBurdenSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())

	scaleSS, scaleOK := ast.general.rareVariantScaleShares(nullRSS)
	if scaleOK {
		finalQSS = ast.general.scaleRareVariantShareStat(finalQSS, scaleSS)
		finalBurdenSS = ast.general.scaleRareVariantShareStat(finalBurdenSS, scaleSS)
	}

	finalQStat := mpcObj.SSToCVec(cryptoParams, finalQSS)
	finalBurdenStat := mpcObj.SSToCVec(cryptoParams, finalBurdenSS)

	return finalQStat, finalBurdenStat
}
