package gwas

import (
	"fmt"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
)

// initSKAT skips the legacy association-test phenotype/covariate encodings, which SKAT never reads.
func (g *ProtocolInfo) initSKAT() *AssocTest { return &AssocTest{general: g} }

// ComputeSKATStatistics returns the whole-genome secure SKAT (Q) and Burden as
// 1-elem CipherVectors (Σ over all blocks, scaled by 1/(2σ̂²)). Caller collectively decrypts.
func (ast *AssocTest) ComputeSKATStatistics() (qStat, qBurden crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()

	null, X, y0 := ast.nullSetup()

	finalQSS := mpc_core.InitRVec(rtype.Zero(), 1)
	finalBurdenSS := mpc_core.InitRVec(rtype.Zero(), 1)
	for b := 0; b < ast.general.config.GenoNumBlocks; b++ {
		nsnps := ast.skatBlockNumSnps(b)
		if nsnps == 0 {
			continue
		}
		gl := ast.computeGeneLocal(b, nsnps, X, y0)
		qRawSS, bLinSS, _ := ast.blockStat(nsnps, null, gl)
		finalQSS.Add(qRawSS)
		finalBurdenSS.Add(bLinSS)
	}

	finalBurdenSS = mpcObj.SSSquareElemVec(finalBurdenSS)
	finalBurdenSS = mpcObj.TruncVec(finalBurdenSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	if scaleSS, ok := ast.general.rareVariantScaleShares(null.rssSS); ok {
		finalQSS = ast.general.scaleRareVariantShareStat(finalQSS, scaleSS)
		finalBurdenSS = ast.general.scaleRareVariantShareStat(finalBurdenSS, scaleSS)
	}
	return mpcObj.SSToCVec(cps, finalQSS), mpcObj.SSToCVec(cps, finalBurdenSS)
}

// ComputeSKATStatisticsPerBlock returns per-block Q and Burden (slot b = block b's
// statistic, scaled by the common 1/(2σ̂²)) — the per-gene secure SKAT statistics.
func (ast *AssocTest) ComputeSKATStatisticsPerBlock() (qPerBlock, burdenPerBlock crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()
	hub := mpcObj.GetPid() == mpcObj.GetHubPid()

	t0 := time.Now()
	null, X, y0 := ast.nullSetup()
	if hub {
		log.LLvl1(fmt.Sprintf(">>> [secure] null model (β̂, RSS) done (%.1fs)", time.Since(t0).Seconds()))
	}

	nB := ast.general.config.GenoNumBlocks
	qBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)
	bLinBlockSS := mpc_core.InitRVec(rtype.Zero(), nB)
	done := 0
	for b := 0; b < nB; b++ {
		nsnps := ast.skatBlockNumSnps(b)
		if nsnps == 0 {
			continue
		}
		tb := time.Now()
		gl := ast.computeGeneLocal(b, nsnps, X, y0)
		q, bl, _ := ast.blockStat(nsnps, null, gl)
		qBlockSS[b] = q[0]
		bLinBlockSS[b] = bl[0]
		done++
		if hub {
			log.LLvl1(fmt.Sprintf(">>> [secure] block %d/%d done (%d snps, %.1fs; elapsed %.1fs)",
				done, nB, nsnps, time.Since(tb).Seconds(), time.Since(t0).Seconds()))
		}
	}

	burdenSqSS := mpcObj.SSSquareElemVec(bLinBlockSS)
	burdenSqSS = mpcObj.TruncVec(burdenSqSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	if scaleSS, ok := ast.general.rareVariantScaleShares(null.rssSS); ok {
		qBlockSS = ast.general.scaleRareVariantShareStat(qBlockSS, scaleSS)
		burdenSqSS = ast.general.scaleRareVariantShareStat(burdenSqSS, scaleSS)
	}
	return mpcObj.SSToCVec(cps, qBlockSS), mpcObj.SSToCVec(cps, burdenSqSS)
}

func (ast *AssocTest) skatNumInds() []int {
	filtNumInds := ast.general.gwasParams.FiltNumInds()
	if len(filtNumInds) == ast.general.config.NumMainParties+1 {
		return filtNumInds
	}

	return ast.general.config.NumInds
}

func (ast *AssocTest) skatTotalNumInds() int {
	nrows := ast.skatNumInds()
	return Sum(nrows[1 : ast.general.config.NumMainParties+1])
}
