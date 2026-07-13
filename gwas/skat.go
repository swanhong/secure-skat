package gwas

import (
	"fmt"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"go.dedis.ch/onet/v3/log"
)

// ComputeSKATStatistics returns the whole-genome secure SKAT (Q) and Burden as
// 1-elem CipherVectors (Σ over all blocks, scaled by 1/(2σ̂²)). Caller collectively decrypts.
func (ast *AssocTest) ComputeSKATStatistics() (qStat, qBurden crypto.CipherVector) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	rtype := mpcObj.GetRType()

	null, nullRSS, X, y0 := ast.nullSetup()

	finalQSS := mpc_core.InitRVec(rtype.Zero(), 1)
	finalBurdenSS := mpc_core.InitRVec(rtype.Zero(), 1)
	for b := 0; b < ast.general.config.GenoNumBlocks; b++ {
		nsnps := ast.skatBlockNumSnps(b)
		if nsnps == 0 {
			continue
		}
		qRawSS, bLinSS, _ := ast.blockStat(b, nsnps, null, X, y0, nil)
		finalQSS.Add(qRawSS)
		finalBurdenSS.Add(bLinSS)
	}

	finalBurdenSS = mpcObj.SSSquareElemVec(finalBurdenSS)
	finalBurdenSS = mpcObj.TruncVec(finalBurdenSS, mpcObj.GetDataBits(), mpcObj.GetFracBits())
	if scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS); ok {
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
	null, nullRSS, X, y0 := ast.nullSetup()
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
		q, bl, _ := ast.blockStat(b, nsnps, null, X, y0, nil)
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
	if scaleSS, ok := ast.general.rareVariantScaleShares(nullRSS); ok {
		scaleVec := mpc_core.InitRVec(rtype.Zero(), nB)
		for b := 0; b < nB; b++ {
			scaleVec[b] = scaleSS[0].Copy()
		}
		qBlockSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(qBlockSS, scaleVec), mpcObj.GetDataBits(), mpcObj.GetFracBits())
		burdenSqSS = mpcObj.TruncVec(mpcObj.SSMultElemVec(burdenSqSS, scaleVec), mpcObj.GetDataBits(), mpcObj.GetFracBits())
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
	total := 0
	for p := 1; p <= ast.general.config.NumMainParties; p++ {
		total += nrows[p]
	}
	return total
}
