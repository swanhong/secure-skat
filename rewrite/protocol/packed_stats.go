package protocol

import (
	mpc_core "github.com/hhcho/mpc-core"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

type geneBatchTerms struct {
	localGtx     []*mat.Dense
	gtxEncoded   securecrypto.PlainVector
	localGtG     []*mat.Dense
	packedWeight securecrypto.CipherVector
	maskEncoded  securecrypto.PlainVector
	gvLocal      []gvLocalGene
	gvGene       []gvGeneShares
}

func ComputePackedStatistics(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	gp []*mat.Dense,
	gv []*mat.Dense,
	x *mat.Dense,
	y0 *mat.Dense,
	beta mpc_core.RMat,
	xtxInv mpc_core.RMat,
	weight []mpc_core.RVec,
	signedWeight []mpc_core.RVec,
	seed int64,
	timingObservers ...TimingObserver,
) (
	gpQ, gpL, gvQ, gvL mpc_core.RMat,
	geneV, geneS1, geneS2, geneS3 mpc_core.RVec,
) {
	timingObserver := selectTimingObserver(timingObservers)
	/*
		Compute all packed phenotype and kernel statistics.

		For every gene g and phenotype t:
		    gpQ[g,t], gpL[g,t] = public-variant Q and L
		    gvQ[g,t], gvL[g,t] = private-variant Q and L

		For every gene g:
		    V[g]  = Burden variance
		    S1[g] = Trace(K)
		    S2[g] = Trace(K²)
		    S3[g] = Trace(K³)

		All outputs remain secret-shared and use global gene order.
	*/
	geneCount := len(dataParams.Genes)
	phenotypeCount := dataParams.PhenotypeCount
	rtype := mpcObj.GetRType()

	gpQ = mpc_core.InitRMat(rtype.Zero(), geneCount, phenotypeCount)
	gpL = mpc_core.InitRMat(rtype.Zero(), geneCount, phenotypeCount)
	gvQ = mpc_core.InitRMat(rtype.Zero(), geneCount, phenotypeCount)
	gvL = mpc_core.InitRMat(rtype.Zero(), geneCount, phenotypeCount)

	geneV = mpc_core.InitRVec(rtype.Zero(), geneCount)
	geneS1 = mpc_core.InitRVec(rtype.Zero(), geneCount)
	geneS2 = mpc_core.InitRVec(rtype.Zero(), geneCount)
	geneS3 = mpc_core.InitRVec(rtype.Zero(), geneCount)

	// 1. Pack each shared beta column once for reuse by every gene batch.
	betaByPhenotype := make([]mpc_core.RVec, phenotypeCount)
	packedBeta := make([]securecrypto.CipherVector, phenotypeCount)

	for phenotype := 0; phenotype < phenotypeCount; phenotype++ {
		finishTiming := startTiming(
			timingObserver, "pack_beta", "packed_statistics", phenotype,
		)
		betaByPhenotype[phenotype] = mpc_core.InitRVec(rtype.Zero(), dataParams.C)
		for covariate := 0; covariate < dataParams.C; covariate++ {
			betaByPhenotype[phenotype][covariate] = beta[covariate][phenotype].Copy()
		}
		packedBeta[phenotype] = PackBeta(mpcObj, heParams, betaByPhenotype[phenotype])
		finishTiming()
	}

	for _, batch := range cryptoParams.Batches {
		// 2. Prepare the phenotype-independent values for this gene batch.
		finishPrepare := startTiming(
			timingObserver,
			"prepare_gene_batch",
			"packed_statistics",
			NoTimingPhenotype,
		)
		terms := PrepareGeneBatch(
			mpcObj, heParams, dataParams, cryptoParams, batch,
			gp, gv, x, signedWeight, timingObserver,
		)
		finishPrepare()

		// 3. Compute and store the phenotype-dependent public/private Q and L.
		for phenotype := 0; phenotype < phenotypeCount; phenotype++ {
			finishPhenotype := startTiming(
				timingObserver,
				"phenotype_pass",
				"packed_statistics",
				phenotype,
			)
			var yColumn mat.Vector
			if mpcObj.GetPid() != auxiliaryPartyID {
				yColumn = y0.ColView(phenotype)
			}

			finishPublic := startTiming(
				timingObserver, "compute_gp_ql", "phenotype_pass", phenotype,
			)
			batchGpQ, batchGpL := ComputeGpQL(
				mpcObj, heParams, dataParams, cryptoParams, batch,
				gp, yColumn, terms.gtxEncoded, packedBeta[phenotype],
				terms.packedWeight, terms.maskEncoded,
			)
			finishPublic()

			finishPrivate := startTiming(
				timingObserver, "compute_gv_ql", "phenotype_pass", phenotype,
			)
			batchGvQ, batchGvL := ComputeGvQL(
				mpcObj, dataParams, batch, gv, yColumn,
				terms.gvLocal, terms.gvGene, betaByPhenotype[phenotype],
			)
			finishPrivate()

			for position, geneIndex := range batch.GeneIndices {
				gpQ[geneIndex][phenotype] = batchGpQ[position].Copy()
				gpL[geneIndex][phenotype] = batchGpL[position].Copy()
				gvQ[geneIndex][phenotype] = batchGvQ[position].Copy()
				gvL[geneIndex][phenotype] = batchGvL[position].Copy()
			}
			finishPhenotype()
		}

		// 4. Compute and store the phenotype-independent Burden variance and kernel traces.
		finishKernel := startTiming(
			timingObserver,
			"kernel_statistics",
			"packed_statistics",
			NoTimingPhenotype,
		)
		batchV, batchS1, batchS2, batchS3 := ComputeGeneBatchKernelStatistics(
			mpcObj, heParams, dataParams, cryptoParams, batch,
			gp, gv, x, terms.localGtx, terms.localGtG,
			terms.gvLocal, terms.gvGene, xtxInv, weight, signedWeight, seed,
			timingObserver,
		)
		finishKernel()

		for position, geneIndex := range batch.GeneIndices {
			geneV[geneIndex] = batchV[position].Copy()
			geneS1[geneIndex] = batchS1[position].Copy()
			geneS2[geneIndex] = batchS2[position].Copy()
			geneS3[geneIndex] = batchS3[position].Copy()
		}
	}

	return gpQ, gpL, gvQ, gvL, geneV, geneS1, geneS2, geneS3
}

func PrepareGeneBatch(
	mpcObj *mpc.MPC,
	heParams *securecrypto.CryptoParams,
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	gp []*mat.Dense,
	gv []*mat.Dense,
	x *mat.Dense,
	signedWeight []mpc_core.RVec,
	timingObservers ...TimingObserver,
) geneBatchTerms {
	timingObserver := selectTimingObserver(timingObservers)
	/*
		Prepare the phenotype-independent values reused by one gene batch.

		Cohort-local:
		    localGtx = Transpose(Gp) * X
		    localGtG = Transpose(Gp) * Gp

		Encrypted:
		    packedWeight

		B-local:
		    gvLocal

		Shared:
		    gvGene
	*/
	var terms geneBatchTerms

	// 1. Encode the public active-slot mask.
	finishMask := startTiming(
		timingObserver,
		"active_mask_encode",
		"prepare_gene_batch",
		NoTimingPhenotype,
	)
	terms.maskEncoded, _ = securecrypto.EncodeFloatVector(
		heParams, ActiveMask(dataParams, cryptoParams, batch),
	)
	finishMask()

	// 2. Compute cohort-local public-variant terms.
	if mpcObj.GetPid() != auxiliaryPartyID {
		finishGtx := startTiming(
			timingObserver,
			"compute_gtx",
			"prepare_gene_batch",
			NoTimingPhenotype,
		)
		terms.localGtx, terms.gtxEncoded = ComputeGtX(
			heParams, dataParams, cryptoParams, batch, gp, x,
		)
		finishGtx()

		finishLocalGtg := startTiming(
			timingObserver,
			"compute_local_gtg",
			"prepare_gene_batch",
			NoTimingPhenotype,
		)
		terms.localGtG = ComputeLocalGtG(dataParams, batch, gp)
		finishLocalGtg()
	}

	// 3. Pack the shared public-variant weights.
	finishWeight := startTiming(
		timingObserver,
		"pack_weight",
		"prepare_gene_batch",
		NoTimingPhenotype,
	)
	terms.packedWeight = PackWeight(
		mpcObj, heParams, dataParams, cryptoParams, batch, signedWeight,
	)
	finishWeight()

	// 4. Prepare the private-variant terms.
	finishPrivate := startTiming(
		timingObserver,
		"prepare_private_gene_terms",
		"prepare_gene_batch",
		NoTimingPhenotype,
	)
	terms.gvLocal, terms.gvGene = PrepareGvGeneTerms(
		mpcObj, dataParams, batch, gv, gp, x,
	)
	finishPrivate()

	return terms
}
