package gwas

import (
	"bufio"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"go.dedis.ch/onet/v3/log"

	"fmt"

	mpc_core "github.com/hhcho/mpc-core"

	"github.com/tuneinsight/lattigo/v6/examples"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"gonum.org/v1/gonum/mat"
)

type ProtocolInfo struct {
	mpcObj mpc.ParallelMPC
	cps    *crypto.CryptoParams
	cpsPar []*crypto.CryptoParams // One per thread

	// Application-layer traffic used by collective HE key setup. The remaining
	// pre-run communication is derived when skat_fed starts.
	fedSetupComm mpc.CommunicationStats

	// Input files
	genoBlocks     []*GenoFileStream
	genoBlockSizes []int
	pheno          *mat.Dense
	cov            *mat.Dense
	pos            []uint64

	gwasParams *GWASParams

	config *Config
}

type Config struct {
	NumMainParties int `toml:"num_main_parties"`
	HubPartyId     int `toml:"hub_party_id"`

	CkksParams string `toml:"ckks_params"`

	DivSqrtMaxLen int `toml:"div_sqrt_max_len"`

	NumInds    []int `toml:"num_inds"`
	NumSnps    int   `toml:"num_snps"`
	NumCovs    int   `toml:"num_covs"`
	CovAllOnes bool  `toml:"cov_all_ones"`

	ItersPerEval  int `toml:"iter_per_eigenval"`
	NumPCs        int `toml:"num_pcs_to_remove"`
	NumOversample int `toml:"num_oversampling"`
	NumPowerIters int `toml:"num_power_iters"`

	SkipQC             bool `toml:"skip_qc"`
	SkipPCA            bool `toml:"skip_pca"`
	UseCachedQC        bool `toml:"use_cached_qc"`
	UseCachedPCA       bool `toml:"use_cached_pca"`
	UseCachedCombinedQ bool `toml:"use_cached_combined_q"`
	SkipPowerIter      bool `toml:"skip_power_iter"`
	PCARestartIter     int  `toml:"restart_pca_from_iter"`

	IndMissUB    float64 `toml:"imiss_ub"`
	HetLB        float64 `toml:"het_lb"`
	HetUB        float64 `toml:"het_ub"`
	SnpMissUB    float64 `toml:"gmiss"`
	MafLB        float64 `toml:"maf_lb"`
	HweUB        float64 `toml:"hwe_ub"`
	SnpDistThres int     `toml:"snp_dist_thres"`

	BindingIP string `toml:"binding_ipaddr"`
	Servers   map[string]mpc.Server

	SharedKeysPath string `toml:"shared_keys_path"`

	GenoFileFormat string `toml:"geno_file_format"`        // 'blocks' or 'pgen'
	GenoFilePrefix string `toml:"geno_binary_file_prefix"` // If 'pgen' expects a '%d' placeholder for chrom, e.g. "ukb_imp_chr%d_v3" (.pgen/.psam/.pvar)

	GenoNumBlocks     int    `toml:"geno_num_blocks"`
	GenoBlockSizeFile string `toml:"geno_block_size_file"`

	PhenoFile  string `toml:"pheno_file"`
	CovFile    string `toml:"covar_file"`
	SnpPosFile string `toml:"snp_position_file"`

	// skat_fed: private party id + fmt prefix of its per-gene private-variant block files.
	PrivatePid        int    `toml:"private_pid"`
	PrivateOnlyPrefix string `toml:"private_only_prefix"`

	// skat_fed SKAT p-value (Wilson-Hilferty, screening): public trace-column budget. When the public
	// variant count is no larger, the standard basis gives exact traces; otherwise use Hutchinson.
	// Zero disables the SKAT p-value path.
	// The private-variant moment contribution is block-contracted (no padding, count hidden).
	SkatPValueProbes int `toml:"skat_pvalue_probes"`

	// RotKeyPow2Only generates only power-of-two rotation keys (InnerSumAll), skipping the
	// baby-step-giant-step/matmul keys — ~8× less rot-key RAM. Safe for skat_fed (no matmul/PCA).
	RotKeyPow2Only bool `toml:"rotkey_pow2only"`

	UsePrecomputedGenoCount bool   `toml:"use_precomputed_geno_count"`
	GenoCountFile           string `toml:"geno_count_file"`
	SampleKeepFile          string `toml:"sample_keep_file"`
	SnpIdsFile              string `toml:"snp_ids_file"`

	OutDir   string `toml:"output_dir"`
	CacheDir string `toml:"cache_dir"`

	MpcFieldSize                int    `toml:"mpc_field_size"`
	MpcDataBits                 int    `toml:"mpc_data_bits"`
	MpcFracBits                 int    `toml:"mpc_frac_bits"`
	MpcNumThreads               int    `toml:"mpc_num_threads"`
	MpcBooleanShares            bool   `toml:"mpc_boolean_shares"`
	LocalNumThreads             int    `toml:"local_num_threads"`
	LocalAssocNumBlocksParallel int    `toml:"assoc_num_blocks_parallel"`
	MemoryLimit                 uint64 `toml:"memory_limit"`

	// Logistic regression specific
	UseLogistic     bool    `toml:"use_logistic"`
	InverseMatScale float64 `toml:"inverse_mat_scale"`
	A               float64 `toml:"A"`
	B               float64 `toml:"B"`
	Degree          int     `toml:"degree"`
	Epochs          int     `toml:"epochs"`

	// BinaryPheno marks a binary {0,1} (case-control) phenotype; the SKAT null model
	// then centers y by 0.5. A public study-design property, not derived from data.
	BinaryPheno bool `toml:"binary_pheno"`

	Debug          bool  `toml:"debug"`
	BlocksForAssoc []int `toml:"blocks_for_assoc_test"`
	PgenBatchSize  int   `toml:"pgen_batch_nsnp"`
}

func (prot *ProtocolInfo) IsBlockForAssocTest(blockId int) bool {
	if len(prot.config.BlocksForAssoc) == 0 { // if empty test all blocks
		return true
	}
	for _, v := range prot.config.BlocksForAssoc {
		if v == blockId {
			return true
		}
	}
	return false
}

func (prot *ProtocolInfo) IsPgen() bool {
	if prot.config.GenoFileFormat == "pgen" {
		return true
	} else if prot.config.GenoFileFormat == "blocks" {
		return false
	} else {
		panic(fmt.Sprint("Unsupported geno_file_format:", prot.config.GenoFileFormat))
	}
}

func (prot *ProtocolInfo) GetCryptoParams() *crypto.CryptoParams {
	return prot.cps
}

func (prot *ProtocolInfo) GetMpc() mpc.ParallelMPC {
	return prot.mpcObj
}

func (prot *ProtocolInfo) GetConfig() *Config {
	return prot.config
}

func (prot *ProtocolInfo) GetGwasParams() *GWASParams {
	return prot.gwasParams
}

func (prot *ProtocolInfo) GetGenoBlocks() []*GenoFileStream {
	return prot.genoBlocks
}

func InitializeGWASProtocol(config *Config, pid int, mpcOnly bool) (gwasProt *ProtocolInfo) {
	tInit := time.Now()
	defer func() { fedTimings.initTotal = time.Since(tInit) }() // network + keygen + geno/cov load
	var params ckks.Parameters
	if !mpcOnly {
		var paramsLiteral ckks.ParametersLiteral
		switch config.CkksParams {
		case "PN12QP109":
			paramsLiteral = examples.CKKSComplexParamsN12QP109
		case "PN13QP218":
			paramsLiteral = examples.CKKSComplexParamsN13QP218
		case "PN14QP438":
			paramsLiteral = examples.CKKSComplexParamsN14QP438
		case "PN15QP880":
			paramsLiteral = examples.CKKSComplexParamsN15QP881
		case "PN16QP1761":
			paramsLiteral = examples.CKKSComplexParamsPN16QP1761
		default:
			panic("Undefined value of CKKS params in config")
		}
		var err error
		params, err = ckks.NewParametersFromLiteral(paramsLiteral)
		if err != nil {
			panic(err)
		}
	}

	prec := uint(config.MpcFieldSize)
	networks := mpc.ParallelNetworks(mpc.InitCommunication(config.BindingIP, config.Servers, pid, config.NumMainParties+1, config.MpcNumThreads, config.SharedKeysPath))
	// Start the reported protocol traffic after TCP channels/PRG initialization,
	// but before collective HE key setup. This yields setup+compute+decrypt bytes.
	networks.ResetNetworkLog()

	var rtype mpc_core.RElem
	switch config.MpcFieldSize {
	case 256:
		rtype = mpc_core.LElem256Zero
	case 128:
		rtype = mpc_core.LElem128Zero
	default:
		panic("Unsupported value of MPC field size")
	}

	log.LLvl1(fmt.Sprintf("MPC parameters: bit length %d, data bits %d, frac bits %d",
		config.MpcFieldSize, config.MpcDataBits, config.MpcFracBits))
	mpcEnv := mpc.InitParallelMPCEnv(networks, rtype, config.MpcDataBits, config.MpcFracBits)
	for thread := range mpcEnv {
		mpcEnv[thread].SetHubPid(config.HubPartyId)
		mpcEnv[thread].SetBooleanShareFlag(config.MpcBooleanShares)
		mpcEnv[thread].SetDivSqrtMaxLen(config.DivSqrtMaxLen)
	}

	var cps *crypto.CryptoParams
	setupCommStart := networks.GetCommunicationStats()
	if !mpcOnly {
		cps = networks.CollectiveInit(&params, prec, config.RotKeyPow2Only)
	}
	fedSetupComm := networks.GetCommunicationStats().Sub(setupCommStart)

	var pheno, cov *mat.Dense
	var pos []uint64
	var genofs []*GenoFileStream
	var genoBlockSizes []int

	isPgen := config.GenoFileFormat == "pgen"

	genofs = make([]*GenoFileStream, config.GenoNumBlocks)
	genoBlockSizes = make([]int, config.GenoNumBlocks)

	if pid > 0 {
		// Read geno block size file
		file, err := os.Open(config.GenoBlockSizeFile)

		if err != nil {
			log.Fatalf("failed to open: %v", config.GenoBlockSizeFile)
		}
		scanner := bufio.NewScanner(file)
		scanner.Split(bufio.ScanLines)

		for i := 0; i < config.GenoNumBlocks; i++ {
			if !scanner.Scan() {
				log.Fatalf("not enough lines in %v", config.GenoBlockSizeFile)
			}

			genoBlockSizes[i], err = strconv.Atoi(scanner.Text())
			if err != nil {
				log.Fatalf("parse error: %v", config.GenoBlockSizeFile)
			}
		}

		if scanner.Scan() {
			log.Fatalf("too many lines in %v", config.GenoBlockSizeFile)
		}

		file.Close()

		totalSize := 0
		for _, v := range genoBlockSizes {
			totalSize += v
		}
		if totalSize != config.NumSnps {
			log.Fatalf("Sum of block sizes does not match number of snps")
		}

		if !isPgen {
			// Create file streams for geno block files
			for i := range genofs {
				filename := fmt.Sprintf("%s.%d.bin", config.GenoFilePrefix, i)
				log.LLvl1(time.Now().Format(time.RFC3339), "Opening geno file:", filename)
				genofs[i] = NewGenoFileStream(filename, uint64(config.NumInds[pid]), uint64(genoBlockSizes[i]), false)
			}
		}

		tab := '\t'
		pheno = LoadMatrixFromFile(config.PhenoFile, tab)
		cov = LoadMatrixFromFile(config.CovFile, tab)
		pos = LoadSNPPositionFile(config.SnpPosFile, tab)
		log.LLvl1(time.Now().Format(time.RFC3339), "First few SNP positions:", pos[:5])
	}

	gwasParams := InitGWASParams(config.NumInds, config.NumSnps, config.NumCovs, config.NumPCs, config.SnpDistThres)

	return &ProtocolInfo{
		mpcObj:       mpcEnv, // One MPC object for each thread
		cps:          cps,
		fedSetupComm: fedSetupComm,

		genoBlocks:     genofs,
		genoBlockSizes: genoBlockSizes,
		pheno:          pheno,
		cov:            cov,
		pos:            pos,

		gwasParams: gwasParams,
		config:     config,
	}
}

func (g *ProtocolInfo) Phase1() {
	net := g.mpcObj.GetNetworks()

	net.ResetNetworkLog()

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting QC")

	filterParams := InitFilteringSettings(g.config.MafLB, g.config.HweUB, g.config.SnpMissUB, g.config.IndMissUB, g.config.HetLB, g.config.HetUB)
	qc := g.InitQC(filterParams)
	if g.config.SkipQC && !g.config.UseCachedQC {

		qc.filtNumSnps = g.gwasParams.NumSNP()
		qc.filtNumInds = g.gwasParams.NumInds()
		g.gwasParams.SetFiltCounts(qc.filtNumInds, qc.filtNumSnps)
		log.LLvl1(time.Now().Format(time.RFC3339), "Individual and SNP filters skipped")

	} else {

		if g.config.UsePrecomputedGenoCount { // Use precomputed geno count file
			qc.QualityControlProtocolWithPrecomputedGenoStats(g.config.UseCachedQC)
		} else {
			qc.QualityControlProtocol(g.config.UseCachedQC)
		}

	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Finished QC")

	net.PrintNetworkLog()
}

func (g *ProtocolInfo) Phase2() crypto.CipherMatrix {
	net := g.mpcObj.GetNetworks()
	pid := g.mpcObj[0].GetPid()

	net.ResetNetworkLog()

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting PCA")

	var Qpca crypto.CipherMatrix
	pcaCacheFile := g.CachePath("Qpc.txt")
	if g.config.UseCachedPCA {

		if pid > 0 {
			// TODO cache ciphertexts instead
			mat := LoadMatrixFromFileFloat(pcaCacheFile, ',')
			Qpca, _, _, _ = crypto.EncryptFloatMatrixRow(g.cps, mat)
		} else {
			Qpca = make(crypto.CipherMatrix, g.config.NumPCs)
		}

	} else if g.config.SkipPCA { // Qpca = nil

		g.config.NumPCs = 0
		g.gwasParams.SetNumPC(0)
		log.LLvl1(time.Now().Format(time.RFC3339), "PCA skipped")

	} else {

		g.gwasParams.SetPopStratMethod(true)
		Qpca = g.PopulationStratification()
		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("PCA complete: calculated %d PCs", len(Qpca)))

		// TODO cache ciphertexts instead of decrypting
		for p := 1; p <= g.config.NumMainParties; p++ {
			SaveMatrixToFile(g.cps, g.mpcObj[0], Qpca, g.gwasParams.numFiltInds[p], p, pcaCacheFile)
		}

	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Finished PCA")

	net.PrintNetworkLog()

	return Qpca
}

func (g *ProtocolInfo) Phase3(Qpca crypto.CipherMatrix) {
	net := g.mpcObj.GetNetworks()

	net.ResetNetworkLog()

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting association tests")

	assoc, outFilter := g.ComputeAssocStatistics(Qpca)

	log.LLvl1(time.Now().Format(time.RFC3339), "Finished association tests")

	net.PrintNetworkLog()

	// Collective decrypt and save to file
	if g.mpcObj[0].GetPid() > 0 {
		assocDec := g.mpcObj[0].Network.CollectiveDecryptVec(g.cps, assoc, -1)
		out := crypto.DecodeFloatVector(g.cps, assocDec)

		outFinal := make([]float64, SumBool(outFilter))
		index := 0
		for i := range outFilter {
			if outFilter[i] {
				outFinal[index] = out[i]
				index++
			}
		}

		SaveFloatVectorToFile(g.OutPath("assoc.txt"), outFinal)
	}
	log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Output collectively decrypted and saved to: %s", g.OutPath("assoc.txt")))
}

func (g *ProtocolInfo) GWAS() {

	log.LLvl1(time.Now().Format(time.RFC3339), "Starting GWAS protocol")

	g.Phase1()
	Qpca := g.Phase2()
	g.Phase3(Qpca)
}

// SKAT is the main entry point to run the Secure SKAT protocol independently of single-variant GWAS
// runRareVariant runs QC + the per-block low-rank SKAT, returning the decrypted per-block
// SKAT (Q) and Burden statistics (slot b = block b). nil on pid 0.
func (g *ProtocolInfo) runRareVariant() (skat, burden []float64) {
	mpcObj := g.mpcObj[0]
	pid := mpcObj.GetPid()

	log.LLvl1(time.Now().Format(time.RFC3339), "Running rare-variant test (low-rank per-block)")

	// skip-QC sets full counts directly, otherwise run Phase1 (QC) to populate them.
	if len(g.gwasParams.FiltNumInds()) != g.config.NumMainParties+1 {
		if g.config.SkipQC && !g.config.UseCachedQC {
			g.gwasParams.SetFiltCounts(g.gwasParams.NumInds(), g.gwasParams.NumSNP())
			log.LLvl1(time.Now().Format(time.RFC3339), "Individual and SNP filters skipped")
		} else {
			g.Phase1()
		}
	}

	assoc, burdenEnc := g.ComputeSKATStatistics()
	g.mpcObj.GetNetworks().PrintNetworkLog()

	if pid > 0 {
		nB := g.config.GenoNumBlocks
		skat = crypto.DecodeFloatVector(g.cps, mpcObj.Network.CollectiveDecryptVec(g.cps, assoc, -1))[:nB]
		burden = crypto.DecodeFloatVector(g.cps, mpcObj.Network.CollectiveDecryptVec(g.cps, burdenEnc, -1))[:nB]
	}
	return
}

// runFederatedPrivate runs skat_fed (private party loads its blocks); per-gene SKAT statistic and the
// Burden p-value (p = erfc(√(T/2)) ~ χ²₁, R::SKAT r.corr=1). The Burden statistic and zᵀPz stay
// secret-shared (only √(T/2) is revealed), so neither is disclosed. nil on pid 0.
func (g *ProtocolInfo) runFederatedPrivate() (skatOut, burdenPOut, skatPOut []float64) {
	mpcObj := g.mpcObj[0]
	pid := mpcObj.GetPid()
	privatePid := g.config.PrivatePid
	networks := g.mpcObj.GetNetworks()
	commStart := networks.GetCommunicationStats()
	defer networks.PrintCommunicationDelta("skat_fed_run", commStart)
	mode := "raw_q"
	if g.config.SkatPValueProbes > 0 {
		mode = "skat_p"
	}
	metrics := newFedRunMetrics(networks, mode, g.fedSetupComm)

	var privateOnly []*mat.Dense
	loadMark := metrics.mark()
	if pid == privatePid && g.config.PrivateOnlyPrefix != "" {
		n, _ := g.pheno.Dims()
		blocks, err := loadDenseBlocks(g.config.PrivateOnlyPrefix, g.config.GenoNumBlocks, n)
		if err != nil {
			panic(err)
		}
		privateOnly = blocks
		metrics.end("load_private_only", loadMark)
	} else if pid == privatePid {
		log.LLvl1("[skat_fed] WARNING: private_pid set but private_only_prefix empty -> PART B (private variants) is EMPTY")
	}

	assocMark := metrics.mark()
	assocTest := g.initSKAT()
	assocTest.fedMetrics = metrics
	metrics.end("assoc_init", assocMark)
	skat, sqrtT2Enc, skatZEnc := assocTest.ComputeSKATFederatedPrivate(privateOnly, privatePid)
	decryptMark := metrics.mark()
	if pid == 0 {
		secureRun := metrics.finishRun()
		printFedTimingRecords(metrics, secureRun)
		metrics.printRecords()
		return nil, nil, nil
	}
	nB := g.config.GenoNumBlocks
	// Q_skat is nil in SKAT-p mode (only the equivalent z is released) → skip its decrypt.
	if len(skat) > 0 {
		skatOut = crypto.DecodeFloatVector(g.cps, mpcObj.Network.CollectiveDecryptVec(g.cps, skat, -1))[:nB]
	}
	sqrtT2 := crypto.DecodeFloatVector(g.cps, mpcObj.Network.CollectiveDecryptVec(g.cps, sqrtT2Enc, -1))[:nB]
	burdenPOut = make([]float64, nB)
	for b := 0; b < nB; b++ {
		burdenPOut[b] = 1.0 // p=1 at √(T/2)≤0 (CKKS noise on a ~0 statistic); erfc(0)=1
		if sqrtT2[b] > 0 {
			burdenPOut[b] = math.Erfc(sqrtT2[b]) // ∈(0,1)
		}
	}
	// SKAT p-value (screening): only z is revealed; p = ½erfc(z/√2). nil if SKAT p disabled (no probes).
	if len(skatZEnc) > 0 {
		z := crypto.DecodeFloatVector(g.cps, mpcObj.Network.CollectiveDecryptVec(g.cps, skatZEnc, -1))[:nB]
		skatPOut = make([]float64, nB)
		for b := 0; b < nB; b++ {
			skatPOut[b] = 0.5 * math.Erfc(z[b]/math.Sqrt2)
		}
	}
	metrics.end("decrypt_outputs", decryptMark)
	secureRun := metrics.finishRun()
	log.LLvl1(fmt.Sprintf("[skat_fed] decrypt: %v | total run: %v",
		metrics.stageDuration("decrypt_outputs").Round(time.Millisecond), secureRun.Round(time.Millisecond)))
	logFedTimingTree(metrics, secureRun)
	printFedTimingRecords(metrics, secureRun)
	metrics.printRecords()
	return skatOut, burdenPOut, skatPOut
}

func logFedTimingTree(metrics *fedRunMetrics, secureRun time.Duration) {
	ms := func(d time.Duration) time.Duration { return d.Round(time.Millisecond) }
	st := mpc.SetupTiming
	setup := st.PubKey + st.RelinKey + st.RotKey
	netload := nonNegativeDuration(fedTimings.initTotal - setup)
	nullAgg := metrics.sumDuration("null_aggregate_xtx", "null_aggregate_xty", "null_aggregate_yty")
	public := metrics.sumDuration("gene_public_score_gtx_gty", "gene_public_weight_dosage_maf", "gene_public_stat_share")
	privateLocal := metrics.stageDuration("gene_local_private_gtg_gtx_gty")
	privateShare := metrics.stageDuration("gene_private_stat_share")
	burden := metrics.sumDuration("gene_burden_public_gtg", "gene_burden_public_gtx",
		"gene_burden_private_cross", "gene_burden_projection")
	moments := metrics.sumDuration("gene_moments_setup_gtx", "gene_moments_public_trace_gtg_gtx",
		"gene_moments_private_cross")
	finalize := metrics.parentLeafDuration("post_block_finalize", "")
	grand := fedTimings.initTotal + secureRun
	for _, ln := range []string{
		"[skat_fed] timing tree:",
		fmt.Sprintf("  ├─ init                    %v", ms(fedTimings.initTotal)),
		fmt.Sprintf("  │    ├─ crypto setup       %v", ms(setup)),
		fmt.Sprintf("  │    │    ├─ PubKeyGen      %v", ms(st.PubKey)),
		fmt.Sprintf("  │    │    ├─ RelinKeyGen    %v", ms(st.RelinKey)),
		fmt.Sprintf("  │    │    └─ RotKeyGen      %v", ms(st.RotKey)),
		fmt.Sprintf("  │    └─ network+geno/cov   %v", ms(netload)),
		fmt.Sprintf("  ├─ load privateOnly        %v", ms(metrics.stageDuration("load_private_only"))),
		fmt.Sprintf("  ├─ assoc init              %v", ms(metrics.stageDuration("assoc_init"))),
		fmt.Sprintf("  ├─ compute                 %v", ms(fedTimings.total)),
		fmt.Sprintf("  │    ├─ null model         %v", ms(fedTimings.nullTotal)),
		fmt.Sprintf("  │    │    ├─ aggregate      %v", ms(nullAgg)),
		fmt.Sprintf("  │    │    ├─ inverse        %v", ms(metrics.stageDuration("null_solve"))),
		fmt.Sprintf("  │    │    ├─ beta           %v", ms(metrics.stageDuration("null_beta_pack"))),
		fmt.Sprintf("  │    │    └─ RSS            %v", ms(metrics.stageDuration("null_rss"))),
		fmt.Sprintf("  │    ├─ pre-block setup    %v", ms(metrics.stageDuration("pre_block_setup"))),
		fmt.Sprintf("  │    ├─ blocks (%d)         %v  [per-block %s]", len(fedTimings.blockSecs), ms(fedTimings.blocks), blockSecStats(fedTimings.blockSecs)),
		fmt.Sprintf("  │    │    ├─ shape sync      %v", ms(metrics.stageDuration("gene_shape_sync"))),
		fmt.Sprintf("  │    │    ├─ local contract  %v", ms(metrics.stageDuration("gene_local_public_gtg_gtx_gty"))),
		fmt.Sprintf("  │    │    ├─ public stat     %v", ms(public)),
		fmt.Sprintf("  │    │    ├─ private total   %v", ms(privateLocal+privateShare)),
		fmt.Sprintf("  │    │    │    ├─ local       %v", ms(privateLocal)),
		fmt.Sprintf("  │    │    │    └─ stat/share  %v", ms(privateShare)),
		fmt.Sprintf("  │    │    ├─ burden variance %v", ms(burden)),
		fmt.Sprintf("  │    │    ├─ moments         %v", ms(moments)),
		fmt.Sprintf("  │    │    └─ other overhead  %v", ms(metrics.stageDuration("gene_other"))),
		fmt.Sprintf("  │    └─ post-block finalize %v", ms(finalize)),
		fmt.Sprintf("  ├─ decrypt                 %v", ms(metrics.stageDuration("decrypt_outputs"))),
		fmt.Sprintf("  ├─ run overhead            %v", ms(metrics.stageDuration("run_other"))),
		fmt.Sprintf("  └─ TOTAL                   %v", ms(grand)),
	} {
		log.LLvl1(ln)
	}
}

// SKATFederatedPrivate saves the per-gene SKAT statistic and Burden p-value (runFederatedPrivate
// details the reveal set; Burden statistic and zᵀPz are never released).
func (g *ProtocolInfo) SKATFederatedPrivate() {
	skat, burdenP, skatP := g.runFederatedPrivate()
	// Cumulative application-layer traffic from collective key setup through
	// secure computation and output decryption. The final shutdown barrier is excluded.
	mode := "raw_q"
	if g.config.SkatPValueProbes > 0 {
		mode = "skat_p"
	}
	g.mpcObj.GetNetworks().PrintCommunicationSummaryWithMode("skat_fed_total", mode)
	if g.mpcObj[0].GetPid() > 0 {
		if len(skat) > 0 { // nil/empty in SKAT-p mode (Q_skat withheld; only the z-derived p is released)
			SaveFloatVectorToFile(g.OutPath("skat_fed_out.txt"), skat)
		}
		SaveFloatVectorToFile(g.OutPath("skat_fed_burden_p_out.txt"), burdenP)
		if len(skatP) > 0 {
			SaveFloatVectorToFile(g.OutPath("skat_fed_skat_p_out.txt"), skatP)
		}
	}
}

// SKAT runs the per-block SKAT test, saving per-block Q and Burden.
func (g *ProtocolInfo) SKAT() {
	skat, burden := g.runRareVariant()
	if g.mpcObj[0].GetPid() > 0 {
		SaveFloatVectorToFile(g.OutPath("skat_out.txt"), skat)
		SaveFloatVectorToFile(g.OutPath("burden_out.txt"), burden)
	}
}

// Burden runs the per-block Burden test, saving per-block Burden.
func (g *ProtocolInfo) Burden() {
	_, burden := g.runRareVariant()
	if g.mpcObj[0].GetPid() > 0 {
		SaveFloatVectorToFile(g.OutPath("burden_out.txt"), burden)
	}
}

// SKATO runs the per-block SKAT-O: skato_b = (1-rho)*skat_b + rho*burden_b, blended from
// the public per-block statistics.
func (g *ProtocolInfo) SKATO(rho float64) {
	skat, burden := g.runRareVariant()
	if g.mpcObj[0].GetPid() > 0 {
		skato := make([]float64, len(skat))
		for b := range skato {
			skato[b] = (1-rho)*skat[b] + rho*burden[b]
		}
		SaveFloatVectorToFile(g.OutPath("skat_out.txt"), skat)
		SaveFloatVectorToFile(g.OutPath("burden_out.txt"), burden)
		SaveFloatVectorToFile(g.OutPath("skato_out.txt"), skato)
	}
}

// SetPhenoAndCov allows deterministic testing scripts to override File I/O datasets in memory
func (g *ProtocolInfo) SetPhenoAndCov(pheno, cov *mat.Dense) {
	g.pheno = pheno
	g.cov = cov
}

// ComputeSKATStatistics runs the low-rank per-block secure SKAT: returns per-block Q and
// Burden (slot b = block b's per-gene statistic).
func (g *ProtocolInfo) ComputeSKATStatistics() (crypto.CipherVector, crypto.CipherVector) {
	assocTest := g.initSKAT()
	return assocTest.ComputeSKATStatisticsPerBlock()
}

// rareVariantScaleShares returns the shared 1/(2σ̂²) = (dof/2)/RSS scale factor (dof = N−c),
// and false if RSS/dof is unavailable. RSS stays secret-shared from the null solve through division.
func (g *ProtocolInfo) rareVariantScaleShares(rssSS mpc_core.RElem) (mpc_core.RVec, bool) {
	mpcObj := g.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()
	scaleSS := mpc_core.InitRVec(rtype.Zero(), 1)
	if rssSS == nil {
		return scaleSS, false
	}

	nrowsAll := g.gwasParams.FiltNumInds()
	if len(nrowsAll) != g.config.NumMainParties+1 {
		nrowsAll = g.config.NumInds
	}
	totalInds := 0
	for p := 1; p <= g.config.NumMainParties; p++ {
		totalInds += nrowsAll[p]
	}
	dof := totalInds - (g.gwasParams.NumCov() + 1)
	if dof <= 0 {
		return scaleSS, false
	}

	numerSS := mpc_core.InitRVec(rtype.Zero(), 1)
	if pid == mpcObj.GetHubPid() {
		numerSS[0] = rtype.FromFloat64(float64(dof)/2.0, mpcObj.GetFracBits())
	}
	return mpcObj.Divide(numerSS, mpc_core.RVec{rssSS}, false), true
}

// scaleRareVariantShareStat multiplies a shared statistic by the shared scale and truncates.
func (g *ProtocolInfo) scaleRareVariantShareStat(stat, scale mpc_core.RVec) mpc_core.RVec {
	if len(stat) == 0 || len(scale) == 0 {
		return stat
	}
	mpcObj := g.mpcObj[0]
	return mpcObj.TruncVec(mpcObj.SSMultElemVecScalar(stat, scale[0]), mpcObj.GetDataBits(), mpcObj.GetFracBits())
}

func (g *ProtocolInfo) CZeroTest() {
	//pid := g.mpcObj[0].GetPid()
	//mpc := g.mpcObj[0]
	//params := g.cps

	//cm := crypto.CZeroMat(params, 2, 2)
	//M := mpc.Network.CollectiveDecryptMat(params, cm, -1)
	//log.LLvl1(crypto.PlaintextToDense(params, M, 3))

	//pid := g.mpcObj[0].GetPid()
	mpc := g.mpcObj[0]
	rtype := mpc.GetRType()

	a := rtype.Zero()
	for i := 0; i < 100; i++ {
		t := time.Now()
		_, _ = mpc.SqrtAndSqrtInverse(mpc_core.RVec{a}, false)
		log.LLvl1(time.Since(t))
	}

}

func (g *ProtocolInfo) ConversionTest() {

	pid := g.mpcObj[0].GetPid()
	params := g.cps
	//slots := params.GetSlots()

	x := make([]float64, 10)
	for i := range x {
		x[i] = float64(i + 1)
	}

	cv, _ := crypto.EncryptFloatVector(params, x)

	if pid > 0 {
		pv := g.mpcObj[0].Network.CollectiveDecryptVec(params, cv, 2)
		f := crypto.DecodeFloatVector(params, pv)[:5]
		log.LLvl1(f)
	}

	ss := g.mpcObj[0].CiphertextToSS(params, g.mpcObj[0].GetRType(), cv[0], 1, 10)
	log.LLvl1("rtype", g.mpcObj[0].GetRType().TypeID())

	ssRev := g.mpcObj[0].RevealSymVec(ss)
	for i := range ssRev {
		log.LLvl1("Conv reveal", ssRev[i].Float64(g.mpcObj[0].GetFracBits()))
	}

	_, ssInv := g.mpcObj[0].SqrtAndSqrtInverse(ss, false)

	ssRev = g.mpcObj[0].RevealSymVec(ssInv)
	for i := range ssRev {
		log.LLvl1("SqrtInv reveal", ssRev[i].Float64(g.mpcObj[0].GetFracBits()))
	}

	out := g.mpcObj[0].SStoCiphertext(params, ssInv)
	if pid > 0 {
		pv := g.mpcObj[0].Network.CollectiveDecryptVec(params, crypto.CipherVector{out}, -1)
		log.LLvl1(crypto.DecodeFloatVector(params, pv)[:5])
	}

}

func (g *ProtocolInfo) Test() {
	if g.mpcObj[0].GetPid() == 0 {
		return
	}

	params := g.cps
	slots := params.GetSlots()

	x := make([]float64, slots)
	for i := range x {
		x[i] = 1.0
	}
	cv, _ := crypto.EncryptFloatVector(params, x)
	d := cv[0].Value[0].Coeffs
	log.LLvl1(time.Now().Format(time.RFC3339), "Enc check2:", d[0][0], d[1][1], d[2][2])

	start := time.Now()
	//filename := "data/gwas-toy-block/geno_party1.0.bin"
	//gfs := NewGenoFileStream(filename, uint64(500), uint64(500), true)
	filename := "data/tmp.bin"
	gfs := NewGenoFileStream(filename, uint64(10), uint64(76600), true)

	out := make(crypto.CipherMatrix, 1)
	if g.mpcObj[0].GetPid() == 2 {
		out, _, _ = MatMult4Stream(params, crypto.CipherMatrix{cv}, gfs, 5, true, false, 0)

		d = out[0][0].Value[0].Coeffs
		log.LLvl1(time.Now().Format(time.RFC3339), "Out check:", d[0][0], d[1][1], d[2][2])
	}
	//out, _ := MatMult4Stream(params, crypto.CipherMatrix{cv}, g.genoBlocks[0], 5, true)

	log.LLvl1(time.Now().Format(time.RFC3339), "MatMult: elapsed time", time.Since(start))

	pv := g.mpcObj[0].Network.CollectiveDecryptVec(params, out[0], 2)
	f := crypto.DecodeFloatVector(params, pv)[:100]
	log.LLvl1(f)
}

func (g *ProtocolInfo) SyncAndTerminate(closeChannelFlag bool) {
	mainMPCObj := g.mpcObj[0]
	pid := mainMPCObj.GetPid()

	var dummy mpc_core.RElem = mainMPCObj.GetRType().Zero()
	if pid == 0 {
		for p := 1; p < mainMPCObj.GetNParty(); p++ {
			dummy = mainMPCObj.Network.ReceiveRElem(dummy, p)
			mainMPCObj.Network.SendRData(dummy, p)
		}
	} else {
		mainMPCObj.Network.SendRData(dummy, 0)
		dummy = mainMPCObj.Network.ReceiveRElem(dummy, 0)
	}

	if closeChannelFlag {
		// Close all threads
		for t := range g.mpcObj {
			g.mpcObj[t].Network.CloseAll()
		}
	}

}

func (g *ProtocolInfo) OutPath(filename string) string {
	return filepath.Join(g.config.OutDir, filename)
}

func (g *ProtocolInfo) CachePath(filename string) string {
	return filepath.Join(g.config.CacheDir, filename)
}

func (g *ProtocolInfo) GeneratePCAInput(numSnpsPCA int, snpFiltPCA []bool, isPgen bool) (*GenoFileStream, *GenoFileStream) {
	mergedFile := g.CachePath("geno_pca.bin")
	outTransFile := g.CachePath("geno_pca_transpose.bin")

	pid := g.mpcObj[0].GetPid()
	// numIndsPCA := int(g.genoBlocks[0].NumRowsToKeep())
	numIndsPCA := g.gwasParams.numFiltInds[pid]
	if pid > 0 {
		log.LLvl1(fmt.Sprintf("GeneratePCAInput: filtered local data has %d snps, %d samples (pid %d)", numSnpsPCA, numIndsPCA, pid))
	}

	if _, err := os.Stat(mergedFile); os.IsNotExist(err) {
		numSnpsPCAPerBlock := make([]int, g.config.GenoNumBlocks)
		if isPgen { // pgen format
			shift := 0
			for chr := 0; chr < g.config.GenoNumBlocks; chr++ { // Iterate over chromosomes
				pgenPrefix := fmt.Sprintf(g.config.GenoFilePrefix, chr+1)   // 1-based
				outFile := g.CachePath(fmt.Sprintf("geno_pca.%d.bin", chr)) // 0-based

				m := g.genoBlockSizes[chr]
				numSnpsPCAPerBlock[chr] = SumBool(snpFiltPCA[shift : shift+m])

				FilterMatrixFilePgen(pgenPrefix, numIndsPCA, numSnpsPCAPerBlock[chr], g.config.SampleKeepFile, g.config.SnpIdsFile, shift, snpFiltPCA[shift:shift+m], outFile, false)

				shift += m
			}
		} else { // blocks format
			indFilt := g.genoBlocks[0].RowFilt()
			if indFilt == nil {
				indFilt = OnesBool(int(g.genoBlocks[0].NumRows()))
			}

			shift := 0
			for i := range g.genoBlocks { // TODO parallelize
				genoFile := fmt.Sprintf("%s.%d.bin", g.config.GenoFilePrefix, i)
				outFile := g.CachePath(fmt.Sprintf("geno_pca.%d.bin", i))

				n, m := int(g.genoBlocks[i].NumRows()), int(g.genoBlocks[i].NumCols())

				FilterMatrixFile(genoFile, n, m, indFilt, snpFiltPCA[shift:shift+m], outFile)

				numSnpsPCAPerBlock[i] = SumBool(snpFiltPCA[shift : shift+m])

				shift += m
			}
		}

		MergeBlockFiles(g.CachePath("geno_pca"), numIndsPCA, numSnpsPCAPerBlock, mergedFile)
	} else {
		log.LLvl1("Cache file found:", mergedFile)
	}

	if _, err := os.Stat(outTransFile); os.IsNotExist(err) {
		TransposeMatrixFile(mergedFile, numIndsPCA, numSnpsPCA, outTransFile)
	} else {
		log.LLvl1("Cache file found:", outTransFile)
	}

	genoFs1 := NewGenoFileStream(mergedFile, uint64(numIndsPCA), uint64(numSnpsPCA), true)
	genoFs2 := NewGenoFileStream(outTransFile, uint64(numSnpsPCA), uint64(numIndsPCA), true)

	return genoFs1, genoFs2
}

func snpDistanceFiltering(snpPositions []uint64, snpFilt []bool, snpDistThreshold uint64) (int, []bool) {

	numSnpsPCA := 0
	prevPos := uint64(0)
	snpFiltPCA := make([]bool, len(snpFilt))

	for i := range snpFilt {
		if snpFilt[i] {
			if numSnpsPCA == 0 || snpPositions[i] >= prevPos+snpDistThreshold {
				snpFiltPCA[i] = true
				prevPos = snpPositions[i]
				numSnpsPCA++
			}
		}
	}

	return numSnpsPCA, snpFiltPCA
}

func (g *ProtocolInfo) PopulationStratification() crypto.CipherMatrix {
	params := g.gwasParams
	mpc := g.mpcObj[0]
	pid := mpc.GetPid()
	isPgen := g.IsPgen()

	if params.GetPopStratMethod() {

		log.LLvl1(time.Now().Format(time.RFC3339), "Starting distributed PCA routine")

		var genoReduced, genoReducedT *GenoFileStream
		var numSnpsPCA int
		var snpFiltPCA []bool
		if pid > 0 {
			snpFiltQC := OnesBool(params.numSnps)
			if !g.config.SkipQC || g.config.UseCachedQC {
				if isPgen {
					snpFiltQC = readFilterFromFile(g.CachePath("gkeep.txt"), params.numSnps, false)
				} else {
					// Concatenate snp filters
					shift := 0
					for i := range g.genoBlocks {
						m := int(g.genoBlocks[i].NumCols())
						copy(snpFiltQC[shift:shift+m], g.genoBlocks[i].ColFilt())
						shift += m
					}
				}
			}

			log.LLvl1(time.Now().Format(time.RFC3339), "SNP pruning for PCA")

			numSnpsPCA, snpFiltPCA = snpDistanceFiltering(g.pos, snpFiltQC, params.MinSnpDistThreshold())

			if pid == mpc.GetHubPid() {
				mpc.Network.SendInt(numSnpsPCA, 0)
			}

			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Number of SNPs selected for PCA: %d", numSnpsPCA))

			log.LLvl1(time.Now().Format(time.RFC3339), "Generating reduced input file for PCA")

			genoReduced, genoReducedT = g.GeneratePCAInput(numSnpsPCA, snpFiltPCA, isPgen)

		} else { // Party 0
			numSnpsPCA = mpc.Network.ReceiveInt(mpc.GetHubPid())

			log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Number of SNPs selected for PCA: %d", numSnpsPCA))
		}

		g.gwasParams.SetNumSnpsPCA(numSnpsPCA)

		pca := g.InitPCA(genoReduced, genoReducedT)

		start := time.Now()

		log.LLvl1(time.Now().Format(time.RFC3339), "AssertSync")
		g.mpcObj[0].AssertSync()

		pca.Q = pca.DistributedPCA()

		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Finished distributed PCA, %s", time.Since(start)))

		return pca.Q
	}

	//TODO: add LMM option here
	return make(crypto.CipherMatrix, g.gwasParams.NumPC())

}

func (g *ProtocolInfo) ComputeAssocStatistics(Qpca crypto.CipherMatrix) (crypto.CipherVector, []bool) {
	assocTest := g.InitAssociationTests(Qpca)
	return assocTest.GetAssociationStats()
}
