package gwas

import (
	"bufio"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"go.dedis.ch/onet/v3/log"

	"fmt"

	mpc_core "github.com/hhcho/mpc-core"

	"github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/examples"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/mat"
)

type ProtocolInfo struct {
	mpcObj mpc.ParallelMPC
	cps    *crypto.CryptoParams
	cpsPar []*crypto.CryptoParams // One per thread

	// Input files
	genoBlocks           []*GenoFileStream
	genoBlockSizes       []int
	hiddenGenoBlockSizes []int
	pheno                *mat.Dense
	cov                  *mat.Dense
	pos                  []uint64

	gwasParams *GWASParams

	config *Config
}

type RareVariantMode string

const (
	RareVariantModeGWAS   RareVariantMode = "gwas"
	RareVariantModeSKAT   RareVariantMode = "skat"
	RareVariantModeBurden RareVariantMode = "burden"
	RareVariantModeSKATO  RareVariantMode = "skato"
)

type RareVariantSetMode string

const (
	RareVariantSetModeShared         RareVariantSetMode = "shared"
	RareVariantSetModeMVPPublicUnion RareVariantSetMode = "mvp_public_union"
)

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

	RareVariantSetMode      string `toml:"rare_variant_set_mode"`
	PublicVariantPartyID    int    `toml:"public_variant_party_id"`
	PrivateVariantPartyID   int    `toml:"private_variant_party_id"`
	HiddenGenoFileFormat    string `toml:"hidden_geno_file_format"`
	HiddenGenoFilePrefix    string `toml:"hidden_geno_binary_file_prefix"`
	HiddenGenoNumBlocks     int    `toml:"hidden_geno_num_blocks"`
	HiddenGenoBlockSizeFile string `toml:"hidden_geno_block_size_file"`
	HiddenSnpIdsFile        string `toml:"hidden_snp_ids_file"`
	HiddenSnpPosFile        string `toml:"hidden_snp_position_file"`
	HiddenSampleKeepFile    string `toml:"hidden_sample_keep_file"`

	PhenoFile  string `toml:"pheno_file"`
	CovFile    string `toml:"covar_file"`
	SnpPosFile string `toml:"snp_position_file"`

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

	Debug          bool  `toml:"debug"`
	BlocksForAssoc []int `toml:"blocks_for_assoc_test"`
	PgenBatchSize  int   `toml:"pgen_batch_nsnp"`
}

func (config *Config) VariantSetMode() RareVariantSetMode {
	mode := strings.TrimSpace(strings.ToLower(config.RareVariantSetMode))
	if mode == "" {
		return RareVariantSetModeShared
	}
	return RareVariantSetMode(mode)
}

func (config *Config) IsMVPPublicUnionMode() bool {
	return config.VariantSetMode() == RareVariantSetModeMVPPublicUnion
}

func (prot *ProtocolInfo) VariantSetMode() RareVariantSetMode {
	return prot.config.VariantSetMode()
}

func (prot *ProtocolInfo) IsMVPPublicUnionMode() bool {
	return prot.config.IsMVPPublicUnionMode()
}

func readBlockSizesFile(filename string, expectedBlocks int) []int {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatalf("failed to open: %v", filename)
	}
	defer file.Close()

	var sizes []int
	scanner := bufio.NewScanner(file)
	scanner.Split(bufio.ScanLines)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}
		size, err := strconv.Atoi(line)
		if err != nil {
			log.Fatalf("parse error: %v", filename)
		}
		sizes = append(sizes, size)
	}
	if err := scanner.Err(); err != nil {
		log.Fatalf("scan error: %v", filename)
	}

	if expectedBlocks > 0 && len(sizes) != expectedBlocks {
		log.Fatalf("block count mismatch in %v: expected %d, got %d", filename, expectedBlocks, len(sizes))
	}
	return sizes
}

func readSampleKeepIIDs(filename string) []string {
	file, err := os.Open(filename)
	if err != nil {
		log.Fatalf("failed to open sample keep file: %v", filename)
	}
	defer file.Close()

	var out []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) == 0 {
			continue
		}
		if len(fields) >= 2 {
			out = append(out, fields[1])
		} else {
			out = append(out, fields[0])
		}
	}
	if err := scanner.Err(); err != nil {
		log.Fatalf("scan error: %v", filename)
	}
	return out
}

func assertSameSampleOrder(leftPath, rightPath string) {
	left := readSampleKeepIIDs(leftPath)
	right := readSampleKeepIIDs(rightPath)
	if len(left) != len(right) {
		log.Fatalf("sample keep order mismatch: %s has %d samples, %s has %d samples", leftPath, len(left), rightPath, len(right))
	}
	for i := range left {
		if left[i] != right[i] {
			log.Fatalf("sample keep order mismatch at row %d: %s has %q, %s has %q", i+1, leftPath, left[i], rightPath, right[i])
		}
	}
}

func validateMVPPublicUnionConfig(config *Config, pid int) {
	if !config.IsMVPPublicUnionMode() {
		return
	}
	if config.NumMainParties != 2 {
		log.Fatalf("mvp_public_union v1 supports exactly two data parties, got %d", config.NumMainParties)
	}
	if config.PrivateVariantPartyID != 1 || config.PublicVariantPartyID != 2 {
		log.Fatalf("mvp_public_union v1 requires private_variant_party_id=1 (AoU) and public_variant_party_id=2 (MVP)")
	}
	if config.HiddenGenoFileFormat != "pgen" {
		log.Fatalf("mvp_public_union v1 supports hidden_geno_file_format=\"pgen\" only")
	}
	if config.HiddenGenoFilePrefix == "" || config.HiddenGenoBlockSizeFile == "" ||
		config.HiddenSnpIdsFile == "" || config.HiddenSampleKeepFile == "" {
		log.Fatalf("mvp_public_union requires hidden genotype prefix, block size file, SNP ids file, and sample keep file")
	}
	if pid == config.PrivateVariantPartyID {
		assertSameSampleOrder(config.SampleKeepFile, config.HiddenSampleKeepFile)
	}
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

	if !mpcOnly {
		for thread := range networks {
			networks[thread].SetMHEParams(&params)
		}
	}

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
	if !mpcOnly {
		cps = networks.CollectiveInit(&params, prec)
	}

	var pheno, cov *mat.Dense
	var pos []uint64
	var genofs []*GenoFileStream
	var genoBlockSizes []int
	var hiddenGenoBlockSizes []int

	isPgen := config.GenoFileFormat == "pgen"
	validateMVPPublicUnionConfig(config, pid)

	genofs = make([]*GenoFileStream, config.GenoNumBlocks)
	genoBlockSizes = make([]int, config.GenoNumBlocks)

	if pid > 0 {
		// Read geno block size file
		genoBlockSizes = readBlockSizesFile(config.GenoBlockSizeFile, config.GenoNumBlocks)

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

	if config.IsMVPPublicUnionMode() {
		hiddenGenoBlockSizes = readBlockSizesFile(config.HiddenGenoBlockSizeFile, config.HiddenGenoNumBlocks)
		config.HiddenGenoNumBlocks = len(hiddenGenoBlockSizes)
	}

	gwasParams := InitGWASParams(config.NumInds, config.NumSnps, config.NumCovs, config.NumPCs, config.SnpDistThres)

	return &ProtocolInfo{
		mpcObj: mpcEnv, // One MPC object for each thread
		cps:    cps,

		genoBlocks:           genofs,
		genoBlockSizes:       genoBlockSizes,
		hiddenGenoBlockSizes: hiddenGenoBlockSizes,
		pheno:                pheno,
		cov:                  cov,
		pos:                  pos,

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

func (g *ProtocolInfo) decryptRareVariantScalar(stat crypto.CipherVector) []float64 {
	mpcObj := g.mpcObj[0]
	pid := mpcObj.GetPid()
	if pid == 0 {
		return nil
	}

	statDec := mpcObj.Network.CollectiveDecryptVec(g.cps, stat, -1)
	out := crypto.DecodeFloatVector(g.cps, statDec)
	return []float64{out[0]}
}

func (g *ProtocolInfo) rareVariantScaleShares(rssEnc crypto.CipherVector) (mpc_core.RVec, bool) {
	mpcObj := g.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()
	scaleSS := mpc_core.InitRVec(rtype.Zero(), 1)
	sourcePid := mpcObj.GetHubPid()

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

	var rssCt *rlwe.Ciphertext
	if pid == sourcePid && rssEnc != nil && len(rssEnc) > 0 {
		rssCt = rssEnc[0]
	}
	rssSS := mpcObj.CiphertextToSS(g.cps, rtype, rssCt, sourcePid, 1)

	numerSS := mpc_core.InitRVec(rtype.Zero(), 1)
	if pid == sourcePid {
		numerSS[0] = rtype.FromFloat64(float64(dof)/2.0, mpcObj.GetFracBits())
	}

	return mpcObj.Divide(numerSS, rssSS, false), true
}

func (g *ProtocolInfo) scaleRareVariantShareStat(stat, scale mpc_core.RVec) mpc_core.RVec {
	if len(stat) == 0 || len(scale) == 0 {
		return stat
	}
	mpcObj := g.mpcObj[0]
	out := mpcObj.SSMultElemVec(stat, scale)
	return mpcObj.TruncVec(out, mpcObj.GetDataBits(), mpcObj.GetFracBits())
}

func (g *ProtocolInfo) saveRareVariantScalar(filename string, stat crypto.CipherVector) {
	if out := g.decryptRareVariantScalar(stat); out != nil {
		SaveFloatVectorToFile(g.OutPath(filename), out)
	}
}

func (g *ProtocolInfo) combineSKATOStatistic(skat, burden crypto.CipherVector, rho float64) crypto.CipherVector {
	if g.mpcObj[0].GetPid() == 0 || skat == nil || burden == nil || len(skat) == 0 || len(burden) == 0 || skat[0] == nil || burden[0] == nil {
		return crypto.CipherVector{nil}
	}

	cryptoParams := g.cps

	skatPart := crypto.CMultConstRescale(cryptoParams, skat, 1.0-rho, false)
	burdenPart := crypto.CMultConstRescale(cryptoParams, burden, rho, false)
	return crypto.CAdd(cryptoParams, skatPart, burdenPart)
}

func (g *ProtocolInfo) RunRareVariantTest(mode RareVariantMode, skatoRho float64) {
	log.LLvl1(time.Now().Format(time.RFC3339), "Running rare-variant protocol in mode:", string(mode))

	// Rare-variant tests still rely on QC-filtered sample/SNP counts and filters.
	// Reuse the standard QC phase to populate gwasParams before SKAT/Burden/SKAT-O.
	g.Phase1()

	skat, burden := g.ComputeSKATStatistics()

	log.LLvl1(time.Now().Format(time.RFC3339), "Finished rare-variant statistic computation")

	net := g.mpcObj.GetNetworks()
	net.PrintNetworkLog()

	switch mode {
	case RareVariantModeSKAT:
		g.saveRareVariantScalar("skat_out.txt", skat)
		g.saveRareVariantScalar("burden_out.txt", burden)
		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Output collectively decrypted and saved to: %s and %s", g.OutPath("skat_out.txt"), g.OutPath("burden_out.txt")))
	case RareVariantModeBurden:
		g.saveRareVariantScalar("burden_out.txt", burden)
		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Output collectively decrypted and saved to: %s", g.OutPath("burden_out.txt")))
	case RareVariantModeSKATO:
		skato := g.combineSKATOStatistic(skat, burden, skatoRho)
		g.saveRareVariantScalar("skat_out.txt", skat)
		g.saveRareVariantScalar("burden_out.txt", burden)
		g.saveRareVariantScalar("skato_out.txt", skato)
		log.LLvl1(time.Now().Format(time.RFC3339), fmt.Sprintf("Output collectively decrypted and saved to: %s, %s, and %s", g.OutPath("skat_out.txt"), g.OutPath("burden_out.txt"), g.OutPath("skato_out.txt")))
	default:
		panic(fmt.Sprintf("unsupported rare-variant mode: %s", mode))
	}
}

// SKAT is the main entry point to run the Secure SKAT protocol independently of single-variant GWAS.
func (g *ProtocolInfo) SKAT() {
	g.RunRareVariantTest(RareVariantModeSKAT, 0.0)
}

func (g *ProtocolInfo) Burden() {
	g.RunRareVariantTest(RareVariantModeBurden, 0.0)
}

func (g *ProtocolInfo) SKATO(rho float64) {
	g.RunRareVariantTest(RareVariantModeSKATO, rho)
}

// SetPhenoAndCov allows deterministic testing scripts to override File I/O datasets in memory
func (g *ProtocolInfo) SetPhenoAndCov(pheno, cov *mat.Dense) {
	g.pheno = pheno
	g.cov = cov
}

func (g *ProtocolInfo) ComputeSKATStatistics() (crypto.CipherVector, crypto.CipherVector) {
	assocTest := g.InitAssociationTests(nil) // SKAT does not use PCA
	return assocTest.ComputeSKATStatistics()
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
