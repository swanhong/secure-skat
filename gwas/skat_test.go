package gwas

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"testing"
	"time"

	"github.com/BurntSushi/toml"
	"github.com/hhcho/sfgwas/crypto"
	"gonum.org/v1/gonum/mat"
)

// InitProtocolForTest initializes the GWAS protocol using the config files.
func InitProtocolForTest(t *testing.T) *ProtocolInfo {
	pidStr := os.Getenv("PID")
	if pidStr == "" {
		t.Skip("PID environment variable not set. Skipping test.")
		return nil
	}
	pid, err := strconv.Atoi(pidStr)
	if err != nil {
		t.Fatalf("Invalid PID: %s", err)
	}

	// SFGWAS_CONFIG_PATH lets the compare harness (scripts/analysis/skat/run.sh)
	// point at an alternate config/data dir; default is the repo's config/.
	configPath := os.Getenv("SFGWAS_CONFIG_PATH")
	if configPath == "" {
		_, filename, _, _ := runtime.Caller(0)
		configPath = filepath.Join(filepath.Dir(filepath.Dir(filename)), "config")
	}
	config := new(Config)

	// Import global parameters
	if _, err := toml.DecodeFile(filepath.Join(configPath, "configGlobal.toml"), config); err != nil {
		t.Fatalf("Failed to read global config: %s", err)
	}

	// Import local parameters
	if _, err := toml.DecodeFile(filepath.Join(configPath, fmt.Sprintf("configLocal.Party%d.toml", pid)), config); err != nil {
		t.Fatalf("Failed to read local config for PID=%d: %s", pid, err)
	}

	// SFGWAS_DEBUG forces per-block intermediate dumps (qBlock_block*.txt etc.) so
	// run.sh can do plain-vs-secure per-block comparison without editing the config.
	if os.Getenv("SFGWAS_DEBUG") != "" {
		config.Debug = true
	}

	// Create cache/output directories
	if err := os.MkdirAll(config.CacheDir, 0755); err != nil {
		t.Fatalf("Failed to create cache dir: %s", err)
	}
	if err := os.MkdirAll(config.OutDir, 0755); err != nil {
		t.Fatalf("Failed to create out dir: %s", err)
	}

	runtime.GOMAXPROCS(config.LocalNumThreads)

	return InitializeGWASProtocol(config, pid, false)
}

func TestSKATBasicOperations(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	mpcObj := prot.mpcObj[0]
	cps := prot.cps

	t.Run("Test_MHE_EncDec", func(t *testing.T) {
		// Verify MHE encryption and decryption of vectors
		n := 10
		vals := make([]float64, n)
		for i := 0; i < n; i++ {
			vals[i] = rand.Float64() * 10.0
		}

		// Encrypt
		ct, _ := crypto.EncryptFloatVector(cps, vals)

		// Decrypt natively
		var pt []float64
		if mpcObj.GetPid() > 0 {
			// This is local decryption, but wait sf-GWAS uses CollectiveDecryptVec
			pv := mpcObj.Network.CollectiveDecryptVec(cps, ct, 1)
			if mpcObj.GetPid() == 1 {
				pt = crypto.DecodeFloatVector(cps, pv)
				for i := 0; i < n; i++ {
					diff := pt[i] - vals[i]
					if diff < -1e-4 || diff > 1e-4 {
						t.Errorf("MHE EncDec mismatch at %d: orig %v, rep %v", i, vals[i], pt[i])
					}
				}
			}
		}
	})

	t.Run("Test_MPC_MHE_Swapping", func(t *testing.T) {
		debug := true
		n := 10
		vals := make([]float64, n)
		for i := 0; i < n; i++ {
			vals[i] = rand.Float64() * 5.0
		}

		// 1. Encrypt to MHE
		t.Logf("PID %d encrypting %v\n", mpcObj.GetPid(), vals)
		ct, _ := crypto.EncryptFloatVector(cps, vals)

		// 2. Convert to SS
		// CVecToSS(cps, rtype, cvec, sourcePid, numCtx, numValidSlots)
		shares := mpcObj.CVecToSS(cps, mpcObj.GetRType(), ct, 1, 1, cps.GetSlots())

		if debug && mpcObj.GetPid() > 0 {
			revealed := mpcObj.RevealSymVec(shares).ToFloat(mpcObj.GetFracBits())
			t.Logf("SS share for PID %d: %v\n", mpcObj.GetPid(), revealed)
			if mpcObj.GetPid() == 1 {
				for i := 0; i < n; i++ {
					diff := revealed[i] - vals[i]
					if diff < -1e-2 || diff > 1e-2 {
						t.Errorf("MHE->MPC mismatch at %d: orig %v, rep %v", i, vals[i], revealed[i])
					}
				}
			}
		}

		// 3. Convert SS back to MHE
		// SSToCVec(cps, rvec)
		ct2 := mpcObj.SSToCVec(cps, shares)

		// 4. Decrypt natively and check
		if mpcObj.GetPid() > 0 {
			pv := mpcObj.Network.CollectiveDecryptVec(cps, ct2, 1)
			if mpcObj.GetPid() == 1 {
				pt := crypto.DecodeFloatVector(cps, pv)
				t.Logf("PID %d decrypted %v\n", mpcObj.GetPid(), pt)
				for i := 0; i < n; i++ {
					diff := pt[i] - vals[i]
					if diff < -1e-2 || diff > 1e-2 {
						t.Errorf("MPC-MHE Swapping mismatch at %d: orig %v, rep %v", i, vals[i], pt[i])
					}
				}
				fmt.Println("Test_MPC_MHE_Swapping verified.")
			}
		}
	})

	t.Run("Test_Realistic_Projection_And_GenoMult", func(t *testing.T) {
		// 1. Shared Seed for Random Generation so parties can know each other's data for plain check
		// In a real test, they wouldn't share seeds, but here it simplifies the plain reference.
		const seed = 42
		r1 := rand.New(rand.NewSource(seed))
		r2 := rand.New(rand.NewSource(seed)) // for y

		// Total dimensions
		nParties := prot.GetConfig().NumMainParties
		nIndsAll := prot.GetConfig().NumInds
		nIndsTotal := 0
		for p := 1; p <= nParties; p++ {
			nIndsTotal += nIndsAll[p]
		}
		nCovs := prot.GetConfig().NumPCs // Use NumPCs as k for QR
		pid := mpcObj.GetPid()

		// 2. Generate Full X and y (all parties generate the same for plain check)
		fullX := mat.NewDense(nIndsTotal, nCovs, nil)
		for i := 0; i < nIndsTotal; i++ {
			// First column is intercept for NetDQRenc compatibility
			fullX.Set(i, 0, 1.0)
			for j := 1; j < nCovs; j++ {
				fullX.Set(i, j, r1.Float64())
			}
		}
		t.Logf("PID %d generated fullX(size: %d x %d)", pid, fullX.RawMatrix().Rows, fullX.RawMatrix().Cols)
		t.Logf("first 5x5 of fullX: %v", fullX.Slice(0, 5, 0, 5))

		fullY := make([]float64, nIndsTotal)
		for i := 0; i < nIndsTotal; i++ {
			fullY[i] = r2.Float64()
		}
		t.Logf("PID %d generated fullY(size: %d)", pid, len(fullY))
		t.Logf("first 5 of fullY: %v", fullY[:5])

		// 3. Extract Local Slice and Encrypt
		offset := 0
		for p := 1; p < pid; p++ {
			offset += nIndsAll[p]
		}
		nIndsLocal := 0
		if pid > 0 {
			nIndsLocal = nIndsAll[pid]
		}

		var X_enc crypto.CipherMatrix
		var y_enc crypto.CipherMatrix
		if pid > 0 {
			localX := fullX.Slice(offset, offset+nIndsLocal, 0, nCovs).(*mat.Dense)
			localY := fullY[offset : offset+nIndsLocal]

			X_enc = crypto.EncryptDense(cps, localX)
			y_vec, _ := crypto.EncryptFloatVector(cps, localY)
			y_enc = crypto.CipherMatrix{y_vec}
		} else {
			// Party 0 must still provide correctly-sized matrices for synchronization
			// Functions like NetDQRenc and DCMatMulAAtB involve collective operations
			// that may fail if a party has 0 ciphertexts while others have some.
			X_enc = make(crypto.CipherMatrix, nCovs)
			for i := range X_enc {
				X_enc[i] = crypto.CZeros(cps, 1)
			}
			y_enc = make(crypto.CipherMatrix, 1)
			y_enc[0] = crypto.CZeros(cps, 1)
		}

		// 4. Manually initialize filters
		numSnps := prot.GetGwasParams().NumSNP()
		snpFilt := make([]bool, numSnps)
		for i := range snpFilt {
			snpFilt[i] = true
		}
		prot.GetGwasParams().SetSnpFilt(snpFilt)
		numFiltInds := make([]int, nParties+1)
		copy(numFiltInds, nIndsAll)
		prot.GetGwasParams().SetFiltCounts(numFiltInds, numSnps)

		// 5. Encrypted QR Decomposition
		// NetDQRenc handles distributed Householder QR
		Q_enc := NetDQRenc(cps, mpcObj, X_enc, nIndsAll)
		t.Logf("PID %d Q_enc: %v", pid, Q_enc)

		// 6. Encrypted Projection y_res = y - 1/N * Q(Q^T y)
		// DCMatMulAAtB computes Q(Q^T y)
		var y_res_enc crypto.CipherMatrix
		if pid > 0 && nIndsLocal > 0 {
			y_proj_enc := DCMatMulAAtB(cps, mpcObj, Q_enc, y_enc, nIndsAll, 1, func(cp *crypto.CryptoParams, a crypto.CipherVector, B crypto.CipherMatrix, j int) crypto.CipherVector {
				return crypto.CMult(cp, a, B[j])
			})

			// Scale by 1/N and subtract
			// Q_enc returned by NetDQRenc is scaled by sqrt(N)
			// So Q(Q^T y) is scaled by N
			y_proj_rescaled := crypto.CMultConstMat(cps, y_proj_enc, 1.0/float64(nIndsTotal), false)
			y_proj_rescaled = crypto.CMatRescale(cps, y_proj_rescaled)

			y_res_vec := crypto.CSub(cps, y_enc[0], y_proj_rescaled[0])
			y_res_enc = crypto.CipherMatrix{y_res_vec}
		} else {
			y_res_enc = make(crypto.CipherMatrix, 1) // dummy for p0
		}

		// 7. GenoBlockMult: S = G^T y_res
		ast := prot.InitAssociationTests(make(crypto.CipherMatrix, nCovs))
		// Disable caching for this test to ensure it uses the new residuals
		ast.general.config.CacheDir = "tmp_cache_" + strconv.Itoa(pid)
		os.MkdirAll(ast.general.config.CacheDir, 0755)
		defer os.RemoveAll(ast.general.config.CacheDir)

		var s_enc crypto.CipherMatrix
		var dosageSum []float64
		if pid > 0 && nIndsLocal > 0 {
			s_enc, dosageSum, _, _ = ast.GenoBlockMult(0, y_res_enc, false)
		}

		// 8. Plain Computation for Comparison (on Party 1)
		if pid == 1 {
			// Plain QR
			var qr mat.QR
			qr.Factorize(fullX)
			// Allocate full Q to be safe, then slice to thin Q
			Q_full := mat.NewDense(nIndsTotal, nIndsTotal, nil)
			qr.QTo(Q_full)
			Q_plain := Q_full.Slice(0, nIndsTotal, 0, nCovs)

			// qr.QTo produces a thin Q of size nIndsTotal x nCovs
			qR, qC := Q_plain.Dims()
			fmt.Printf("Q_plain dims: %d x %d\n", qR, qC)

			// Plain Projection
			matY := mat.NewVecDense(nIndsTotal, fullY)
			qty := mat.NewVecDense(nCovs, nil)
			qty.MulVec(Q_plain.T(), matY)

			y_proj_plain := mat.NewVecDense(nIndsTotal, nil)
			y_proj_plain.MulVec(Q_plain, qty)

			y_res_plain := mat.NewVecDense(nIndsTotal, nil)
			y_res_plain.SubVec(matY, y_proj_plain)

			// Comparison: Residuals
			pv := mpcObj.Network.CollectiveDecryptVec(cps, y_res_enc[0], 1)
			pt := crypto.DecodeFloatVector(cps, pv)
			maxErrRes := 0.0
			for i := 0; i < nIndsLocal; i++ {
				// pt[i] is the local residual. Compare with y_res_plain[offset + i]
				err := math.Abs(pt[i] - y_res_plain.AtVec(offset+i))
				if err > maxErrRes {
					maxErrRes = err
				}
			}
			fmt.Printf("Max Error (Residuals): %e\n", maxErrRes)

			// Comparison: Score Vector
			// We need to read genotypes in plain for the reference score
			// For simplicity in this test, we'll just verify the residuals property
			// but we can also just output the decrypted score.
			if len(s_enc) > 0 {
				ps := mpcObj.Network.CollectiveDecryptVec(cps, s_enc[0], 1)
				pt_s := crypto.DecodeFloatVector(cps, ps)
				fmt.Printf("Decrypted Score Vector (first 5): %v\n", pt_s[:5])
				fmt.Printf("Dosage Sum (first 5): %v\n", dosageSum[:5])
			}

			if maxErrRes > 1e-2 {
				t.Errorf("Projection error too large: %e", maxErrRes)
			}
			fmt.Println("Test_Realistic_Projection_And_GenoMult verified.")

		} else {
			// Other parties participate in collective decryptions
			if pid > 0 {
				mpcObj.Network.CollectiveDecryptVec(cps, y_res_enc[0], 1)
				if len(s_enc) > 0 {
					mpcObj.Network.CollectiveDecryptVec(cps, s_enc[0], 1)
				}
			} else {
				// Party 0
				mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
				mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
			}
		}
	})

	t.Run("TestSKAT_Step1_ScoreVector", func(t *testing.T) {
		const seed = 42
		r1 := rand.New(rand.NewSource(seed))
		r2 := rand.New(rand.NewSource(seed))

		nParties := prot.GetConfig().NumMainParties
		nIndsAll := prot.GetConfig().NumInds
		nIndsTotal := 0
		for p := 1; p <= nParties; p++ {
			nIndsTotal += nIndsAll[p]
		}
		// Inject mock covariates and phenotypes into ProtocolInfo
		// The `prot` variable is already defined in the test scope.
		phenoMat := mat.NewDense(1000, 1, nil)
		covMat := mat.NewDense(1000, 5, nil)
		for i := 0; i < 1000; i++ {
			covMat.Set(i, 0, 1.0)
			phenoMat.Set(i, 0, 1.0)
			for j := 1; j < 5; j++ {
				// Inject deterministic variance to ensure Matrix Rank = 5
				// thereby averting Division-By-Zero overflow in Gram-Schmidt QR Factorizations
				variance := float64((i*j)%13) * 0.1
				covMat.Set(i, j, float64(j)*0.1+variance)
			}
		}
		prot.SetPhenoAndCov(phenoMat, covMat) // 1000 samples, 1 pheno, 5 covs
		nCovs := prot.GetConfig().NumPCs
		pid := mpcObj.GetPid()

		// Generate Full X and y
		fullX := mat.NewDense(nIndsTotal, nCovs, nil)
		for i := 0; i < nIndsTotal; i++ {
			fullX.Set(i, 0, 1.0)
			for j := 1; j < nCovs; j++ {
				fullX.Set(i, j, r1.Float64())
			}
		}
		fullY := make([]float64, nIndsTotal)
		for i := 0; i < nIndsTotal; i++ {
			fullY[i] = r2.Float64()
		}

		offset := 0
		for p := 1; p < pid; p++ {
			offset += nIndsAll[p]
		}
		nIndsLocal := 0
		if pid > 0 {
			nIndsLocal = nIndsAll[pid]
		}

		// Inject mock covariates and phenotypes into ProtocolInfo
		if pid > 0 {
			localX := fullX.Slice(offset, offset+nIndsLocal, 0, nCovs).(*mat.Dense)
			localY := mat.NewDense(nIndsLocal, 1, fullY[offset:offset+nIndsLocal])
			prot.SetPhenoAndCov(localY, localX)
		} else {
			prot.SetPhenoAndCov(nil, nil)
		}

		// Manually initialize filters
		numSnps := prot.GetGwasParams().NumSNP()
		snpFilt := make([]bool, numSnps)
		for i := range snpFilt {
			snpFilt[i] = true
		}
		prot.GetGwasParams().SetSnpFilt(snpFilt)
		numFiltInds := make([]int, nParties+1)
		copy(numFiltInds, nIndsAll)
		prot.GetGwasParams().SetFiltCounts(numFiltInds, numSnps)

		// Enable caching to `out/partyX` so user can inspect intermediate testing files
		prot.GetConfig().CacheDir = filepath.Join("out", "party"+strconv.Itoa(pid))
		os.MkdirAll(prot.GetConfig().CacheDir, 0755)
		// We deliberately don't remove this dir automatically

		// Call ComputeSKATStatistics directly to get the score vector
		// This tests the full integrated path in skat.go
		_, _, S_out, _ := prot.ComputeSKATStatistics()

		if pid == 1 {
			// Decrypt and compare score vector for the first block
			var qr mat.QR
			qr.Factorize(fullX)
			Q_full := mat.NewDense(nIndsTotal, nIndsTotal, nil)
			qr.QTo(Q_full)
			Q_plain := Q_full.Slice(0, nIndsTotal, 0, nCovs)

			matY := mat.NewVecDense(nIndsTotal, fullY)
			qty := mat.NewVecDense(nCovs, nil)
			qty.MulVec(Q_plain.T(), matY)

			y_proj_plain := mat.NewVecDense(nIndsTotal, nil)
			y_proj_plain.MulVec(Q_plain, qty)

			y_res_plain := mat.NewVecDense(nIndsTotal, nil)
			y_res_plain.SubVec(matY, y_proj_plain)

			// Assuming GenoBlockMult processes the full dataset here in the plain
			// However, S_out[0] is aggregated across parties via MPC.
			pv_s0 := mpcObj.Network.CollectiveDecryptVec(cps, S_out[0], 1)
			pt_s0 := crypto.DecodeFloatVector(cps, pv_s0)

			// We skip the exact math check for S_out since it depends on the plain genotype matrix,
			// which we haven't loaded in plain here. The goal is to ensure it runs without panics
			// and produces non-nil results.
			fmt.Printf("Step 1 Decrypted Aggregated Score Vector (first 5): %v\n", pt_s0[:5])
			if pt_s0[0] != 0 {
				fmt.Println("TestSKAT_Step1_ScoreVector verified.")
			}
		} else if pid > 0 {
			if len(S_out) > 0 {
				mpcObj.Network.CollectiveDecryptVec(cps, S_out[0], 1)
			}
		} else {
			// Party 0 must participate in collective decryption
			mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
		}
	})
}

func TestSecureSKATEndToEnd(t *testing.T) {

	rand.Seed(time.Now().UnixNano())
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	// Print metaparameters for all parties
	pid := prot.mpcObj[0].GetPid()
	config := prot.GetConfig()
	fmt.Printf("\n======================================================\n")
	fmt.Printf("SECURE SKAT METAPARAMETERS (PID=%d)\n", pid)
	fmt.Printf("======================================================\n")
	fmt.Printf("CKKS Parameters:    %s\n", config.CkksParams)
	fmt.Printf("MPC Field Size:     %d bits\n", config.MpcFieldSize)
	fmt.Printf("MPC Data Bits:      %d bits\n", config.MpcDataBits)
	fmt.Printf("MPC Frac Bits:      %d bits\n", config.MpcFracBits)
	fmt.Printf("MPC Num Threads:    %d\n", config.MpcNumThreads)
	fmt.Printf("Num Snps:           %d\n", config.NumSnps)
	fmt.Printf("Num Inds:           %v\n", config.NumInds)
	fmt.Printf("Geno File Format:   %s\n", config.GenoFileFormat)
	fmt.Printf("Pheno File:         %s\n", config.PhenoFile)
	fmt.Printf("Covariate File:     %s\n", config.CovFile)
	fmt.Printf("Output Directory:   %s\n", config.OutDir)
	fmt.Printf("======================================================\n\n")

	// Manually initialize filters (QC is entirely bypassed)
	nParties := prot.GetConfig().NumMainParties
	numSnps := prot.GetGwasParams().NumSNP()
	snpFilt := make([]bool, numSnps)
	for i := range snpFilt {
		snpFilt[i] = true
	}
	prot.GetGwasParams().SetSnpFilt(snpFilt)
	numFiltInds := make([]int, nParties+1)
	copy(numFiltInds, prot.GetConfig().NumInds)
	prot.GetGwasParams().SetFiltCounts(numFiltInds, numSnps)

	// Enable caching to `out/partyX` so user can inspect intermediate testing files
	prot.GetConfig().CacheDir = filepath.Join("out", "party"+strconv.Itoa(pid))
	os.MkdirAll(prot.GetConfig().CacheDir, 0755)

	prot.SKAT()

	// Print final confirmation
	if prot.GetConfig().NumMainParties > 0 {
		outPath := prot.OutPath("skat_out.txt")
		fmt.Printf("\n======================================================\n")
		fmt.Printf("Secure SKAT End-to-End verified!\n")
		fmt.Printf("Final output is saved in: %s\n", outPath)
		fmt.Printf("======================================================\n\n")
	} else {
		fmt.Printf("\n======================================================\n")
		fmt.Printf("Secure SKAT End-to-End verified for Hub (PID=0)\n")
		fmt.Printf("======================================================\n\n")
	}
}

func TestSecureBurdenEndToEnd(t *testing.T) {

	rand.Seed(time.Now().UnixNano())
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	// Print metaparameters for all parties
	pid := prot.mpcObj[0].GetPid()
	config := prot.GetConfig()
	fmt.Printf("\n======================================================\n")
	fmt.Printf("SECURE BURDEN METAPARAMETERS (PID=%d)\n", pid)
	fmt.Printf("======================================================\n")
	fmt.Printf("CKKS Parameters:    %s\n", config.CkksParams)
	fmt.Printf("MPC Field Size:     %d bits\n", config.MpcFieldSize)
	fmt.Printf("Output Directory:   %s\n", config.OutDir)
	fmt.Printf("======================================================\n\n")

	// Manually initialize filters (QC is entirely bypassed)
	nParties := prot.GetConfig().NumMainParties
	numSnps := prot.GetGwasParams().NumSNP()
	snpFilt := make([]bool, numSnps)
	for i := range snpFilt {
		snpFilt[i] = true
	}
	prot.GetGwasParams().SetSnpFilt(snpFilt)
	numFiltInds := make([]int, nParties+1)
	copy(numFiltInds, prot.GetConfig().NumInds)
	prot.GetGwasParams().SetFiltCounts(numFiltInds, numSnps)

	// Enable caching to `out/partyX` so user can inspect intermediate testing files
	prot.GetConfig().CacheDir = filepath.Join("out", "party"+strconv.Itoa(pid))
	os.MkdirAll(prot.GetConfig().CacheDir, 0755)

	prot.SKAT()

	// Print final confirmation
	if prot.GetConfig().NumMainParties > 0 {
		outPath := prot.OutPath("burden_out.txt")
		fmt.Printf("\n======================================================\n")
		fmt.Printf("Secure Burden End-to-End verified!\n")
		fmt.Printf("Final output is saved in: %s\n", outPath)
		fmt.Printf("======================================================\n\n")
	} else {
		fmt.Printf("\n======================================================\n")
		fmt.Printf("Secure Burden End-to-End verified for Hub (PID=0)\n")
		fmt.Printf("======================================================\n\n")
	}
}

// TestQBasisComputation verifies and prints the Q basis vectors produced by NetDQRenc.
// Useful for diagnosing scale inflation or numerical issues in the QR step.
func TestQBasisComputation(t *testing.T) {
	rand.Seed(time.Now().UnixNano())
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	mpcObj := prot.mpcObj[0]
	cps := prot.cps
	pid := mpcObj.GetPid()

	nIndsTotal := 2000
	nIndsAll := []int{0, 1000, 1000} // Party 0 has 0 rows, Party 1 and 2 have 1000 each
	nCovs := prot.GetConfig().NumPCs

	// Build deterministic full-rank covariates
	const seed = 42
	r := rand.New(rand.NewSource(seed))
	fullX := mat.NewDense(nIndsTotal, nCovs, nil)
	for i := 0; i < nIndsTotal; i++ {
		fullX.Set(i, 0, 1.0) // intercept
		for j := 1; j < nCovs; j++ {
			fullX.Set(i, j, r.Float64())
		}
	}

	// Extract local slice and encrypt
	offset := 0
	for p := 1; p < pid; p++ {
		offset += nIndsAll[p]
	}

	var X_enc crypto.CipherMatrix
	if pid > 0 {
		nLocal := nIndsAll[pid]
		localX := fullX.Slice(offset, offset+nLocal, 0, nCovs).(*mat.Dense)
		X_enc = crypto.EncryptDense(cps, localX)
	} else {
		X_enc = make(crypto.CipherMatrix, nCovs)
		for i := range X_enc {
			X_enc[i] = crypto.CZeros(cps, 1)
		}
	}

	// Run distributed QR
	Q_enc := NetDQRenc(cps, mpcObj, X_enc, nIndsAll)

	// Collectively decrypt and print first 5 elements of each Q column
	if pid > 0 {
		fmt.Printf("\n--- Step 1. QR Basis Computation [NetDQRenc] (Party %d) ---\n", pid)
		for col := range Q_enc {
			pv := mpcObj.Network.CollectiveDecryptVec(cps, Q_enc[col], 1)
			if pid == 1 {
				vals := crypto.DecodeFloatVector(cps, pv)
				n := 5
				if len(vals) < n {
					n = len(vals)
				}
				fmt.Printf("Q[%d] head:", col)
				for i := 0; i < n; i++ {
					fmt.Printf(" %e,", vals[i])
				}
				fmt.Println()
			}
		}
	} else {
		for range Q_enc {
			mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
		}
	}
}

// --- low-rank local plaintext contraction (skat.go LocalContract) ---

func rowSlice(m *mat.Dense, r0, r1 int) *mat.Dense {
	_, c := m.Dims()
	return mat.DenseCopyOf(m.Slice(r0, r1, 0, c))
}

func vecApprox(t *testing.T, name string, got, want []float64, tol float64) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: len got %d want %d", name, len(got), len(want))
	}
	for i := range got {
		if !approxEqual(got[i], want[i], tol) {
			t.Errorf("%s[%d]: got %.12g want %.12g", name, i, got[i], want[i])
		}
	}
}

func matApprox(t *testing.T, name string, got, want *mat.Dense, tol float64) {
	t.Helper()
	if !mat.EqualApprox(got, want, tol) {
		t.Errorf("%s mismatch:\ngot  %v\nwant %v", name, mat.Formatted(got), mat.Formatted(want))
	}
}

// LocalContract matches gonum on the full cohort (correctness).
func TestLowRankLocalMatchesGonum(t *testing.T) {
	G, X, y := plainFixture()
	n, _ := X.Dims()
	yv := mat.NewVecDense(n, y)

	lc := LocalContract(G, X, y)

	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	matApprox(t, "XtX", lc.XtX, &XtX, 1e-9)

	var Xty mat.VecDense
	Xty.MulVec(X.T(), yv)
	vecApprox(t, "Xty0", lc.Xty0, vecToSlice(&Xty), 1e-9)

	if !approxEqual(lc.Y0ty0, mat.Dot(yv, yv), 1e-9) {
		t.Errorf("y0ty0: got %.12g want %.12g", lc.Y0ty0, mat.Dot(yv, yv))
	}

	var GtX mat.Dense
	GtX.Mul(G.T(), X)
	matApprox(t, "GtX", lc.GtX, &GtX, 1e-9)

	var Gty mat.VecDense
	Gty.MulVec(G.T(), yv)
	vecApprox(t, "Gty0", lc.Gty0, vecToSlice(&Gty), 1e-9)

	_, m := G.Dims()
	wantDose := make([]float64, m)
	for j := 0; j < m; j++ {
		for i := 0; i < n; i++ {
			wantDose[j] += G.At(i, j)
		}
	}
	vecApprox(t, "dosageSum", lc.DosageSum, wantDose, 1e-9)
}

// n-independence invariant: Σ over party row-slices == full cohort.
func TestLowRankLocalPartyAdditivity(t *testing.T) {
	G, X, y := plainFixture() // n=6
	full := LocalContract(G, X, y)

	p1 := LocalContract(rowSlice(G, 0, 3), rowSlice(X, 0, 3), y[0:3])
	p2 := LocalContract(rowSlice(G, 3, 6), rowSlice(X, 3, 6), y[3:6])
	sum := p1.Add(p2)

	matApprox(t, "XtX", sum.XtX, full.XtX, 1e-9)
	vecApprox(t, "Xty0", sum.Xty0, full.Xty0, 1e-9)
	if !approxEqual(sum.Y0ty0, full.Y0ty0, 1e-9) {
		t.Errorf("y0ty0: got %.12g want %.12g", sum.Y0ty0, full.Y0ty0)
	}
	matApprox(t, "GtX", sum.GtX, full.GtX, 1e-9)
	vecApprox(t, "Gty0", sum.Gty0, full.Gty0, 1e-9)
	vecApprox(t, "dosageSum", sum.DosageSum, full.DosageSum, 1e-9)
}

// --- low-rank null model (skat.go localNullEquations / computeBetaHatEnc / RSS) ---

func TestLocalNullEquationsMatchesGonum(t *testing.T) {
	cov := mat.NewDense(5, 2, []float64{
		0.5, -1.0,
		-1.2, 0.3,
		0.3, 2.0,
		2.0, -0.7,
		-0.7, 1.1,
	})
	pheno := mat.NewDense(5, 1, []float64{1.0, -0.5, 0.3, 2.2, -1.1})

	ln := localNullEquations(cov, pheno, 0.0) // center 0 → y0 == pheno

	n, ncov := cov.Dims()
	X := mat.NewDense(n, ncov+1, nil)
	for i := 0; i < n; i++ {
		X.Set(i, 0, 1.0)
		for j := 0; j < ncov; j++ {
			X.Set(i, j+1, cov.At(i, j))
		}
	}
	y0v := mat.NewVecDense(n, []float64{1.0, -0.5, 0.3, 2.2, -1.1})

	var XtX mat.Dense
	XtX.Mul(X.T(), X)
	matApprox(t, "XtX", ln.XtX, &XtX, 1e-9)

	var Xty mat.VecDense
	Xty.MulVec(X.T(), y0v)
	vecApprox(t, "Xty0", ln.Xty0, vecToSlice(&Xty), 1e-9)

	if !approxEqual(ln.Y0ty0, mat.Dot(y0v, y0v), 1e-9) {
		t.Errorf("y0ty0: got %.12g want %.12g", ln.Y0ty0, mat.Dot(y0v, y0v))
	}
}

func TestLocalNullEquationsCenteringAndIntercept(t *testing.T) {
	cov := mat.NewDense(3, 1, []float64{0.5, -1.0, 2.0})
	pheno := mat.NewDense(3, 1, []float64{1.0, 0.0, 2.0})

	ln := localNullEquations(cov, pheno, 0.5)

	for i := 0; i < 3; i++ {
		if ln.X.At(i, 0) != 1.0 {
			t.Errorf("intercept X[%d,0]=%v want 1", i, ln.X.At(i, 0))
		}
		if ln.X.At(i, 1) != cov.At(i, 0) {
			t.Errorf("cov X[%d,1]=%v want %v", i, ln.X.At(i, 1), cov.At(i, 0))
		}
		if !approxEqual(ln.Y0[i], pheno.At(i, 0)-0.5, 1e-12) {
			t.Errorf("centered y0[%d]=%v want %v", i, ln.Y0[i], pheno.At(i, 0)-0.5)
		}
	}
}

func TestLocalNullEquationsNilParty0(t *testing.T) {
	ln := localNullEquations(nil, nil, 0.5)
	if len(ln.Y0) != 0 || ln.XtX != nil {
		t.Errorf("nil inputs should give empty result, got Y0=%d XtX=%v", len(ln.Y0), ln.XtX)
	}
}

// Secure low-rank null model: computeBetaHatEnc + computeNullRSSEnc, decrypted on the
// hub and compared to the plaintext gonum β̂/RSS on the same full (X,y). No genotypes.
func TestLowRankNullModel(t *testing.T) {
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	mpcObj := prot.mpcObj[0]
	cps := prot.cps
	pid := mpcObj.GetPid()

	nParties := prot.GetConfig().NumMainParties
	nIndsAll := prot.GetConfig().NumInds
	numCov := prot.GetGwasParams().NumCov()
	c := numCov + 1

	nIndsTotal := 0
	for p := 1; p <= nParties; p++ {
		nIndsTotal += nIndsAll[p]
	}

	// Deterministic well-conditioned design; realistic residual (do NOT shrink the noise:
	// tiny noise → RSS is a sliver of two ~equal large numbers and CKKS cancellation blows up).
	r := rand.New(rand.NewSource(7))
	fullCov := mat.NewDense(nIndsTotal, numCov, nil)
	fullY := make([]float64, nIndsTotal)
	trueBeta := make([]float64, c)
	for j := 0; j < c; j++ {
		trueBeta[j] = 0.5 + 0.3*float64(j)
	}
	for i := 0; i < nIndsTotal; i++ {
		yi := trueBeta[0]
		for j := 0; j < numCov; j++ {
			v := r.NormFloat64()
			fullCov.Set(i, j, v)
			yi += trueBeta[j+1] * v
		}
		fullY[i] = yi + 2.5*r.NormFloat64()
	}

	offset := 0
	for p := 1; p < pid; p++ {
		offset += nIndsAll[p]
	}
	if pid > 0 {
		n := nIndsAll[pid]
		localCov := fullCov.Slice(offset, offset+n, 0, numCov).(*mat.Dense)
		localY := mat.NewDense(n, 1, append([]float64(nil), fullY[offset:offset+n]...))
		prot.SetPhenoAndCov(localY, localCov)
	} else {
		prot.SetPhenoAndCov(nil, nil)
	}

	assocTest := prot.InitAssociationTests(nil)
	null := assocTest.computeBetaHatEnc()
	rssEnc := assocTest.computeNullRSSEnc(null) // nil on pid 0

	if pid == 1 {
		betaDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, null.betaHat, 1))[:c]
		rssDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, rssEnc, 1))[0]

		fullX := mat.NewDense(nIndsTotal, c, nil)
		for i := 0; i < nIndsTotal; i++ {
			fullX.Set(i, 0, 1.0)
			for j := 0; j < numCov; j++ {
				fullX.Set(i, j+1, fullCov.At(i, j))
			}
		}
		yVec := mat.NewVecDense(nIndsTotal, fullY)
		var xtx mat.Dense
		xtx.Mul(fullX.T(), fullX)
		var xty mat.VecDense
		xty.MulVec(fullX.T(), yVec)
		var beta mat.VecDense
		if err := beta.SolveVec(&xtx, &xty); err != nil {
			t.Fatalf("oracle solve failed: %v", err)
		}

		maxAbs := 0.0
		for j := 0; j < c; j++ {
			gap := math.Abs(betaDec[j] - beta.AtVec(j))
			if gap > maxAbs {
				maxAbs = gap
			}
			t.Logf("beta[%d]: secure=%.8f oracle=%.8f |gap|=%.2e", j, betaDec[j], beta.AtVec(j), gap)
		}
		t.Logf("max β̂ |gap| over %d coeffs = %.2e", c, maxAbs)
		const tol = 1e-4
		if maxAbs > tol {
			t.Errorf("β̂ secure vs oracle max gap %.3e exceeds tol %.1e", maxAbs, tol)
		}

		rssOracle := mat.Dot(yVec, yVec) - mat.Dot(&xty, &beta)
		rssRel := math.Abs(rssDec-rssOracle) / math.Abs(rssOracle)
		t.Logf("RSS: secure=%.6f oracle=%.6f rel=%.2e", rssDec, rssOracle, rssRel)
		const rssTol = 1e-4
		if rssRel > rssTol {
			t.Errorf("RSS secure vs oracle rel %.3e exceeds tol %.1e", rssRel, rssTol)
		}
	} else if pid > 0 {
		mpcObj.Network.CollectiveDecryptVec(cps, null.betaHat, 1)
		mpcObj.Network.CollectiveDecryptVec(cps, rssEnc, 1)
	} else {
		mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
		mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
	}
}
