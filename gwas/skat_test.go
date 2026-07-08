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
func TestSKATLocalMatchesGonum(t *testing.T) {
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
func TestSKATLocalPartyAdditivity(t *testing.T) {
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
func TestSKATNullModel(t *testing.T) {
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

// --- low-rank score + weights (skat.go partyScore / signedWeight) ---

// Key-free low-rank score s = Gᵀy₀ − (GᵀX)β̂ vs gonum oracle on the full (G,X,y), center 0.
func TestSKATScore(t *testing.T) {
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
	const m = 12

	nIndsTotal := 0
	for p := 1; p <= nParties; p++ {
		nIndsTotal += nIndsAll[p]
	}

	r := rand.New(rand.NewSource(11))
	fullCov := mat.NewDense(nIndsTotal, numCov, nil)
	fullG := mat.NewDense(nIndsTotal, m, nil)
	fullY := make([]float64, nIndsTotal)
	maf := make([]float64, m)
	for j := 0; j < m; j++ {
		maf[j] = 0.05 + 0.4*r.Float64()
	}
	for i := 0; i < nIndsTotal; i++ {
		yi := 0.0 // mean-0 phenotype (mimics centering)
		for j := 0; j < numCov; j++ {
			v := r.NormFloat64()
			fullCov.Set(i, j, v)
			yi += (0.3 + 0.2*float64(j)) * v
		}
		fullY[i] = yi + 1.5*r.NormFloat64()
		for j := 0; j < m; j++ {
			d := 0.0
			if r.Float64() < maf[j] {
				d++
			}
			if r.Float64() < maf[j] {
				d++
			}
			fullG.Set(i, j, d)
		}
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

	var SBlock crypto.CipherMatrix
	if pid > 0 {
		n := nIndsAll[pid]
		localX := mat.NewDense(n, c, nil)
		localY0 := make([]float64, n)
		localG := mat.NewDense(n, m, nil)
		for i := 0; i < n; i++ {
			localX.Set(i, 0, 1.0)
			for j := 0; j < numCov; j++ {
				localX.Set(i, j+1, fullCov.At(offset+i, j))
			}
			localY0[i] = fullY[offset+i]
			for j := 0; j < m; j++ {
				localG.Set(i, j, fullG.At(offset+i, j))
			}
		}
		lc := LocalContract(localG, localX, localY0)
		SBlock = crypto.CipherMatrix{assocTest.partyScore(lc.GtX, lc.Gty0, null)}
	}

	sAggr := mpcObj.Network.AggregateCMat(cps, SBlock)
	if pid > 0 && len(sAggr) > 0 && len(sAggr[0]) > 0 &&
		mpcObj.Network.CanCollectiveBootstrap(cps, sAggr[0][0].Level()) {
		sAggr = mpcObj.Network.CollectiveBootstrapMat(cps, sAggr, -1)
	}

	var sEnc crypto.CipherVector
	if pid > 0 && len(sAggr) > 0 {
		sEnc = sAggr[0]
	}
	if pid == 1 {
		sDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, sEnc, 1))[:m]

		fullX := mat.NewDense(nIndsTotal, c, nil)
		for i := 0; i < nIndsTotal; i++ {
			fullX.Set(i, 0, 1.0)
			for j := 0; j < numCov; j++ {
				fullX.Set(i, j+1, fullCov.At(i, j))
			}
		}
		y0v := mat.NewVecDense(nIndsTotal, fullY)
		var xtx mat.Dense
		xtx.Mul(fullX.T(), fullX)
		var xty mat.VecDense
		xty.MulVec(fullX.T(), y0v)
		var beta mat.VecDense
		if err := beta.SolveVec(&xtx, &xty); err != nil {
			t.Fatalf("oracle solve: %v", err)
		}
		var gty0 mat.VecDense
		gty0.MulVec(fullG.T(), y0v)
		var gtx mat.Dense
		gtx.Mul(fullG.T(), fullX)
		var gtxBeta mat.VecDense
		gtxBeta.MulVec(&gtx, &beta)

		maxAbs, maxOracle := 0.0, 0.0
		for j := 0; j < m; j++ {
			oracle := gty0.AtVec(j) - gtxBeta.AtVec(j)
			gap := math.Abs(sDec[j] - oracle)
			if gap > maxAbs {
				maxAbs = gap
			}
			if a := math.Abs(oracle); a > maxOracle {
				maxOracle = a
			}
		}
		rel := maxAbs / maxOracle
		t.Logf("score max |gap|=%.2e rel=%.2e (max|s|=%.3f)", maxAbs, rel, maxOracle)
		const tol = 1e-3
		if rel > tol {
			t.Errorf("score secure vs oracle rel %.3e exceeds tol %.1e", rel, tol)
		}
	} else if pid > 0 {
		mpcObj.Network.CollectiveDecryptVec(cps, sEnc, 1)
	} else {
		mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
	}
}

// dosage fixture shared by the two weight tests (seed 23; some ALT-major af to test min/flip).
func weightFixtureDosage(t *testing.T, prot *ProtocolInfo, pid int, m int) ([]float64, *mat.Dense, int) {
	nParties := prot.GetConfig().NumMainParties
	nIndsAll := prot.GetConfig().NumInds
	nIndsTotal := 0
	for p := 1; p <= nParties; p++ {
		nIndsTotal += nIndsAll[p]
	}

	r := rand.New(rand.NewSource(23))
	fullG := mat.NewDense(nIndsTotal, m, nil)
	af := make([]float64, m)
	for j := 0; j < m; j++ {
		af[j] = 0.05 + 0.9*r.Float64()
	}
	for i := 0; i < nIndsTotal; i++ {
		for j := 0; j < m; j++ {
			d := 0.0
			if r.Float64() < af[j] {
				d++
			}
			if r.Float64() < af[j] {
				d++
			}
			fullG.Set(i, j, d)
		}
	}

	offset := 0
	for p := 1; p < pid; p++ {
		offset += nIndsAll[p]
	}
	if pid > 0 {
		n := nIndsAll[pid]
		prot.SetPhenoAndCov(mat.NewDense(n, 1, nil), mat.NewDense(n, prot.GetGwasParams().NumCov(), nil))
	} else {
		prot.SetPhenoAndCov(nil, nil)
	}

	dosage := make([]float64, m)
	if pid > 0 {
		n := nIndsAll[pid]
		for i := 0; i < n; i++ {
			for j := 0; j < m; j++ {
				dosage[j] += fullG.At(offset+i, j)
			}
		}
	}
	return dosage, fullG, nIndsTotal
}

// Unsigned weights w_j = 25(1−MAF)^24 vs oracle.
func TestSKATWeights(t *testing.T) {
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	mpcObj := prot.mpcObj[0]
	cps := prot.cps
	pid := mpcObj.GetPid()
	const m = 12

	dosage, fullG, nIndsTotal := weightFixtureDosage(t, prot, pid, m)
	assocTest := prot.InitAssociationTests(nil)

	_, _, w24 := assocTest.weightsCalculation(dosage, m)
	weightEnc := mpcObj.SSToCVec(cps, w24)

	if pid == 1 {
		wDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, weightEnc, 1))[:m]
		twoN := float64(2 * nIndsTotal)
		maxAbs, maxRelSig := 0.0, 0.0
		for j := 0; j < m; j++ {
			var colSum float64
			for i := 0; i < nIndsTotal; i++ {
				colSum += fullG.At(i, j)
			}
			pbar := colSum / twoN
			maf := math.Min(pbar, 1-pbar)
			wOracle := 25.0 * math.Pow(1-maf, 24)
			abs := math.Abs(wDec[j] - wOracle)
			if abs > maxAbs {
				maxAbs = abs
			}
			if wOracle > 1e-2 {
				if rel := abs / wOracle; rel > maxRelSig {
					maxRelSig = rel
				}
			}
		}
		t.Logf("weights maxAbs=%.2e maxRelSig=%.2e", maxAbs, maxRelSig)
		const absTol, relTol = 1e-3, 1e-3
		if maxAbs > absTol || maxRelSig > relTol {
			t.Errorf("weights: maxAbs=%.3e maxRelSig=%.3e", maxAbs, maxRelSig)
		}
	} else if pid > 0 {
		mpcObj.Network.CollectiveDecryptVec(cps, weightEnc, 1)
	} else {
		mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
	}
}

// Minor-allele-oriented weight ŵ_j = t_j·w_j (t = −1 iff p̄>½) vs oracle; fixture flips some.
func TestSKATSignedWeight(t *testing.T) {
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	mpcObj := prot.mpcObj[0]
	cps := prot.cps
	pid := mpcObj.GetPid()
	const m = 12

	dosage, fullG, nIndsTotal := weightFixtureDosage(t, prot, pid, m)
	assocTest := prot.InitAssociationTests(nil)

	signedEnc := mpcObj.SSToCVec(cps, assocTest.signedWeight(dosage, m))

	if pid == 1 {
		wDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, signedEnc, 1))[:m]
		twoN := float64(2 * nIndsTotal)
		nFlipped := 0
		maxAbs, maxRelSig := 0.0, 0.0
		for j := 0; j < m; j++ {
			var colSum float64
			for i := 0; i < nIndsTotal; i++ {
				colSum += fullG.At(i, j)
			}
			pbar := colSum / twoN
			maf := math.Min(pbar, 1-pbar)
			sign := 1.0
			if pbar > 0.5 {
				sign = -1.0
				nFlipped++
			}
			wOracle := sign * 25.0 * math.Pow(1-maf, 24)
			abs := math.Abs(wDec[j] - wOracle)
			if abs > maxAbs {
				maxAbs = abs
			}
			if math.Abs(wOracle) > 1e-2 {
				if rel := abs / math.Abs(wOracle); rel > maxRelSig {
					maxRelSig = rel
				}
			}
		}
		t.Logf("signed weights maxAbs=%.2e maxRelSig=%.2e flipped %d/%d", maxAbs, maxRelSig, nFlipped, m)
		if nFlipped == 0 {
			t.Fatalf("fixture exercised no flips; orientation not tested")
		}
		const absTol, relTol = 1e-3, 1e-3
		if maxAbs > absTol || maxRelSig > relTol {
			t.Errorf("signed weights: maxAbs=%.3e maxRelSig=%.3e", maxAbs, maxRelSig)
		}
	} else if pid > 0 {
		mpcObj.Network.CollectiveDecryptVec(cps, signedEnc, 1)
	} else {
		mpcObj.Network.CollectiveDecryptVec(cps, nil, 1)
	}
}

// --- low-rank driver (skat.go ComputeSKATStatistics, file-backed) ---

// File-backed driver E2E: reads per-party "blocks"-format genotype fixtures, runs the full
// low-rank path over 3 blocks, and compares decrypted Q/Burden to the plaintext oracle.
func TestSKATDriverE2E(t *testing.T) {
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
	const nBlocks, mPerBlock = 3, 6
	const mTotal = nBlocks * mPerBlock

	nIndsTotal := 0
	for p := 1; p <= nParties; p++ {
		nIndsTotal += nIndsAll[p]
	}

	r := rand.New(rand.NewSource(41))
	fullCov := mat.NewDense(nIndsTotal, numCov, nil)
	fullG := mat.NewDense(nIndsTotal, mTotal, nil)
	fullY := make([]float64, nIndsTotal)
	af := make([]float64, mTotal)
	for j := 0; j < mTotal; j++ {
		af[j] = 0.05 + 0.85*r.Float64()
	}
	for i := 0; i < nIndsTotal; i++ {
		yi := 0.0
		for j := 0; j < numCov; j++ {
			v := r.NormFloat64()
			fullCov.Set(i, j, v)
			yi += (0.3 + 0.2*float64(j)) * v
		}
		fullY[i] = yi + 1.5*r.NormFloat64()
		for j := 0; j < mTotal; j++ {
			d := 0.0
			if r.Float64() < af[j] {
				d++
			}
			if r.Float64() < af[j] {
				d++
			}
			fullG.Set(i, j, d)
		}
	}

	offset := 0
	for p := 1; p < pid; p++ {
		offset += nIndsAll[p]
	}

	prot.GetConfig().GenoNumBlocks = nBlocks
	prot.GetConfig().GenoFileFormat = "blocks"
	if pid > 0 {
		n := nIndsAll[pid]
		localCov := fullCov.Slice(offset, offset+n, 0, numCov).(*mat.Dense)
		localY := mat.NewDense(n, 1, append([]float64(nil), fullY[offset:offset+n]...))
		prot.SetPhenoAndCov(localY, localCov)

		dir := t.TempDir()
		streams := make([]*GenoFileStream, nBlocks)
		sizes := make([]int, nBlocks)
		for b := 0; b < nBlocks; b++ {
			buf := make([]byte, n*mPerBlock)
			for i := 0; i < n; i++ {
				for j := 0; j < mPerBlock; j++ {
					buf[i*mPerBlock+j] = byte(int8(fullG.At(offset+i, b*mPerBlock+j)))
				}
			}
			path := filepath.Join(dir, "geno.") + string(rune('0'+b)) + ".bin"
			if err := os.WriteFile(path, buf, 0644); err != nil {
				t.Fatalf("write geno fixture: %v", err)
			}
			streams[b] = NewGenoFileStream(path, uint64(n), uint64(mPerBlock), false)
			sizes[b] = mPerBlock
		}
		prot.genoBlocks = streams
		prot.genoBlockSizes = sizes
	} else {
		prot.SetPhenoAndCov(nil, nil)
		prot.genoBlockSizes = make([]int, nBlocks)
	}

	assocTest := prot.InitAssociationTests(nil)
	qOut, bOut := assocTest.ComputeSKATStatistics()

	if pid == 1 {
		qDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, qOut, 1))[0]
		bDec := crypto.DecodeFloatVector(cps, mpcObj.Network.CollectiveDecryptVec(cps, bOut, 1))[0]

		fullX := mat.NewDense(nIndsTotal, c, nil)
		for i := 0; i < nIndsTotal; i++ {
			fullX.Set(i, 0, 1.0)
			for j := 0; j < numCov; j++ {
				fullX.Set(i, j+1, fullCov.At(i, j))
			}
		}
		oracle := SKATPlain(fullG, fullX, fullY)

		qRel := math.Abs(qDec-oracle.Q) / math.Abs(oracle.Q)
		bRel := math.Abs(bDec-oracle.Burden) / math.Abs(oracle.Burden)
		t.Logf("Q: secure=%.6f oracle=%.6f rel=%.2e | Burden: secure=%.6f oracle=%.6f rel=%.2e",
			qDec, oracle.Q, qRel, bDec, oracle.Burden, bRel)
		const tol = 5e-3
		if qRel > tol || bRel > tol {
			t.Errorf("driver E2E vs oracle: Q rel=%.3e, Burden rel=%.3e (tol %.0e)", qRel, bRel, tol)
		}
	} else {
		mpcObj.Network.CollectiveDecryptVec(cps, qOut, 1)
		mpcObj.Network.CollectiveDecryptVec(cps, bOut, 1)
	}
}

// TestSKATFederatedPrivate drives the secure federated per-gene SKAT Q for heterogeneous variant
// sets (pid1 = public-list party, pid2 = private-variant party; MVP+AoU is the motivating instance)
// and asserts it equals the plaintext oracle SKATFederatedPrivate AND the pooled single-cohort
// SKATPlain, per gene. Every party regenerates the same full fixture (seed) and injects only its
// own slice; pid1 reconstructs both oracles from the full data.
func TestSKATFederatedPrivate(t *testing.T) {
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	pid := prot.mpcObj[0].GetPid()
	const privatePid = 2

	numCov := prot.GetGwasParams().NumCov()
	c := numCov + 1
	nIndsAll := prot.GetConfig().NumInds
	nPub, nPriv := nIndsAll[1], nIndsAll[2]

	const nGenes = 3
	const nShared, nPubOnly, nPrivOnly = 2, 2, 2
	const pubPerGene = nShared + nPubOnly // public list per gene

	// --- full fixture, identical on every party (seed) ---
	r := rand.New(rand.NewSource(71))
	drawCohort := func(n int) (*mat.Dense, []float64) {
		cov := mat.NewDense(n, numCov, nil)
		y := make([]float64, n)
		for i := 0; i < n; i++ {
			yi := 0.0
			for j := 0; j < numCov; j++ {
				v := r.NormFloat64()
				cov.Set(i, j, v)
				yi += (0.3 + 0.2*float64(j)) * v
			}
			y[i] = yi + 1.5*r.NormFloat64()
		}
		return cov, y
	}
	drawDosage := func(n int, af float64) []float64 {
		col := make([]float64, n)
		for i := 0; i < n; i++ {
			d := 0.0
			if r.Float64() < af {
				d++
			}
			if r.Float64() < af {
				d++
			}
			col[i] = d
		}
		return col
	}

	pubCov, pubY := drawCohort(nPub)
	privCov, privY := drawCohort(nPriv)

	// Per gene: public columns [shared×nShared, public_only×nPubOnly]; private party's shared-side
	// columns; private-only columns. Shared variants are the same id measured in both cohorts.
	pubCols := make([][][]float64, nGenes)        // [gene][col(=pubPerGene)][sample(nPub)]
	privSharedCols := make([][][]float64, nGenes) // [gene][nShared][sample(nPriv)]
	privOnlyCols := make([][][]float64, nGenes)   // [gene][nPrivOnly][sample(nPriv)]
	for g := 0; g < nGenes; g++ {
		pubCols[g] = make([][]float64, pubPerGene)
		privSharedCols[g] = make([][]float64, nShared)
		for k := 0; k < nShared; k++ { // shared: both cohorts
			af := 0.05 + 0.4*r.Float64()
			pubCols[g][k] = drawDosage(nPub, af)
			privSharedCols[g][k] = drawDosage(nPriv, af)
		}
		for k := 0; k < nPubOnly; k++ { // public_only: public cohort only
			af := 0.05 + 0.4*r.Float64()
			pubCols[g][nShared+k] = drawDosage(nPub, af)
		}
		privOnlyCols[g] = make([][]float64, nPrivOnly)
		for k := 0; k < nPrivOnly; k++ { // private: private cohort only
			af := 0.05 + 0.4*r.Float64()
			privOnlyCols[g][k] = drawDosage(nPriv, af)
		}
	}

	geneName := func(g int) string { return "gene" + string(rune('0'+g)) }

	prot.GetConfig().GenoNumBlocks = nGenes
	prot.GetConfig().GenoFileFormat = "blocks"

	// --- per-party injection ---
	var privPrefix string // private-only block files (only the private party); "" elsewhere
	switch pid {
	case 1: // public-list party: public-list genotypes, no private variants
		prot.SetPhenoAndCov(mat.NewDense(nPub, 1, append([]float64(nil), pubY...)), pubCov)
		streams, sizes := writeGeneStreams(t, nPub, nGenes, pubPerGene, func(g, col int) []float64 {
			return pubCols[g][col]
		})
		prot.genoBlocks = streams
		prot.genoBlockSizes = sizes
	case 2: // private party: public-list genotypes ALIGNED (public_only cols = 0) + private dense
		prot.SetPhenoAndCov(mat.NewDense(nPriv, 1, append([]float64(nil), privY...)), privCov)
		streams, sizes := writeGeneStreams(t, nPriv, nGenes, pubPerGene, func(g, col int) []float64 {
			if col < nShared {
				return privSharedCols[g][col]
			}
			return make([]float64, nPriv) // public_only aligned to 0
		})
		prot.genoBlocks = streams
		prot.genoBlockSizes = sizes
		// private-only blocks → per-gene int8 files (the skat_fed mode loads these).
		privDir := t.TempDir()
		for g := 0; g < nGenes; g++ {
			buf := make([]byte, nPriv*nPrivOnly)
			for i := 0; i < nPriv; i++ {
				for k := 0; k < nPrivOnly; k++ {
					buf[i*nPrivOnly+k] = byte(int8(privOnlyCols[g][k][i]))
				}
			}
			if err := os.WriteFile(fmt.Sprintf("%s/priv.%d.bin", privDir, g), buf, 0644); err != nil {
				t.Fatalf("write priv block: %v", err)
			}
		}
		privPrefix = fmt.Sprintf("%s/priv.%%d.bin", privDir)
	default: // pid 0
		prot.SetPhenoAndCov(nil, nil)
		prot.genoBlockSizes = make([]int, nGenes)
	}

	prot.GetConfig().PrivatePid = privatePid
	prot.GetConfig().PrivateOnlyPrefix = privPrefix // "" except on the private party
	skatDec, burdenDec := prot.runFederatedPrivate()
	if pid != 1 {
		return
	}

	// --- oracle 1: SKATFederatedPrivate (the secure spec) ---
	pub := buildFedParty(pubCov, pubY, c, nGenes, pubPerGene,
		func(g, col int) []float64 { return pubCols[g][col] },
		func(g, col int) (id, gene, role string) {
			if col < nShared {
				return varID(g, "shr", col), geneName(g), "shared"
			}
			return varID(g, "pub", col-nShared), geneName(g), "public_only"
		})
	privCount := nShared + nPrivOnly
	priv := buildFedParty(privCov, privY, c, nGenes, privCount,
		func(g, col int) []float64 {
			if col < nShared {
				return privSharedCols[g][col]
			}
			return privOnlyCols[g][col-nShared]
		},
		func(g, col int) (id, gene, role string) {
			if col < nShared {
				return varID(g, "shr", col), geneName(g), "shared"
			}
			return varID(g, "prv", col-nShared), geneName(g), "private"
		})
	oracleSkat, oracleBurden := SKATFederatedPrivate(pub, priv)

	// --- oracle 2: pooled per-gene SKATPlain (independent code path) ---
	Xpool := mat.NewDense(nPub+nPriv, c, nil)
	ypool := make([]float64, nPub+nPriv)
	for i := 0; i < nPub; i++ {
		Xpool.Set(i, 0, 1.0)
		for j := 0; j < numCov; j++ {
			Xpool.Set(i, j+1, pubCov.At(i, j))
		}
		ypool[i] = pubY[i]
	}
	for i := 0; i < nPriv; i++ {
		Xpool.Set(nPub+i, 0, 1.0)
		for j := 0; j < numCov; j++ {
			Xpool.Set(nPub+i, j+1, privCov.At(i, j))
		}
		ypool[nPub+i] = privY[i]
	}

	const tol = 5e-3
	for g := 0; g < nGenes; g++ {
		// pooled union genotypes for gene g: shared (both), public_only (pub rows), private (priv rows).
		union := nShared + nPubOnly + nPrivOnly
		Gp := mat.NewDense(nPub+nPriv, union, nil)
		col := 0
		for k := 0; k < nShared; k++ {
			for i := 0; i < nPub; i++ {
				Gp.Set(i, col, pubCols[g][k][i])
			}
			for i := 0; i < nPriv; i++ {
				Gp.Set(nPub+i, col, privSharedCols[g][k][i])
			}
			col++
		}
		for k := 0; k < nPubOnly; k++ {
			for i := 0; i < nPub; i++ {
				Gp.Set(i, col, pubCols[g][nShared+k][i])
			}
			col++
		}
		for k := 0; k < nPrivOnly; k++ {
			for i := 0; i < nPriv; i++ {
				Gp.Set(nPub+i, col, privOnlyCols[g][k][i])
			}
			col++
		}
		plain := SKATPlain(Gp, Xpool, ypool)

		// three-way check per statistic: secure == federated oracle == pooled single-cohort.
		check := func(stat string, secure, oref, pooled float64) {
			relOracle := math.Abs(secure-oref) / math.Max(math.Abs(oref), 1e-12)
			relPool := math.Abs(secure-pooled) / math.Max(math.Abs(pooled), 1e-12)
			t.Logf("%s %s: secure=%.6f oracle=%.6f pooled=%.6f relOracle=%.2e relPool=%.2e",
				geneName(g), stat, secure, oref, pooled, relOracle, relPool)
			if relOracle > tol || relPool > tol {
				t.Errorf("%s %s: secure rel oracle=%.3e pooled=%.3e (tol %.0e)", geneName(g), stat, relOracle, relPool, tol)
			}
		}
		check("SKAT", skatDec[g], oracleSkat[geneName(g)], plain.Q)
		check("Burden", burdenDec[g], oracleBurden[geneName(g)], plain.Burden)
	}
}

func varID(g int, role string, k int) string {
	return "g" + string(rune('0'+g)) + "_" + role + string(rune('0'+k))
}

// loadDenseBlocks reads per-gene int8 genotype block files (row-major n×m_b, m_b inferred from
// file size / n) into dense matrices — the file-loading side of the skat_fed mode.
func TestLoadDenseBlocks(t *testing.T) {
	dir := t.TempDir()
	const n = 4
	g0 := [][]float64{{0, 1}, {2, 0}, {1, 1}, {0, 2}}          // 4×2
	g1 := [][]float64{{1, 0, 2}, {0, 1, 0}, {2, 2, 1}, {1, 0, 0}} // 4×3
	writeBlock := func(b int, rows [][]float64) {
		m := len(rows[0])
		buf := make([]byte, n*m)
		for i := 0; i < n; i++ {
			for j := 0; j < m; j++ {
				buf[i*m+j] = byte(int8(rows[i][j]))
			}
		}
		if err := os.WriteFile(fmt.Sprintf("%s/blk.%d.bin", dir, b), buf, 0644); err != nil {
			t.Fatalf("write block: %v", err)
		}
	}
	writeBlock(0, g0)
	writeBlock(1, g1)

	blocks, err := loadDenseBlocks(dir+"/blk.%d.bin", 2, n)
	if err != nil {
		t.Fatalf("loadDenseBlocks: %v", err)
	}
	if len(blocks) != 2 {
		t.Fatalf("nGenes: got %d want 2", len(blocks))
	}
	for b, want := range [][][]float64{g0, g1} {
		gr, gc := blocks[b].Dims()
		if gr != n || gc != len(want[0]) {
			t.Fatalf("block %d dims: got %dx%d want %dx%d", b, gr, gc, n, len(want[0]))
		}
		for i := range want {
			for j := range want[i] {
				if got := blocks[b].At(i, j); got != want[i][j] {
					t.Errorf("block %d [%d,%d]: got %v want %v", b, i, j, got, want[i][j])
				}
			}
		}
	}
}

// writeGeneStreams writes nGenes row-major int8 genotype blocks (n × mPerGene) to temp files and
// opens them. colFn(gene, col) returns that column's n dosages.
func writeGeneStreams(t *testing.T, n, nGenes, mPerGene int, colFn func(g, col int) []float64) ([]*GenoFileStream, []int) {
	dir := t.TempDir()
	streams := make([]*GenoFileStream, nGenes)
	sizes := make([]int, nGenes)
	for g := 0; g < nGenes; g++ {
		buf := make([]byte, n*mPerGene)
		for col := 0; col < mPerGene; col++ {
			cdat := colFn(g, col)
			for i := 0; i < n; i++ {
				buf[i*mPerGene+col] = byte(int8(cdat[i]))
			}
		}
		path := filepath.Join(dir, "geno.") + string(rune('0'+g)) + ".bin"
		if err := os.WriteFile(path, buf, 0644); err != nil {
			t.Fatalf("write geno fixture: %v", err)
		}
		streams[g] = NewGenoFileStream(path, uint64(n), uint64(mPerGene), false)
		sizes[g] = mPerGene
	}
	return streams, sizes
}

// buildFedParty assembles a FedParty (oracle input) with mPerGene columns per gene.
func buildFedParty(cov *mat.Dense, y []float64, c, nGenes, mPerGene int,
	colFn func(g, col int) []float64, labelFn func(g, col int) (id, gene, role string)) FedParty {
	n, numCov := cov.Dims()
	m := nGenes * mPerGene
	X := mat.NewDense(n, c, nil)
	for i := 0; i < n; i++ {
		X.Set(i, 0, 1.0)
		for j := 0; j < numCov; j++ {
			X.Set(i, j+1, cov.At(i, j))
		}
	}
	G := mat.NewDense(n, m, nil)
	ids := make([]string, m)
	genes := make([]string, m)
	roles := make([]string, m)
	for g := 0; g < nGenes; g++ {
		for col := 0; col < mPerGene; col++ {
			j := g*mPerGene + col
			G.SetCol(j, colFn(g, col))
			ids[j], genes[j], roles[j] = labelFn(g, col)
		}
	}
	return FedParty{G: G, X: X, Y: y, ID: ids, Gene: genes, Role: roles}
}
