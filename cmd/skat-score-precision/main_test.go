package main

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"testing"

	"github.com/tuneinsight/lattigo/v6/multiparty/mpckks"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func TestChooseGenesEvenly(t *testing.T) {
	if got, want := chooseGenes(10, 4), []int{0, 3, 6, 9}; !reflect.DeepEqual(got, want) {
		t.Fatalf("chooseGenes: got %v want %v", got, want)
	}
	if got := chooseGenes(3, 0); !reflect.DeepEqual(got, []int{0, 1, 2}) {
		t.Fatalf("choose all: %v", got)
	}
}

func TestSyntheticFedPrepInputsEndToEnd(t *testing.T) {
	root := t.TempDir()
	for _, d := range []string{"A", "B", "config"} {
		if err := os.Mkdir(filepath.Join(root, d), 0o700); err != nil {
			t.Fatal(err)
		}
	}
	m := manifest{NGenes: 1, GeneNames: []string{"SYNTH"}, PublicM: []int{2}, PrivateM: []int{1}}
	b, _ := json.Marshal(m)
	writeTestFile(t, filepath.Join(root, "manifest.json"), b)
	writeTestFile(t, filepath.Join(root, "config", "configGlobal.toml"), []byte(
		"ckks_params = \"PN14QP438\"\nmpc_frac_bits = 30\nmpc_field_size = 128\nnum_covs = 1\nbinary_pheno = false\n"))
	writeTestFile(t, filepath.Join(root, "A", "pheno.txt"), []byte("1\n2\n3\n"))
	writeTestFile(t, filepath.Join(root, "B", "pheno.txt"), []byte("2\n1\n4\n"))
	writeTestFile(t, filepath.Join(root, "A", "cov.txt"), []byte("-1\n0\n1\n"))
	writeTestFile(t, filepath.Join(root, "B", "cov.txt"), []byte("1\n0\n-1\n"))
	writeTestFile(t, filepath.Join(root, "A", "geno.0.bin"), []byte{0, 1, 1, 0, 2, 1})
	writeTestFile(t, filepath.Join(root, "B", "geno.0.bin"), []byte{1, 0, 0, 2, 1, 1})
	writeTestFile(t, filepath.Join(root, "B", "priv.0.bin"), []byte{0, 1, 2})

	loaded, err := loadManifest(root)
	if err != nil || loaded.NGenes != 1 {
		t.Fatalf("manifest: %+v %v", loaded, err)
	}
	cfg, err := loadGlobalConfig(root)
	if err != nil {
		t.Fatal(err)
	}
	a, err := loadCohort(root, "A", cfg.BinaryPheno, cfg.NumCovs)
	if err != nil {
		t.Fatal(err)
	}
	c, err := loadCohort(root, "B", cfg.BinaryPheno, cfg.NumCovs)
	if err != nil {
		t.Fatal(err)
	}
	beta, err := fitNull(&a, &c, cfg.MPCFracBits)
	if err != nil {
		t.Fatal(err)
	}
	ga, err := contractBlock(filepath.Join(root, "A", "geno.0.bin"), &a, 2)
	if err != nil {
		t.Fatal(err)
	}
	gb, err := contractBlock(filepath.Join(root, "B", "geno.0.bin"), &c, 2)
	if err != nil {
		t.Fatal(err)
	}
	records, _ := recordsForGene(0, ga, gb, beta, a.N+c.N, cfg.MPCFracBits)
	metrics := make([]categoryMetrics, 1)
	engine, err := newHEEngine(cfg.CKKSParams, beta, uint(cfg.MPCFieldSize), 1)
	if err != nil {
		t.Fatal(err)
	}
	chunk := newScoreChunker(engine, metrics, true)
	for _, r := range records {
		metrics[0].addPlain(r)
		if err := chunk.Add(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := chunk.Flush(); err != nil {
		t.Fatal(err)
	}
	if !metrics[0].HasHE || math.IsNaN(metrics[0].heQRel()) {
		t.Fatalf("missing HE metric: %+v", metrics[0])
	}
}

func writeTestFile(t *testing.T, path string, data []byte) {
	t.Helper()
	if err := os.WriteFile(path, data, 0o600); err != nil {
		t.Fatal(err)
	}
}

func TestContractBlockUsesLocalMinorOrientation(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "geno.bin")
	// Rows: [2 0], [2 1], [-1 2]. Missing maps to zero before orientation.
	// Column 0 has dosage 4>n=3 and is recoded to [0,0,2]; column 1 is a
	// tie and stays [0,1,2].
	if err := os.WriteFile(path, []byte{2, 0, 2, 1, 255, 2}, 0o600); err != nil {
		t.Fatal(err)
	}
	co := cohort{
		N: 3, C: 1,
		X:        []float64{1, 1, 1},
		Y:        []float64{1, 2, 3},
		Residual: []float64{1, 2, 3},
	}
	got, err := contractBlock(path, &co, 2)
	if err != nil {
		t.Fatal(err)
	}
	wantDose := []float64{2, 3}
	wantGty := []float64{6, 8}
	for j := range wantDose {
		if got.Dosage[j] != wantDose[j] || got.Gty[j] != wantGty[j] || got.Stable[j] != wantGty[j] {
			t.Fatalf("column %d: dosage/gty/stable = %.1f/%.1f/%.1f", j, got.Dosage[j], got.Gty[j], got.Stable[j])
		}
		if got.Gtx[0][j] != wantDose[j] {
			t.Fatalf("column %d intercept contraction %.1f want %.1f", j, got.Gtx[0][j], wantDose[j])
		}
	}
}

func TestHEMetricFailsClosedOnNonFinite(t *testing.T) {
	var m categoryMetrics
	m.addPlain(scoreRecord{Weight: 1, Reference: 2, Cancel64: 2, SSLike: 2, TermBound: 4})
	m.addHE(1, 2, math.NaN())
	if m.HEInvalid != 1 || !math.IsInf(m.heQRel(), 1) || !math.IsInf(m.heLNorm(), 1) {
		t.Fatalf("non-finite HE result was not fail-closed: %+v", m)
	}
}

func TestQAuditReportsZeroAndInvalidGenes(t *testing.T) {
	zero := categoryMetrics{}
	bad := categoryMetrics{Variants: 1, RefQ: 1, HasHE: true, HEInvalid: 1}
	rows := []geneResult{{Public: zero}, {Public: bad}}
	a := auditQ(rows, func(g geneResult) categoryMetrics { return g.Public }, func(m categoryMetrics) float64 { return m.HEQ }, true)
	if a.Total != 2 || a.Degenerate != 1 || a.Invalid != 1 || len(a.Relative) != 0 {
		t.Fatalf("unexpected audit: %+v", a)
	}
}

func TestQAuditKeepsNoisyZeroQGene(t *testing.T) {
	var m categoryMetrics
	m.addPlain(scoreRecord{Weight: 1, Reference: 0, Cancel64: 0, SSLike: 0, TermBound: 1})
	m.addHE(1, 0, 1e-4)
	a := auditQ([]geneResult{{Public: m}}, func(g geneResult) categoryMetrics { return g.Public }, func(m categoryMetrics) float64 { return m.HEQ }, true)
	if a.Degenerate != 1 || len(a.DegenerateAbs) != 1 || math.Abs(a.DegenerateAbs[0]-1e-8) > 1e-20 {
		t.Fatalf("noisy zero-Q gene was lost: %+v", a)
	}
	if !isFinite(m.heLNorm()) || m.heLNorm() <= 0 {
		t.Fatalf("zero-L normalized error should remain explicit: %g", m.heLNorm())
	}
}

func TestContractBlockCachedSumsMatchLiteralOrientation(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "geno.bin")
	raw := []int8{2, 0, 2, 1, -1, 2, 1, -1}
	encoded := make([]byte, len(raw))
	for i, v := range raw {
		encoded[i] = byte(v)
	}
	if err := os.WriteFile(path, encoded, 0o600); err != nil {
		t.Fatal(err)
	}
	co := cohort{
		N: 4,
		C: 2,
		X: []float64{
			1, 1e6,
			1, -2e6,
			1, 3e6,
			1, -4e6,
		},
		Y:        []float64{1e8 + 1, -1e8 + 2, 5e7 + 3, -5e7 + 4},
		Residual: []float64{0.25, -0.5, 0.75, -1},
	}
	co.SumX = []float64{4, -2e6}
	for _, v := range co.Y {
		co.SumY += v
	}
	for _, v := range co.Residual {
		co.SumResidual += v
	}
	got, err := contractBlock(path, &co, 2)
	if err != nil {
		t.Fatal(err)
	}
	want := literalContract(raw, &co, 2)
	compareContraction(t, got, want)
}

func literalContract(raw []int8, co *cohort, m int) contraction {
	out := emptyContraction(co.C, m)
	dosage := make([]float64, m)
	for i := 0; i < co.N; i++ {
		for j := 0; j < m; j++ {
			if raw[i*m+j] > 0 {
				dosage[j] += float64(raw[i*m+j])
			}
		}
	}
	for i := 0; i < co.N; i++ {
		for j := 0; j < m; j++ {
			v := 0.0
			if raw[i*m+j] > 0 {
				v = float64(raw[i*m+j])
			}
			if dosage[j] > float64(co.N) {
				v = 2 - v
			}
			out.Dosage[j] += v
			out.Gty[j] += v * co.Y[i]
			out.Stable[j] += v * co.Residual[i]
			for k := 0; k < co.C; k++ {
				out.Gtx[k][j] += v * co.X[i*co.C+k]
			}
		}
	}
	return out
}

func compareContraction(t *testing.T, got, want contraction) {
	t.Helper()
	check := func(label string, a, b float64) {
		t.Helper()
		tol := 1e-12 * math.Max(1, math.Abs(b))
		if math.Abs(a-b) > tol {
			t.Fatalf("%s: got %.17g want %.17g (tol %.3g)", label, a, b, tol)
		}
	}
	for j := range want.Gty {
		check("dosage", got.Dosage[j], want.Dosage[j])
		check("gty", got.Gty[j], want.Gty[j])
		check("stable", got.Stable[j], want.Stable[j])
		for k := range want.Gtx {
			check("gtx", got.Gtx[k][j], want.Gtx[k][j])
		}
	}
}

func TestMetricKappaQAndRequiredPrecision(t *testing.T) {
	var m categoryMetrics
	m.addPlain(scoreRecord{Weight: 2, Reference: 3, Cancel64: 3, TermBound: 6})
	// kappa_Q = |w*t|/|w*s| = 2.
	if math.Abs(m.kappaQ()-2) > 1e-12 {
		t.Fatalf("kappaQ %.16g", m.kappaQ())
	}
	want := (math.Sqrt(1.01) - 1) / 2
	if math.Abs(m.requiredEpsilon1PctQ()-want) > 1e-15 {
		t.Fatalf("required epsilon %.16g want %.16g", m.requiredEpsilon1PctQ(), want)
	}
}

func TestPackedCKKSScoreMatchesSameBetaPlaintext(t *testing.T) {
	beta := []float64{0.25, -0.4}
	gty1 := []float64{10, -3, 8, 0.25}
	gty2 := []float64{2, 5, -1, -0.75}
	gtx1 := [][]float64{{4, 2, 7, 1}, {1, -2, 3, 4}}
	gtx2 := [][]float64{{3, 1, 0, 2}, {-1, 5, 2, -3}}
	for _, profile := range []string{"PN14QP438", "PN14QP427S38", "PN14QP431S40", "PN14QP436S45"} {
		t.Run(profile, func(t *testing.T) {
			engine, err := newHEEngine(profile, beta, 128, 1)
			if err != nil {
				t.Fatal(err)
			}
			party := engine.partyScore(gty1, gtx1)
			if got, want := cipherLevel(party), engine.MaxLevel-1; got != want {
				t.Fatalf("party score level=%d want=%d", got, want)
			}
			for j, ct := range party {
				if gap := math.Abs(math.Log2(ct.Scale.Float64()) - float64(engine.LogDefaultScale)); gap > 0.25 {
					t.Fatalf("ciphertext %d log-scale gap %.3f bits", j, gap)
				}
			}
			got := engine.evaluate(gty1, gtx1, gty2, gtx2, true)
			for i := range got {
				want := gty1[i] + gty2[i]
				for k := range beta {
					want -= (gtx1[k][i] + gtx2[k][i]) * beta[k]
				}
				// This is a circuit-wiring regression, not the scientific
				// acceptance threshold. The full Workbench dataset determines
				// whether a profile handles realistic cancellation.
				if gap := math.Abs(got[i] - want); gap > 5e-3 {
					t.Fatalf("slot %d: got %.10f want %.10f gap %.3e", i, got[i], want, gap)
				}
			}
		})
	}
}

func TestN14HighScaleProfileBudgets(t *testing.T) {
	tests := []struct {
		name      string
		logScale  int
		maxLevel  int
		logQ      []int
		logP      []int
		logQP     float64
		depth     int
		parties   int
		security  int
		wantSlots int
		wantLogN  int
	}{
		{name: "PN14QP427S38", logScale: 38, maxLevel: 7, logQ: []int{40, 39, 38, 38, 38, 38, 38, 38}, logP: []int{40, 40, 40}, logQP: 427, depth: 3, parties: 2, security: 128, wantSlots: 8192, wantLogN: 14},
		{name: "PN14QP431S40", logScale: 40, maxLevel: 7, logQ: []int{51, 40, 40, 40, 40, 40, 40, 40}, logP: []int{50, 50}, logQP: 431, depth: 3, parties: 2, security: 128, wantSlots: 8192, wantLogN: 14},
		{name: "PN14QP436S45", logScale: 45, maxLevel: 6, logQ: []int{56, 45, 45, 45, 45, 45, 45}, logP: []int{55, 55}, logQP: 436, depth: 3, parties: 2, security: 128, wantSlots: 8192, wantLogN: 14},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			lit, err := parameterLiteral(tc.name)
			if err != nil {
				t.Fatal(err)
			}
			params, err := ckks.NewParametersFromLiteral(lit)
			if err != nil {
				t.Fatal(err)
			}
			if params.LogN() != tc.wantLogN || params.MaxSlots() != tc.wantSlots {
				t.Fatalf("shape: logN=%d slots=%d", params.LogN(), params.MaxSlots())
			}
			if params.LogDefaultScale() != tc.logScale || params.MaxLevel() != tc.maxLevel {
				t.Fatalf("budget: scale=2^%d maxLevel=%d", params.LogDefaultScale(), params.MaxLevel())
			}
			if params.LevelsConsumedPerRescaling() != 1 {
				t.Fatalf("levels per rescale=%d want=1", params.LevelsConsumedPerRescaling())
			}
			if got := params.LogQi(); !reflect.DeepEqual(got, tc.logQ) {
				t.Fatalf("logQ=%v want=%v", got, tc.logQ)
			}
			if got := params.LogPi(); !reflect.DeepEqual(got, tc.logP) {
				t.Fatalf("logP=%v want=%v", got, tc.logP)
			}
			if logQP := params.LogQP(); math.Abs(logQP-tc.logQP) > 0.1 {
				t.Fatalf("logQP %.3f differs from nominal %.3f", logQP, tc.logQP)
			}

			minRefreshLevel, _, ok := mpckks.GetMinimumLevelForRefresh(
				tc.security, params.DefaultScale(), tc.parties, params.Q())
			if !ok {
				t.Fatal("no collective-refresh level satisfies the requested mask bound")
			}
			postDepthLevel := params.MaxLevel() - tc.depth
			if postDepthLevel < minRefreshLevel {
				t.Fatalf("depth-%d leaves level %d, below refresh minimum %d", tc.depth, postDepthLevel, minRefreshLevel)
			}
		})
	}
}
