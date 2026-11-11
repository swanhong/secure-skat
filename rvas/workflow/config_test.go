package workflow

import (
	"math"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"
)

func TestPrepareRunDir(t *testing.T) {
	directory := t.TempDir()
	runDir := filepath.Join(directory, "run")
	files := map[string]string{
		"configGlobal.toml":  "",
		"configLocal.toml":   "shared_keys_path = \"{run_dir}/secure/keys\"\n",
		"configPrepare.toml": "run_dir = \"" + runDir + "\"\n",
	}
	write := func(name, contents string) {
		t.Helper()
		if err := os.WriteFile(filepath.Join(directory, name), []byte(contents), 0o600); err != nil {
			t.Fatal(err)
		}
	}
	for name, contents := range files {
		write(name, contents)
	}
	prepared, err := LoadPrepareConfig(directory)
	if err != nil {
		t.Fatal(err)
	}
	party, err := LoadPartyConfig(directory, 1)
	if err != nil {
		t.Fatal(err)
	}
	if prepared.RunDir != runDir || party.RunDir != runDir || party.SharedKeysPath != runDir+"/secure/keys" {
		t.Fatal("prepare and party must use run_dir from Prepare")
	}
	if prepareRequestFromConfig(prepared).PreparedCacheDir != filepath.Join(runDir, "cache", "prepared") {
		t.Fatal("cache must be inside run_dir")
	}
	for _, name := range []string{"configGlobal.toml", "configLocal.toml"} {
		write(name, files[name]+"run_dir = \"old-path\"\n")
		if _, err := LoadPrepareConfig(directory); err == nil || !strings.Contains(err.Error(), "run_dir belongs in") {
			t.Fatalf("expected migration error for %s, got %v", name, err)
		}
		write(name, files[name])
	}
}

func TestTraceSeedDerivation(t *testing.T) {
	key := []byte("01234567890123456789012345678901")
	base := deriveTraceSeed(key, "EUR", 1)
	if deriveTraceSeed(key, "EUR", 1) != base {
		t.Fatal("derivation not deterministic")
	}
	for _, other := range []int64{
		deriveTraceSeed(key, "AFR", 1),
		deriveTraceSeed(key, "EUR", 2),
		deriveTraceSeed([]byte("different key"), "EUR", 1),
	} {
		if other == base {
			t.Fatal("derivation did not separate contexts")
		}
	}
}

func TestCovariateColumnValidation(t *testing.T) {
	config, err := LoadPartyConfig(filepath.Join("..", "..", "config", "mvp"), 1)
	if err != nil {
		t.Fatal(err)
	}
	for _, test := range []struct {
		name    string
		columns []string
		valid   bool
	}{
		{name: "selected columns", columns: []string{"AGE", "PC4"}, valid: true},
		{name: "missing"},
		{name: "blank", columns: []string{"PC1", "  "}},
		{name: "duplicate", columns: []string{"PC1", "PC1"}},
	} {
		t.Run(test.name, func(t *testing.T) {
			candidate := *config
			candidate.CovariateColumns = test.columns
			err := validatePartyConfig(&candidate, cohortAPartyID)
			if (err == nil) != test.valid {
				t.Fatalf("valid=%v, got %v", test.valid, err)
			}
			if prepareErr := validatePrepareConfig(&candidate); (prepareErr == nil) != test.valid {
				t.Fatalf("prepare valid=%v, got %v", test.valid, prepareErr)
			}
			if err != nil && !strings.Contains(err.Error(), "covariate_columns") {
				t.Fatalf("unexpected validation error: %v", err)
			}
		})
	}
}

func TestHelperWithoutCovariateColumns(t *testing.T) {
	config, err := LoadPartyConfig(filepath.Join("..", "..", "config", "cp0"), auxiliaryPartyID)
	if err != nil {
		t.Fatal(err)
	}
	if len(config.CovariateColumns) != 0 {
		t.Fatal("cp0 must not configure covariate columns")
	}
	if err := validatePartyConfig(config, auxiliaryPartyID); err != nil {
		t.Fatal(err)
	}
}

func TestShippedCovariateConfigs(t *testing.T) {
	var schema []string
	for _, test := range []struct {
		directory string
		party     int
	}{
		{"1kg-A", 1}, {"1kg-B", 2}, {"1kg-cp0", 0},
		{"mvp", 1}, {"aou", 2}, {"cp0", 0},
		{"mvp-test", 1}, {"aou-test", 2}, {"cp0-test", 0},
	} {
		t.Run(test.directory, func(t *testing.T) {
			directory := filepath.Join("..", "..", "config", test.directory)
			config, err := LoadPartyConfig(directory, test.party)
			if err != nil {
				t.Fatal(err)
			}
			if err := validatePartyConfig(config, test.party); err != nil {
				t.Fatal(err)
			}
			if strings.HasPrefix(test.directory, "1kg-") {
				if len(config.CovariateColumns) != 10 {
					t.Fatalf("expected ten selected PCs, got %v", config.CovariateColumns)
				}
				if schema == nil {
					schema = config.CovariateColumns
				} else if !slices.Equal(schema, config.CovariateColumns) {
					t.Fatalf("party covariate schema differs: %v", config.CovariateColumns)
				}
			}
			if test.party == auxiliaryPartyID {
				return
			}
			prepared, err := LoadPrepareConfig(directory)
			if err != nil {
				t.Fatal(err)
			}
			if err := validatePrepareConfig(prepared); err != nil {
				t.Fatal(err)
			}
			request := prepareRequestFromConfig(prepared)
			if !slices.Equal(request.CovariateColumns, config.CovariateColumns) {
				t.Fatal("preparation and runtime selected different covariates")
			}
			if strings.HasPrefix(test.directory, "1kg-") &&
				(!prepared.IsCovSingleColumn || prepared.CovariateColumn != "pca_features") {
				t.Fatal("1KG must retain the packed PCA input format")
			}
		})
	}
}

func TestPreparedCovariateDimensions(t *testing.T) {
	directory := t.TempDir()
	write := func(name, contents string) string {
		t.Helper()
		path := filepath.Join(directory, name)
		if err := os.WriteFile(path, []byte(contents), 0o600); err != nil {
			t.Fatal(err)
		}
		return path
	}
	config := &Config{
		Chromosomes:       []int{1},
		CovariateColumns:  []string{"AGE", "PC4"},
		PhenotypeColumns:  []string{"trait"},
		CovariateFile:     write("cov.txt", "\n20\t0.25\n\n40\t-0.5\n"),
		PhenotypeFile:     write("pheno.txt", "1\n2\n"),
		GenesFile:         write("genes.txt", "gene1\n"),
		VariantCountsFile: write("block_sizes.txt", "1\n"),
		GenotypeDirectory: filepath.Join(directory, "geno"),
		GlobalMaxMAF:      []float64{0.01},
	}
	write("cohort_size.txt", "2\n")
	write("allele_counts.txt", "1\n")
	input, err := loadPartyInput(config, cohortAPartyID, "EUR")
	if err != nil {
		t.Fatal(err)
	}
	rows, columns := input.X.Dims()
	if rows != 2 || columns != 3 || input.SampleCount != 2 {
		t.Fatalf("wrong shape after adding intercept: %dx%d, samples=%d", rows, columns, input.SampleCount)
	}
	for i, want := range [][]float64{{1, 20, 0.25}, {1, 40, -0.5}} {
		if !slices.Equal(input.X.RawRowView(i), want) {
			t.Fatalf("covariate row %d: got %v want %v", i, input.X.RawRowView(i), want)
		}
	}

	// Six tokens are divisible by the selected width two, but these are two
	// three-column rows; flattening them would silently create a third sample.
	config.CovariateFile = write("old-cov.txt", "\n1 2 3\n4 5 6\n")
	_, err = loadPartyInput(config, cohortAPartyID, "EUR")
	if err == nil || !strings.Contains(err.Error(), config.CovariateFile+":2:") ||
		!strings.Contains(err.Error(), "expected 2 columns, got 3") {
		t.Fatalf("expected row-width error for stale prepared data, got %v", err)
	}
}

func TestGlobalMAFRequired(t *testing.T) {
	config, err := LoadPartyConfig(filepath.Join("..", "..", "config", "cp0"), 0)
	if err != nil {
		t.Fatal(err)
	}
	for _, thresholds := range [][]float64{nil, {}, {0}, {1}, {math.NaN()}, {0.01, 0.01}} {
		candidate := *config
		candidate.GlobalMaxMAF = thresholds
		if err := validateGlobalConfig(&candidate); err == nil {
			t.Fatalf("accepted invalid thresholds: %v", thresholds)
		}
	}
	for _, thresholds := range [][]float64{{0.01, 0.001}, {0.5}} {
		candidate := *config
		candidate.GlobalMaxMAF = thresholds
		if err := validateGlobalConfig(&candidate); err != nil {
			t.Fatal(err)
		}
	}
}
