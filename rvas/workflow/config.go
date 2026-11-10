package workflow

import (
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"

	"github.com/BurntSushi/toml"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
)

const (
	globalConfigFilename  = "configGlobal.toml"
	prepareConfigFilename = "configPrepare.toml"
)

type chromosomeList []int

func (chromosomes *chromosomeList) UnmarshalTOML(value interface{}) error {
	values, ok := value.([]interface{})
	if !ok || len(values) == 0 {
		return fmt.Errorf("chromosomes must be a nonempty integer array or [\"all\"]")
	}
	if len(values) == 1 && values[0] == "all" {
		all := make(chromosomeList, 22)
		for i := range all {
			all[i] = i + 1
		}
		*chromosomes = all
		return nil
	}
	selected := make(chromosomeList, len(values))
	var seen [23]bool
	for i, value := range values {
		chromosome, ok := value.(int64)
		if !ok || chromosome < 1 || chromosome > 22 {
			return fmt.Errorf("chromosomes must contain integers from 1 to 22 or exactly [\"all\"]")
		}
		if seen[chromosome] {
			return fmt.Errorf("duplicate chromosome: %d", chromosome)
		}
		seen[chromosome] = true
		selected[i] = int(chromosome)
	}
	*chromosomes = selected
	return nil
}

type Config struct {
	RunDir        string         `toml:"run_dir"`
	PublicVarList string         `toml:"public_var_list"`
	Chromosomes   chromosomeList `toml:"chromosomes"`

	Genotype   string `toml:"genotype"`
	GenePanel  string `toml:"gene_panel"`
	Annotation string `toml:"annotation"`
	Phenotype  string `toml:"phenotype"`
	Covariate  string `toml:"covariate"`
	Ancestry   string `toml:"ancestry"`

	PhenotypeIDColumn string    `toml:"phenotype_id_column"`
	CovariateIDColumn string    `toml:"covariate_id_column"`
	CovariateColumn   string    `toml:"covariate_column"`
	CovariateColumns  []string  `toml:"covariate_columns"`
	IsCovSingleColumn bool      `toml:"is_cov_single_column"`
	AncestryIDColumn  string    `toml:"ancestry_id_column"`
	AncestryColumn    string    `toml:"ancestry_column"`
	PhenotypeColumns  []string  `toml:"phenotype_columns"`
	Ancestries        []string  `toml:"ancestries"`
	GlobalMaxMAF      []float64 `toml:"global_max_maf"`
	Masks             []string  `toml:"masks"`
	LocalMaxMAF       *float64  `toml:"local_max_maf"`

	SamplesPerCohort int64   `toml:"samples_per_cohort"`
	SharedRate       float64 `toml:"shared_rate"`

	DebugSplitNumParties int `toml:"debug_split_num_parties"`

	CKKS           string `toml:"ckks"`
	MpcNumThreads  int    `toml:"mpc_num_threads"`
	DataBits       int    `toml:"data_bits"`
	FractionalBits int    `toml:"fractional_bits"`
	Probes         int    `toml:"probes"`

	Plink2Bin     string        `toml:"plink2_bin"`
	GeneSelection GeneSelection `toml:"gene_selection"`

	BindingIP string                `toml:"binding_ipaddr"`
	Servers   map[string]mpc.Server `toml:"servers"`

	SharedKeysPath           string `toml:"shared_keys_path"`
	LocalNumThreads          int    `toml:"local_num_threads"`
	GenotypeDirectory        string `toml:"genotype_dir"`
	PrivateGenotypeDirectory string `toml:"private_genotype_dir"`
	PhenotypeFile            string `toml:"phenotype_file"`
	CovariateFile            string `toml:"covariate_file"`
	GenesFile                string `toml:"genes_file"`
	VariantCountsFile        string `toml:"variant_counts_file"`
	PrintLogs                *bool  `toml:"print_logs"`
}

type GeneSelection struct {
	Mode          string `toml:"mode"`
	PerChromosome int    `toml:"per_chromosome"`
	Path          string `toml:"path"`
}

func loadConfig(directory string, filenames ...string) (*Config, error) {
	config := &Config{SharedRate: 0.6, IsCovSingleColumn: true}
	for _, filename := range filenames {
		path := filepath.Join(directory, filename)
		metadata, err := toml.DecodeFile(path, config)
		if err != nil {
			return nil, fmt.Errorf("decode config %q: %w", path, err)
		}
		if filename != prepareConfigFilename && metadata.IsDefined("run_dir") {
			return nil, fmt.Errorf("%s: run_dir belongs in %s", path, prepareConfigFilename)
		}
		if metadata.IsDefined("max_maf") {
			return nil, fmt.Errorf("%s: max_maf was replaced by local_max_maf (cut by the site cohort's genotype MAF)", path)
		}
	}
	for index, ancestry := range config.Ancestries {
		config.Ancestries[index] = strings.ToUpper(strings.TrimSpace(ancestry))
	}
	for _, path := range []*string{
		&config.SharedKeysPath,
		&config.GenotypeDirectory,
		&config.PrivateGenotypeDirectory,
		&config.PhenotypeFile,
		&config.CovariateFile,
		&config.GenesFile,
		&config.VariantCountsFile,
	} {
		*path = strings.ReplaceAll(*path, "{run_dir}", config.RunDir)
	}
	return config, nil
}

func LoadPrepareConfig(directory string) (*Config, error) {
	filenames := []string{globalConfigFilename}
	if _, err := os.Stat(filepath.Join(directory, "configLocal.toml")); err == nil {
		filenames = append(filenames, "configLocal.toml")
	} else if !os.IsNotExist(err) {
		return nil, err
	}
	return loadConfig(directory, append(filenames, prepareConfigFilename)...)
}

func LoadPartyConfig(directory string, partyID int) (*Config, error) {
	localFilename := "configLocal.toml"
	if _, err := os.Stat(filepath.Join(directory, localFilename)); os.IsNotExist(err) {
		localFilename = fmt.Sprintf("configLocal.Party%d.toml", partyID)
	}
	return loadConfig(directory, globalConfigFilename, localFilename, prepareConfigFilename)
}

func requireStrings(fields ...string) error {
	for index := 0; index < len(fields); index += 2 {
		if strings.TrimSpace(fields[index+1]) == "" {
			return fmt.Errorf("%s is required", fields[index])
		}
	}
	return nil
}

func validateGlobalConfig(config *Config) error {
	if err := requireStrings("ckks", config.CKKS); err != nil {
		return err
	}
	if len(config.Chromosomes) == 0 || len(config.Ancestries) == 0 || len(config.PhenotypeColumns) == 0 {
		return fmt.Errorf("chromosomes, ancestries, and phenotype_columns are required")
	}
	if len(config.CovariateColumns) == 0 {
		return fmt.Errorf("covariate_columns is required")
	}
	for index, threshold := range config.GlobalMaxMAF {
		if threshold <= 0 || threshold > 0.5 || slices.Contains(config.GlobalMaxMAF[:index], threshold) {
			return fmt.Errorf("global_max_maf must contain distinct values in (0, 0.5]")
		}
	}
	seen := make(map[string]bool, len(config.CovariateColumns))
	for _, column := range config.CovariateColumns {
		if strings.TrimSpace(column) == "" || seen[column] {
			return fmt.Errorf("covariate_columns must contain distinct nonblank names")
		}
		seen[column] = true
	}
	if config.MpcNumThreads < 1 || config.DataBits < 1 ||
		config.FractionalBits < 1 || config.FractionalBits > config.DataBits || config.Probes < 1 {
		return fmt.Errorf("invalid protocol dimensions")
	}
	if !securecrypto.IsPackedSKATParameters(config.CKKS) {
		return fmt.Errorf("unsupported packed SKAT CKKS parameters %q", config.CKKS)
	}
	return nil
}

func validatePrepareConfig(config *Config) error {
	if err := validateGlobalConfig(config); err != nil {
		return err
	}
	if config.DebugSplitNumParties != 0 && config.DebugSplitNumParties != 2 {
		return fmt.Errorf("debug_split_num_parties must be 0 or 2")
	}
	if config.SharedRate < 0 || config.SharedRate > 1 {
		return fmt.Errorf("shared_rate must be between 0 and 1")
	}
	if config.LocalMaxMAF != nil && (*config.LocalMaxMAF <= 0 || *config.LocalMaxMAF > 0.5) {
		return fmt.Errorf("local_max_maf must be in (0, 0.5]")
	}
	if config.IsCovSingleColumn && strings.TrimSpace(config.CovariateColumn) == "" {
		return fmt.Errorf("covariate_column is required")
	}
	return requireStrings(
		"run_dir", config.RunDir,
		"genotype", config.Genotype,
		"gene_panel", config.GenePanel,
		"annotation", config.Annotation,
		"phenotype", config.Phenotype,
		"covariate", config.Covariate,
		"ancestry", config.Ancestry,
		"phenotype_id_column", config.PhenotypeIDColumn,
		"covariate_id_column", config.CovariateIDColumn,
		"ancestry_id_column", config.AncestryIDColumn,
		"ancestry_column", config.AncestryColumn,
		"plink2_bin", config.Plink2Bin,
	)
}

func validatePartyConfig(config *Config, partyID int) error {
	if err := validateGlobalConfig(config); err != nil {
		return err
	}
	if err := requireStrings(
		"run_dir", config.RunDir,
		"shared_keys_path", config.SharedKeysPath,
	); err != nil {
		return err
	}
	if config.LocalNumThreads < 1 {
		return fmt.Errorf("local_num_threads must be positive")
	}
	if partyID == auxiliaryPartyID {
		return nil
	}
	if partyID != cohortAPartyID && partyID != cohortBPartyID {
		return fmt.Errorf("unknown party %d", partyID)
	}
	if err := requireStrings(
		"genotype_dir", config.GenotypeDirectory,
		"phenotype_file", config.PhenotypeFile,
		"covariate_file", config.CovariateFile,
		"genes_file", config.GenesFile,
		"variant_counts_file", config.VariantCountsFile,
	); err != nil {
		return err
	}
	if partyID == cohortBPartyID {
		return requireStrings(
			"private_genotype_dir",
			config.PrivateGenotypeDirectory,
		)
	}
	return nil
}

func writeConfigFiles(directory, runDir string, filenames ...string) error {
	destination := filepath.Join(runDir, "config")
	if err := os.MkdirAll(destination, 0o755); err != nil {
		return err
	}
	for _, filename := range filenames {
		contents, err := os.ReadFile(filepath.Join(directory, filename))
		if err != nil {
			return err
		}
		if err := os.WriteFile(filepath.Join(destination, filename), contents, 0o644); err != nil {
			return err
		}
	}
	return nil
}
