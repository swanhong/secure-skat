package workflow

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/BurntSushi/toml"
	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/mpc"
)

const (
	globalConfigFilename  = "configGlobal.toml"
	prepareConfigFilename = "configPrepare.toml"
)

type Config struct {
	RunDir      string `toml:"run_dir"`
	Chromosomes []int  `toml:"chromosomes"`

	Genotype   string `toml:"genotype"`
	GenePanel  string `toml:"gene_panel"`
	Annotation string `toml:"annotation"`
	Phenotype  string `toml:"phenotype"`
	Covariate  string `toml:"covariate"`
	Ancestry   string `toml:"ancestry"`

	PhenotypeIDColumn string   `toml:"phenotype_id_column"`
	CovariateIDColumn string   `toml:"covariate_id_column"`
	CovariateColumn   string   `toml:"covariate_column"`
	AncestryIDColumn  string   `toml:"ancestry_id_column"`
	AncestryColumn    string   `toml:"ancestry_column"`
	PhenotypeColumns  []string `toml:"phenotype_columns"`
	Ancestries        []string `toml:"ancestries"`
	NumCov            int      `toml:"num_cov"`
	Masks             []string `toml:"masks"`
	MaxMAF            *float64 `toml:"max_maf"`

	SamplesPerCohort int64 `toml:"samples_per_cohort"`
	SampleSeed       int64 `toml:"sample_seed"`
	RoleSeed         int64 `toml:"role_seed"`

	CKKS           string `toml:"ckks"`
	MpcNumThreads  int    `toml:"mpc_num_threads"`
	DataBits       int    `toml:"data_bits"`
	FractionalBits int    `toml:"fractional_bits"`
	Probes         int    `toml:"probes"`
	Seed           int64  `toml:"seed"`

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
}

type GeneSelection struct {
	Mode          string `toml:"mode"`
	PerChromosome int    `toml:"per_chromosome"`
	Seed          int64  `toml:"seed"`
	Path          string `toml:"path"`
}

func loadConfig(directory string, filenames ...string) (*Config, error) {
	config := new(Config)
	for _, filename := range filenames {
		path := filepath.Join(directory, filename)
		if _, err := toml.DecodeFile(path, config); err != nil {
			return nil, fmt.Errorf("decode config %q: %w", path, err)
		}
	}
	for index, ancestry := range config.Ancestries {
		config.Ancestries[index] = strings.ToUpper(strings.TrimSpace(ancestry))
	}
	return config, nil
}

func LoadPrepareConfig(directory string) (*Config, error) {
	return loadConfig(directory, globalConfigFilename, prepareConfigFilename)
}

func LoadPartyConfig(directory string, partyID int) (*Config, error) {
	return loadConfig(
		directory,
		globalConfigFilename,
		fmt.Sprintf("configLocal.Party%d.toml", partyID),
	)
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
	if config.NumCov < 1 || config.MpcNumThreads < 1 || config.DataBits < 1 ||
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
		"covariate_column", config.CovariateColumn,
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
