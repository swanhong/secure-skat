package workflow

import (
	"fmt"
	"strings"

	"github.com/BurntSushi/toml"
	securecrypto "github.com/hhcho/sfgwas/crypto"
)

type Config struct {
	RunDir      string `toml:"run_dir"`
	Chromosomes []int  `toml:"chromosomes"`

	Genotype   string `toml:"genotype"`
	GenePanel  string `toml:"gene_panel"`
	Annotation string `toml:"annotation"`
	Phenotype  string `toml:"phenotype"`
	Covariate  string `toml:"covariate"`

	PhenotypeIDColumn    string   `toml:"phenotype_id_column"`
	CovariateIDColumn    string   `toml:"covariate_id_column"`
	CovariateArrayColumn string   `toml:"covariate_array_column"`
	PhenotypeColumns     []string `toml:"phenotype_columns"`
	CovariateColumns     []string `toml:"covariate_columns"`
	Masks                []string `toml:"masks"`

	SamplesPerCohort int64 `toml:"samples_per_cohort"`
	SampleSeed       int64 `toml:"sample_seed"`
	RoleSeed         int64 `toml:"role_seed"`

	CKKS           string `toml:"ckks"`
	DataBits       int    `toml:"data_bits"`
	FractionalBits int    `toml:"fractional_bits"`
	Probes         int    `toml:"probes"`
	Seed           int64  `toml:"seed"`
	PortBase       int    `toml:"port_base"`

	Plink2Bin string `toml:"plink2_bin"`

	GeneSelection GeneSelection `toml:"gene_selection"`
}

type GeneSelection struct {
	Mode          string `toml:"mode"`
	PerChromosome int    `toml:"per_chromosome"`
	Seed          int64  `toml:"seed"`
	Path          string `toml:"path"`
}

func LoadConfig(path string) (*Config, error) {
	config := new(Config)
	if _, err := toml.DecodeFile(path, config); err != nil {
		return nil, fmt.Errorf("decode config %q, %w", path, err)
	}
	return config, nil
}

func ValidateConfig(config *Config) error {
	requiredFields := []struct {
		name  string
		value string
	}{
		{"run_dir", config.RunDir},
		{"genotype", config.Genotype},
		{"gene_panel", config.GenePanel},
		{"annotation", config.Annotation},
		{"phenotype", config.Phenotype},
		{"covariate", config.Covariate},
		{"phenotype_id_column", config.PhenotypeIDColumn},
		{"covariate_id_column", config.CovariateIDColumn},
		{"ckks", config.CKKS},
		{"plink2_bin", config.Plink2Bin},
	}
	// Check required fields are non-empty
	for _, field := range requiredFields {
		if strings.TrimSpace(field.value) == "" {
			return fmt.Errorf("%s is required", field.name)
		}
	}

	if len(config.Chromosomes) == 0 {
		return fmt.Errorf("chromosomes is required")
	}

	seenChromosomes := make(map[int]bool, len(config.Chromosomes))
	for _, chromosome := range config.Chromosomes {
		if chromosome < 1 || chromosome > 22 {
			return fmt.Errorf(
				"chromosome %d must be between 1 and 22",
				chromosome,
			)
		}
		if seenChromosomes[chromosome] {
			return fmt.Errorf("duplicate chromosome %d", chromosome)
		}
		seenChromosomes[chromosome] = true
	}

	columnGroups := []struct {
		name    string
		columns []string
	}{
		{"phenotype_columns", config.PhenotypeColumns},
		{"covariate_columns", config.CovariateColumns},
	}
	for _, group := range columnGroups {
		// at least one column is required for pheno/cov
		if len(group.columns) == 0 {
			return fmt.Errorf("%s is required", group.name)
		}

		// check for empty or duplicate values
		seenColumns := make(map[string]bool, len(group.columns))
		for _, column := range group.columns {
			column = strings.TrimSpace(column)
			if column == "" {
				return fmt.Errorf("%s contains an empty value", group.name)
			}
			if seenColumns[column] {
				return fmt.Errorf(
					"duplicate %s value %q",
					group.name,
					column,
				)
			}
			seenColumns[column] = true
		}
	}

	if len(config.Masks) == 0 {
		return fmt.Errorf("masks is required")
	}

	seenMaskColumns := make(map[string]bool, len(config.Masks))
	for _, mask := range config.Masks {
		column, value, found := strings.Cut(mask, "=")
		column = strings.TrimSpace(column)
		value = strings.TrimSpace(value)

		if !found || column == "" || value == "" {
			return fmt.Errorf(
				"mask must be the format \"column=value\", got %q",
				mask,
			)
		}
		if seenMaskColumns[column] {
			return fmt.Errorf("duplicate mask column %q", column)
		}
		seenMaskColumns[column] = true
	}

	if config.SamplesPerCohort < 0 {
		return fmt.Errorf("samples_per_cohort must be non-negative")
	}

	if !securecrypto.IsPackedSKATParameters(config.CKKS) {
		return fmt.Errorf(
			"unsupported packed SKAT CKKS parameters %q",
			config.CKKS,
		)
	}

	if config.DataBits < 1 {
		return fmt.Errorf("data_bits must be positive")
	}
	if config.FractionalBits < 1 ||
		config.FractionalBits > config.DataBits {
		return fmt.Errorf(
			"fractional_bits must be between 1 and data_bits",
		)
	}
	if config.Probes < 1 {
		return fmt.Errorf("probes must be positive")
	}
	if config.PortBase < 1 || config.PortBase > 65533 {
		return fmt.Errorf("port_base must be between 1 and 65533")
	}

	switch config.GeneSelection.Mode {
	case "random":
		if config.GeneSelection.PerChromosome < 1 {
			return fmt.Errorf(
				"GeneSelection random -> gene_selection.per_chromosome must be positive",
			)
		}
	case "file":
		if strings.TrimSpace(config.GeneSelection.Path) == "" {
			return fmt.Errorf(
				"GeneSelection file -> gene_selection.path is required in file mode",
			)
		}
	case "all":
	default:
		return fmt.Errorf(
			"unsupported gene_selection.mode %q",
			config.GeneSelection.Mode,
		)
	}

	return nil
}
