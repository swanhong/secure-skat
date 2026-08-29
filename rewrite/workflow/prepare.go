package workflow

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
)

type prepareGeneSelection struct {
	Mode          string `json:"mode"`
	PerChromosome int    `json:"per_chromosome"`
	Seed          int64  `json:"seed"`
	Path          string `json:"path"`
}

type prepareRequest struct {
	RunDir      string `json:"run_dir"`
	Chromosomes []int  `json:"chromosomes"`

	Genotype   string `json:"genotype"`
	GenePanel  string `json:"gene_panel"`
	Annotation string `json:"annotation"`
	Phenotype  string `json:"phenotype"`
	Covariate  string `json:"covariate"`
	Plink2Bin  string `json:"plink2_bin"`

	PhenotypeIDColumn string            `json:"phenotype_id_column"`
	CovariateIDColumn string            `json:"covariate_id_column"`
	PhenotypeColumns  []string          `json:"phenotype_columns"`
	CovariateColumns  []string          `json:"covariate_columns"`
	Mask              map[string]string `json:"mask"`
	MaxMAF            *float64          `json:"max_maf"`

	SamplesPerCohort int64 `json:"samples_per_cohort"`
	SampleSeed       int64 `json:"sample_seed"`
	RoleSeed         int64 `json:"role_seed"`

	GeneSelection prepareGeneSelection `json:"gene_selection"`
}

func prepareRequestFromConfig(config *Config) prepareRequest {
	mask := make(map[string]string, len(config.Masks))
	for _, item := range config.Masks {
		column, value, _ := strings.Cut(item, "=")
		mask[strings.TrimSpace(column)] = strings.TrimSpace(value)
	}

	return prepareRequest{
		RunDir:            config.RunDir,
		Chromosomes:       config.Chromosomes,
		Genotype:          config.Genotype,
		GenePanel:         config.GenePanel,
		Annotation:        config.Annotation,
		Phenotype:         config.Phenotype,
		Covariate:         config.Covariate,
		Plink2Bin:         config.Plink2Bin,
		PhenotypeIDColumn: config.PhenotypeIDColumn,
		CovariateIDColumn: config.CovariateIDColumn,
		PhenotypeColumns:  config.PhenotypeColumns,
		CovariateColumns:  config.CovariateColumns,
		Mask:              mask,
		MaxMAF:            config.MaxMAF,
		SamplesPerCohort:  config.SamplesPerCohort,
		SampleSeed:        config.SampleSeed,
		RoleSeed:          config.RoleSeed,
		GeneSelection: prepareGeneSelection{
			Mode:          config.GeneSelection.Mode,
			PerChromosome: config.GeneSelection.PerChromosome,
			Seed:          config.GeneSelection.Seed,
			Path:          config.GeneSelection.Path,
		},
	}
}

func clearGeneratedOutputs(runDir string) error {
	root, err := filepath.Abs(runDir)
	if err != nil {
		return fmt.Errorf("resolve run directory: %w", err)
	}
	if resolved, err := filepath.EvalSymlinks(root); err == nil {
		root = resolved
	} else if !os.IsNotExist(err) {
		return fmt.Errorf("resolve run directory symlinks: %w", err)
	}

	workingDirectory, err := os.Getwd()
	if err != nil {
		return fmt.Errorf("get working directory: %w", err)
	}
	homeDirectory, err := os.UserHomeDir()
	if err != nil {
		return fmt.Errorf("get home directory: %w", err)
	}
	volumeRoot := filepath.VolumeName(root) + string(filepath.Separator)
	if root == volumeRoot || root == workingDirectory || root == homeDirectory {
		return fmt.Errorf("refuse to clear unsafe run directory %s", root)
	}

	for _, name := range []string{
		"selected_genes.tsv",
		"prepared",
		"secure",
		"reference",
		"comparison",
		"metrics",
	} {
		if err := os.RemoveAll(filepath.Join(root, name)); err != nil {
			return fmt.Errorf("remove generated output %s: %w", name, err)
		}
	}
	return nil
}

func Prepare(config *Config) error {
	if err := ValidateConfig(config); err != nil {
		return fmt.Errorf("validate config: %w", err)
	}

	payload, err := json.Marshal(
		prepareRequestFromConfig(config),
	)
	if err != nil {
		return fmt.Errorf("encode preprocessing request: %w", err)
	}

	command := exec.Command(
		"python3",
		"-m",
		"rewrite.preprocessing.prepare",
	)
	command.Stdin = bytes.NewReader(payload)
	command.Stdout = os.Stdout
	command.Stderr = os.Stderr

	if err := command.Run(); err != nil {
		return fmt.Errorf("run preprocessing: %w", err)
	}

	return nil
}
