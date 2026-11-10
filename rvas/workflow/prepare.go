package workflow

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
)

type prepareGeneSelection struct {
	Mode          string `json:"mode"`
	PerChromosome int    `json:"per_chromosome"`
	Path          string `json:"path"`
}

type prepareRequest struct {
	Mode             string    `json:"mode"`
	GlobalMaxMAF     []float64 `json:"global_max_maf"`
	Party            int       `json:"party"`
	PublicVarList    string    `json:"public_var_list"`
	RunDir           string    `json:"run_dir"`
	PreparedCacheDir string    `json:"prepared_cache_dir"`
	Chromosomes      []int     `json:"chromosomes"`

	Genotype   string `json:"genotype"`
	GenePanel  string `json:"gene_panel"`
	Annotation string `json:"annotation"`
	Phenotype  string `json:"phenotype"`
	Covariate  string `json:"covariate"`
	Ancestry   string `json:"ancestry"`
	Plink2Bin  string `json:"plink2_bin"`

	PhenotypeIDColumn string   `json:"phenotype_id_column"`
	CovariateIDColumn string   `json:"covariate_id_column"`
	CovariateColumn   string   `json:"covariate_column"`
	CovariateColumns  []string `json:"covariate_columns"`
	IsCovSingleColumn bool     `json:"is_cov_single_column"`
	AncestryIDColumn  string   `json:"ancestry_id_column"`
	AncestryColumn    string   `json:"ancestry_column"`
	PhenotypeColumns  []string `json:"phenotype_columns"`
	Ancestries        []string `json:"ancestries"`
	LocalMaxMAF       *float64 `json:"local_max_maf"`

	SamplesPerCohort int64   `json:"samples_per_cohort"`
	SharedRate       float64 `json:"shared_rate"`

	DebugSplitNumParties int `json:"debug_split_num_parties"`

	GeneSelection prepareGeneSelection `json:"gene_selection"`
}

func prepareRequestFromConfig(config *Config) prepareRequest {
	return prepareRequest{
		Mode:              "joint",
		GlobalMaxMAF:      config.GlobalMaxMAF,
		RunDir:            config.RunDir,
		PreparedCacheDir:  filepath.Join(config.RunDir, "cache", "prepared"),
		Chromosomes:       config.Chromosomes,
		Genotype:          config.Genotype,
		GenePanel:         config.GenePanel,
		Annotation:        config.Annotation,
		Phenotype:         config.Phenotype,
		Covariate:         config.Covariate,
		Ancestry:          config.Ancestry,
		Plink2Bin:         config.Plink2Bin,
		PhenotypeIDColumn: config.PhenotypeIDColumn,
		CovariateIDColumn: config.CovariateIDColumn,
		CovariateColumn:   config.CovariateColumn,
		CovariateColumns:  config.CovariateColumns,
		IsCovSingleColumn: config.IsCovSingleColumn,
		AncestryIDColumn:  config.AncestryIDColumn,
		AncestryColumn:    config.AncestryColumn,
		PhenotypeColumns:  config.PhenotypeColumns,
		Ancestries:        config.Ancestries,
		LocalMaxMAF:       config.LocalMaxMAF,
		SamplesPerCohort:  config.SamplesPerCohort,
		SharedRate:        config.SharedRate,
		GeneSelection: prepareGeneSelection{
			Mode:          config.GeneSelection.Mode,
			PerChromosome: config.GeneSelection.PerChromosome,
			Path:          config.GeneSelection.Path,
		},

		DebugSplitNumParties: config.DebugSplitNumParties,
	}
}

func clearGeneratedOutputs(runDir, mode string) error {
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

	names := []string{"prepared/local", "cache/prepared/local", "local"}
	if mode != "local" {
		if mode == "both" {
			names = []string{"prepared", "cache/prepared", "local"}
		} else {
			names = nil
			for _, parent := range []string{"prepared", "cache/prepared"} {
				entries, err := os.ReadDir(filepath.Join(root, parent))
				if err != nil && !os.IsNotExist(err) {
					return err
				}
				for _, entry := range entries {
					if entry.Name() != "local" {
						names = append(names, filepath.Join(parent, entry.Name()))
					}
				}
			}
		}
		names = append(names, "selected_genes.tsv",
			"secure", "reference", "comparison", "metrics")
	}
	for _, name := range names {
		if err := os.RemoveAll(filepath.Join(root, name)); err != nil {
			return fmt.Errorf("remove generated output %s: %w", name, err)
		}
	}
	return nil
}

func Prepare(config *Config) error {
	return prepareSite(config, 0, "")
}

func prepareSite(config *Config, party int, publicVarList string) error {
	return prepareSiteMode(config, party, publicVarList, "joint")
}

func prepareSiteMode(config *Config, party int, publicVarList, mode string) error {
	request := prepareRequestFromConfig(config)
	request.Mode = mode
	request.Party = party
	request.PublicVarList = publicVarList
	payload, err := json.Marshal(request)
	if err != nil {
		return fmt.Errorf("encode preprocessing request: %w", err)
	}

	command := exec.Command(
		"python3",
		"-m",
		"python.preprocessing.prepare",
	)
	command.Stdin = bytes.NewReader(payload)
	command.Stdout = os.Stdout
	command.Stderr = os.Stderr

	if err := command.Run(); err != nil {
		return fmt.Errorf("run preprocessing: %w", err)
	}

	return nil
}
