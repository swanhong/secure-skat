package workflow

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/hhcho/sfgwas/rewrite/protocol"
	"gonum.org/v1/gonum/mat"
)

const (
	auxiliaryPartyID = iota
	cohortAPartyID
	cohortBPartyID
)

type ChromosomeInput struct {
	Chromosome        int
	GenotypeDirectory string
	Genes             []protocol.Gene
	CryptoParams      protocol.CryptoParams
}

type PartyInput struct {
	X           *mat.Dense
	Y           *mat.Dense
	SampleCount int
	Chromosomes []ChromosomeInput
}

func readTextMatrix(path string, columns int) (*mat.Dense, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}

	fields := strings.Fields(string(data))
	values := make([]float64, len(fields))
	for index, field := range fields {
		values[index], err = strconv.ParseFloat(field, 64)
		if err != nil {
			return nil, err
		}
	}
	return mat.NewDense(len(values)/columns, columns, values), nil
}

func readRows(
	directory string,
	covariateCount, phenotypeCount int,
) (*mat.Dense, *mat.Dense, error) {
	covariates, err := readTextMatrix(directory+"/cov.txt", covariateCount)
	if err != nil {
		return nil, nil, err
	}
	y, err := readTextMatrix(directory+"/pheno.txt", phenotypeCount)
	if err != nil {
		return nil, nil, err
	}

	rows, _ := covariates.Dims()
	x := mat.NewDense(rows, covariateCount+1, nil)
	for row := 0; row < rows; row++ {
		x.Set(row, 0, 1)
		for col := 0; col < covariateCount; col++ {
			x.Set(row, col+1, covariates.At(row, col))
		}
	}
	return x, y, nil
}

func readGenes(directory string) ([]protocol.Gene, error) {
	geneData, err := os.ReadFile(filepath.Join(directory, "genes.txt"))
	if err != nil {
		return nil, err
	}
	countData, err := os.ReadFile(filepath.Join(directory, "block_sizes.txt"))
	if err != nil {
		return nil, err
	}

	geneIDs := strings.Fields(string(geneData))
	counts := strings.Fields(string(countData))
	if len(geneIDs) != len(counts) {
		return nil, fmt.Errorf(
			"genes.txt has %d genes but block_sizes.txt has %d entries",
			len(geneIDs),
			len(counts),
		)
	}
	genes := make([]protocol.Gene, len(geneIDs))
	for index, geneID := range geneIDs {
		variantCount, err := strconv.Atoi(counts[index])
		if err != nil {
			return nil, err
		}
		genes[index] = protocol.Gene{
			GeneID:       geneID,
			VariantCount: variantCount,
		}
	}
	return genes, nil
}

func decodeInt8(data []byte) []float64 {
	values := make([]float64, len(data))
	for index, value := range data {
		values[index] = float64(int8(value))
	}
	return values
}

func readInt8Matrix(
	path string,
	rows, columns int,
) (*mat.Dense, error) {
	if columns == 0 {
		return nil, nil
	}

	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	return mat.NewDense(
		rows,
		columns,
		decodeInt8(data),
	), nil
}

func readPrivateInt8Matrix(
	path string,
	rows int,
) (*mat.Dense, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	if len(data) == 0 {
		return nil, nil
	}

	values := decodeInt8(data)
	return mat.NewDense(
		rows,
		len(values)/rows,
		values,
	), nil
}

func readGeneBatch(
	directory string,
	genes []protocol.Gene,
	batch protocol.GeneBatch,
	rows int,
	includePrivate bool,
) (public, private []*mat.Dense, err error) {
	if directory == "" {
		return nil, nil, nil
	}

	if includePrivate {
		private = make([]*mat.Dense, len(genes))
	}
	public = make([]*mat.Dense, len(genes))

	for _, geneIndex := range batch.GeneIndices {
		gene := genes[geneIndex]
		block := fmt.Sprintf("block.%d.bin", geneIndex)

		public[geneIndex], err = readInt8Matrix(
			filepath.Join(directory, "geno", block),
			rows,
			gene.VariantCount,
		)
		if err != nil {
			return nil, nil, err
		}

		if includePrivate {
			private[geneIndex], err = readPrivateInt8Matrix(
				filepath.Join(directory, "private", block),
				rows,
			)
			if err != nil {
				return nil, nil, err
			}
		}
	}
	return public, private, nil
}

func readPublicDosageSums(
	directory string,
	genes []protocol.Gene,
	rows int,
) (*mat.Dense, error) {
	total := 0
	for _, gene := range genes {
		total += gene.VariantCount
	}
	if total == 0 {
		return nil, nil
	}

	dosageSums := make([]float64, total)
	offset := 0
	for geneIndex, gene := range genes {
		if gene.VariantCount == 0 {
			continue
		}

		data, err := os.ReadFile(filepath.Join(
			directory,
			"geno",
			fmt.Sprintf("block.%d.bin", geneIndex),
		))
		if err != nil {
			return nil, err
		}

		for row := 0; row < rows; row++ {
			rowOffset := row * gene.VariantCount
			for variant := 0; variant < gene.VariantCount; variant++ {
				dosageSums[offset+variant] +=
					float64(int8(data[rowOffset+variant]))
			}
		}
		offset += gene.VariantCount
	}

	return mat.NewDense(1, total, dosageSums), nil
}

func loadPartyInput(
	config *Config,
	partyID int,
	ancestry string,
) (*PartyInput, error) {
	cohort := ""
	switch partyID {
	case auxiliaryPartyID:
	case cohortAPartyID:
		cohort = "A"
	case cohortBPartyID:
		cohort = "B"
	default:
		return nil, fmt.Errorf("unknown party %d", partyID)
	}

	input := &PartyInput{
		Chromosomes: make(
			[]ChromosomeInput,
			len(config.Chromosomes),
		),
	}
	for index, chromosome := range config.Chromosomes {
		input.Chromosomes[index].Chromosome = chromosome
	}

	if partyID == auxiliaryPartyID {
		return input, nil
	}

	firstDirectory := filepath.Join(
		config.RunDir,
		"prepared",
		ancestry,
		fmt.Sprintf("chr%d", config.Chromosomes[0]),
		cohort,
	)
	var err error
	input.X, input.Y, err = readRows(
		firstDirectory,
		config.NumCov,
		len(config.PhenotypeColumns),
	)
	if err != nil {
		return nil, err
	}
	input.SampleCount, _ = input.X.Dims()

	for index, chromosome := range config.Chromosomes {
		directory := filepath.Join(
			config.RunDir,
			"prepared",
			ancestry,
			fmt.Sprintf("chr%d", chromosome),
		)

		genes, err := readGenes(directory)
		if err != nil {
			return nil, err
		}

		input.Chromosomes[index] = ChromosomeInput{
			Chromosome:        chromosome,
			GenotypeDirectory: filepath.Join(directory, cohort),
			Genes:             genes,
		}
	}

	return input, nil
}
