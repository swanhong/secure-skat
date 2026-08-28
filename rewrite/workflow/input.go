package workflow

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/hhcho/sfgwas/rewrite/protocol"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/mat"
)

const (
	auxiliaryPartyID = iota
	cohortAPartyID
	cohortBPartyID
)

type ChromosomeInput struct {
	Chromosome       int
	Genes            []protocol.Gene
	CryptoParams     protocol.CryptoParams
	PublicGenotypes  []*mat.Dense
	PrivateGenotypes []*mat.Dense
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

func readGenotypes(
	directory string,
	genes []protocol.Gene,
	rows int,
	includePrivate bool,
) (public, private []*mat.Dense, err error) {
	if includePrivate {
		private = make([]*mat.Dense, len(genes))
	}
	public = make([]*mat.Dense, len(genes))

	for index, gene := range genes {
		block := fmt.Sprintf("block.%d.bin", index)

		public[index], err = readInt8Matrix(
			filepath.Join(directory, "geno", block),
			rows,
			gene.VariantCount,
		)
		if err != nil {
			return nil, nil, err
		}

		if includePrivate {
			private[index], err = readPrivateInt8Matrix(
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

func loadPartyInput(
	config *Config,
	partyID int,
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

	literal, err := securecrypto.ResolveCKKSParametersLiteral(
		config.CKKS,
	)
	if err != nil {
		return nil, err
	}
	heParameters, err := ckks.NewParametersFromLiteral(literal)
	if err != nil {
		return nil, err
	}

	input := &PartyInput{
		Chromosomes: make(
			[]ChromosomeInput,
			len(config.Chromosomes),
		),
	}

	if cohort != "" {
		firstDirectory := filepath.Join(
			config.RunDir,
			"prepared",
			fmt.Sprintf("chr%d", config.Chromosomes[0]),
			cohort,
		)
		input.X, input.Y, err = readRows(
			firstDirectory,
			len(config.CovariateColumns),
			len(config.PhenotypeColumns),
		)
		if err != nil {
			return nil, err
		}
		input.SampleCount, _ = input.X.Dims()
	}

	for index, chromosome := range config.Chromosomes {
		directory := filepath.Join(
			config.RunDir,
			"prepared",
			fmt.Sprintf("chr%d", chromosome),
		)

		genes, err := readGenes(directory)
		if err != nil {
			return nil, err
		}

		dataParams := protocol.DataParams{
			Genes:          genes,
			C:              len(config.CovariateColumns) + 1,
			PhenotypeCount: len(config.PhenotypeColumns),
		}
		cryptoParams, err := protocol.BuildCryptoParams(
			dataParams,
			config.Probes,
			heParameters.MaxSlots(),
		)
		if err != nil {
			return nil, err
		}

		var public, private []*mat.Dense
		if cohort != "" {
			public, private, err = readGenotypes(
				filepath.Join(directory, cohort),
				genes,
				input.SampleCount,
				partyID == cohortBPartyID,
			)
			if err != nil {
				return nil, err
			}
		}

		input.Chromosomes[index] = ChromosomeInput{
			Chromosome:       chromosome,
			Genes:            genes,
			CryptoParams:     cryptoParams,
			PublicGenotypes:  public,
			PrivateGenotypes: private,
		}
	}

	return input, nil
}
