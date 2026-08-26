package main

import (
	"bufio"
	"encoding/csv"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/hhcho/sfgwas/rewrite/protocol"
	"gonum.org/v1/gonum/mat"
)

const maxPublicVariants = 4096

type preprocessingMetadata struct {
	SampleCountA     int      `json:"sample_count_a"`
	SampleCountB     int      `json:"sample_count_b"`
	PhenotypeColumns []string `json:"phenotype_columns"`
	CovariateColumns []string `json:"covariate_columns"`
}

type geneMetadata struct {
	ID         string
	Symbol     string
	Chromosome string
	Order      int
}

type preprocessedInput struct {
	Root       string
	Metadata   preprocessingMetadata
	Genes      []geneMetadata
	DataParams protocol.DataParams
}

func loadPreprocessedInput(root string) (preprocessedInput, error) {
	var metadata preprocessingMetadata
	encoded, err := os.ReadFile(filepath.Join(root, "metadata.json"))
	if err != nil {
		return preprocessedInput{}, err
	}
	if err := json.Unmarshal(encoded, &metadata); err != nil {
		return preprocessedInput{}, fmt.Errorf("decode metadata.json: %w", err)
	}

	genes, err := readGeneMetadata(filepath.Join(root, "gene_metadata.tsv"))
	if err != nil {
		return preprocessedInput{}, err
	}
	blockSizes, err := readIntegers(filepath.Join(root, "block_sizes.txt"))
	if err != nil {
		return preprocessedInput{}, err
	}
	if len(genes) != len(blockSizes) {
		return preprocessedInput{}, fmt.Errorf(
			"gene metadata has %d genes but block_sizes.txt has %d",
			len(genes), len(blockSizes),
		)
	}

	protocolGenes := make([]protocol.Gene, len(genes))
	for index := range genes {
		if blockSizes[index] > maxPublicVariants {
			return preprocessedInput{}, fmt.Errorf(
				"gene %s has %d public variants, limit is %d",
				genes[index].ID, blockSizes[index], maxPublicVariants,
			)
		}
		protocolGenes[index] = protocol.Gene{
			GeneID:       genes[index].ID,
			VariantCount: blockSizes[index],
		}
	}

	covariateCount := len(metadata.CovariateColumns) + 1
	phenotypeCount := len(metadata.PhenotypeColumns)
	if metadata.SampleCountA < 1 || metadata.SampleCountB < 1 {
		return preprocessedInput{}, fmt.Errorf("both cohorts must contain samples")
	}
	if phenotypeCount < 1 {
		return preprocessedInput{}, fmt.Errorf("metadata has no phenotype columns")
	}

	return preprocessedInput{
		Root:     root,
		Metadata: metadata,
		Genes:    genes,
		DataParams: protocol.DataParams{
			Genes:          protocolGenes,
			NA:             metadata.SampleCountA,
			NB:             metadata.SampleCountB,
			N:              metadata.SampleCountA + metadata.SampleCountB,
			C:              covariateCount,
			PhenotypeCount: phenotypeCount,
		},
	}, nil
}

func readGeneMetadata(path string) ([]geneMetadata, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = '\t'
	header, err := reader.Read()
	if err != nil {
		return nil, err
	}
	columns := make(map[string]int, len(header))
	for index, name := range header {
		columns[name] = index
	}
	required := []string{"gene_id", "gene_symbol", "chromosome", "gene_order"}
	for _, name := range required {
		if _, ok := columns[name]; !ok {
			return nil, fmt.Errorf("gene_metadata.tsv missing column %s", name)
		}
	}

	genes := []geneMetadata{}
	for {
		row, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		geneOrder, err := strconv.Atoi(row[columns["gene_order"]])
		if err != nil {
			return nil, err
		}
		genes = append(genes, geneMetadata{
			ID:         row[columns["gene_id"]],
			Symbol:     row[columns["gene_symbol"]],
			Chromosome: row[columns["chromosome"]],
			Order:      geneOrder,
		})
	}
	return genes, nil
}

func readIntegers(path string) ([]int, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	values := []int{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		text := strings.TrimSpace(scanner.Text())
		if text == "" {
			continue
		}
		value, err := strconv.Atoi(text)
		if err != nil {
			return nil, err
		}
		values = append(values, value)
	}
	return values, scanner.Err()
}

func loadPartyPhenotypeInput(
	input preprocessedInput,
	party int,
) (x, y *mat.Dense, err error) {
	if party == auxiliaryPartyID {
		return nil, nil, nil
	}

	cohort := "A"
	rows := input.Metadata.SampleCountA
	if party == cohortBPartyID {
		cohort = "B"
		rows = input.Metadata.SampleCountB
	}
	covariates, err := readTextMatrix(
		filepath.Join(input.Root, cohort, "cov.txt"),
		rows,
		len(input.Metadata.CovariateColumns),
	)
	if err != nil {
		return nil, nil, err
	}
	phenotypes, err := readTextMatrix(
		filepath.Join(input.Root, cohort, "pheno.txt"),
		rows,
		len(input.Metadata.PhenotypeColumns),
	)
	if err != nil {
		return nil, nil, err
	}
	return addIntercept(covariates), phenotypes, nil
}

func loadPartyGeneBatch(
	input preprocessedInput,
	party int,
	geneIndices []int,
) (gp, gv []*mat.Dense, err error) {
	gp = make([]*mat.Dense, len(geneIndices))
	gv = make([]*mat.Dense, len(geneIndices))
	if party == auxiliaryPartyID {
		return gp, gv, nil
	}

	cohort := "A"
	rows := input.Metadata.SampleCountA
	if party == cohortBPartyID {
		cohort = "B"
		rows = input.Metadata.SampleCountB
	}

	for position, geneIndex := range geneIndices {
		publicVariantCount := input.DataParams.Genes[geneIndex].VariantCount
		blockName := fmt.Sprintf("block.%d.bin", geneIndex)
		gp[position], err = readDosageMatrix(
			filepath.Join(input.Root, cohort, "geno", blockName),
			rows,
			publicVariantCount,
		)
		if err != nil {
			return nil, nil, err
		}

		if party == cohortBPartyID {
			gv[position], err = readPrivateDosageMatrix(
				filepath.Join(input.Root, "B", "private", blockName), rows,
			)
			if err != nil {
				return nil, nil, err
			}
		}
	}
	return gp, gv, nil
}

func readTextMatrix(path string, rows, columns int) (*mat.Dense, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	values := make([]float64, 0, rows*columns)
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		if len(fields) != columns {
			return nil, fmt.Errorf("%s has row with %d columns, want %d", path, len(fields), columns)
		}
		for _, field := range fields {
			value, err := strconv.ParseFloat(field, 64)
			if err != nil {
				return nil, err
			}
			values = append(values, value)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if len(values) != rows*columns {
		return nil, fmt.Errorf("%s has %d values, want %d", path, len(values), rows*columns)
	}
	return mat.NewDense(rows, columns, values), nil
}

func readDosageMatrix(path string, rows, columns int) (*mat.Dense, error) {
	encoded, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	if len(encoded) != rows*columns {
		return nil, fmt.Errorf("%s has %d bytes, want %d", path, len(encoded), rows*columns)
	}
	if columns == 0 {
		return nil, nil
	}

	values := make([]float64, len(encoded))
	for index, value := range encoded {
		values[index] = float64(int8(value))
	}
	return mat.NewDense(rows, columns, values), nil
}

func readPrivateDosageMatrix(path string, rows int) (*mat.Dense, error) {
	info, err := os.Stat(path)
	if err != nil {
		return nil, err
	}
	if info.Size()%int64(rows) != 0 {
		return nil, fmt.Errorf("%s size is not divisible by %d rows", path, rows)
	}
	return readDosageMatrix(path, rows, int(info.Size()/int64(rows)))
}

func addIntercept(covariates *mat.Dense) *mat.Dense {
	rows, columns := covariates.Dims()
	x := mat.NewDense(rows, columns+1, nil)
	for row := 0; row < rows; row++ {
		x.Set(row, 0, 1)
		for column := 0; column < columns; column++ {
			x.Set(row, column+1, covariates.At(row, column))
		}
	}
	return x
}
