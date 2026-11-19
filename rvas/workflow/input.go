package workflow

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/hhcho/sfgwas/rvas/protocol"
	"gonum.org/v1/gonum/mat"
)

const (
	auxiliaryPartyID = iota
	cohortAPartyID
	cohortBPartyID
)

type ChromosomeInput struct {
	Chromosome               int
	GenotypeDirectory        string
	PrivateGenotypeDirectory string
	Genes                    []protocol.Gene
	CryptoParams             protocol.CryptoParams
	// Site-cohort ALT counts per block column.
	CohortCounts        [][]float64
	PrivateCohortCounts [][]float64
}

type PartyInput struct {
	X           *mat.Dense
	Y           *mat.Dense
	SampleCount int
	CohortSize  int
	Chromosomes []ChromosomeInput
}

// readCountLines reads one line of space-separated counts per gene.
func readCountLines(path string, geneCount int) ([][]float64, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
	if len(lines) != geneCount {
		return nil, fmt.Errorf("%s has %d lines, expected %d genes", path, len(lines), geneCount)
	}
	counts := make([][]float64, geneCount)
	for index, line := range lines {
		for _, field := range strings.Fields(line) {
			value, err := strconv.ParseFloat(field, 64)
			if err != nil {
				return nil, fmt.Errorf("%s:%d: %w", path, index+1, err)
			}
			counts[index] = append(counts[index], value)
		}
	}
	return counts, nil
}

// readCohortCounts reads the files prepare writes next to the genotype directory.
func readCohortCounts(input *PartyInput, partyID int) error {
	for index := range input.Chromosomes {
		chromosome := &input.Chromosomes[index]
		siteDirectory := filepath.Dir(chromosome.GenotypeDirectory)
		if index == 0 {
			data, err := os.ReadFile(filepath.Join(siteDirectory, "cohort_size.txt"))
			if err != nil {
				return err
			}
			if input.CohortSize, err = strconv.Atoi(strings.TrimSpace(string(data))); err != nil {
				return err
			}
		}
		counts, err := readCountLines(filepath.Join(siteDirectory, "allele_counts.txt"), len(chromosome.Genes))
		if err != nil {
			return err
		}
		for geneIndex, gene := range chromosome.Genes {
			if len(counts[geneIndex]) != gene.VariantCount {
				return fmt.Errorf("chromosome %d gene %s: %d allele counts for %d variants",
					chromosome.Chromosome, gene.GeneID, len(counts[geneIndex]), gene.VariantCount)
			}
		}
		chromosome.CohortCounts = counts
		if partyID == cohortBPartyID {
			chromosome.PrivateCohortCounts, err = readCountLines(
				filepath.Join(siteDirectory, "private_allele_counts.txt"), len(chromosome.Genes))
			if err != nil {
				return err
			}
		}
	}
	return nil
}

func readTextMatrix(path string, columns int) (*mat.Dense, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}

	var values []float64
	for lineIndex, line := range strings.Split(string(data), "\n") {
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		if len(fields) != columns {
			return nil, fmt.Errorf("%s:%d: expected %d columns, got %d", path, lineIndex+1, columns, len(fields))
		}
		for _, field := range fields {
			value, err := strconv.ParseFloat(field, 64)
			if err != nil {
				return nil, err
			}
			values = append(values, value)
		}
	}
	return mat.NewDense(len(values)/columns, columns, values), nil
}

func readRows(
	covariatePath, phenotypePath string,
	covariateCount, phenotypeCount int,
) (*mat.Dense, *mat.Dense, error) {
	covariates, err := readTextMatrix(covariatePath, covariateCount)
	if err != nil {
		return nil, nil, err
	}
	y, err := readTextMatrix(phenotypePath, phenotypeCount)
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

func readGenes(genePath, countPath string) ([]protocol.Gene, error) {
	geneData, err := os.ReadFile(genePath)
	if err != nil {
		return nil, err
	}
	countData, err := os.ReadFile(countPath)
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

// genotypeBlockRows bounds how many samples' genotypes are decoded to float64 at once.
const genotypeBlockRows = 4096

// readGeneBatch returns each batch gene's genotype products with X, Y and itself,
// reading one gene at a time and decoding its rows in blocks. keepPrivate returns
// which of a gene's private variants to keep; dropped ones become zero columns.
func readGeneBatch(
	genotypeDirectory, privateGenotypeDirectory string,
	genes []protocol.Gene,
	batch protocol.GeneBatch,
	x, y *mat.Dense,
	keepPrivate func(geneIndex, columns int) ([]bool, error),
) ([]protocol.Genotypes, error) {
	if genotypeDirectory == "" {
		return nil, nil
	}
	rows, xColumns := x.Dims()
	_, yColumns := y.Dims()
	genotypes := make([]protocol.Genotypes, len(genes))

	for _, geneIndex := range batch.GeneIndices {
		block := fmt.Sprintf("block.%d.bin", geneIndex)
		p := genes[geneIndex].VariantCount
		var public, private []byte
		var err error
		if p > 0 {
			if public, err = os.ReadFile(filepath.Join(genotypeDirectory, block)); err != nil {
				return nil, err
			}
		}
		if privateGenotypeDirectory != "" {
			if private, err = os.ReadFile(filepath.Join(privateGenotypeDirectory, block)); err != nil {
				return nil, err
			}
		}
		if len(public) != rows*p || len(private)%rows != 0 {
			return nil, fmt.Errorf("%s: genotype file sizes do not match %d samples", block, rows)
		}
		q := len(private) / rows
		var keep []bool
		if q > 0 {
			if keep, err = keepPrivate(geneIndex, q); err != nil {
				return nil, err
			}
		}
		m := p + q
		if m == 0 {
			continue
		}

		gtg, gtx, gty := mat.NewDense(m, m, nil), mat.NewDense(m, xColumns, nil), mat.NewDense(m, yColumns, nil)
		buffer := mat.NewDense(min(genotypeBlockRows, rows), m, nil)
		var product mat.Dense
		accumulate := func(sum *mat.Dense, left, right mat.Matrix) {
			product.Reset()
			product.Mul(left.T(), right)
			sum.Add(sum, &product)
		}
		for start := 0; start < rows; start += genotypeBlockRows {
			end := min(start+genotypeBlockRows, rows)
			g := buffer.Slice(0, end-start, 0, m).(*mat.Dense)
			for row := start; row < end; row++ {
				values := g.RawRowView(row - start)
				for j := 0; j < p; j++ {
					values[j] = float64(int8(public[row*p+j]))
				}
				for j := 0; j < q; j++ {
					if keep[j] {
						values[p+j] = float64(int8(private[row*q+j]))
					}
				}
			}
			accumulate(gtg, g, g)
			accumulate(gtx, g, x.Slice(start, end, 0, xColumns))
			accumulate(gty, g, y.Slice(start, end, 0, yColumns))
		}

		sub := func(a *mat.Dense, i, k, j, l int) *mat.Dense { return a.Slice(i, k, j, l).(*mat.Dense) }
		gene := &genotypes[geneIndex]
		if p > 0 {
			gene.GpX, gene.GpY, gene.GpGp = sub(gtx, 0, p, 0, xColumns), sub(gty, 0, p, 0, yColumns), sub(gtg, 0, p, 0, p)
		}
		if q > 0 {
			gene.GvX, gene.GvY, gene.GvGv = sub(gtx, p, m, 0, xColumns), sub(gty, p, m, 0, yColumns), sub(gtg, p, m, p, m)
		}
		if p > 0 && q > 0 {
			gene.GpGv = sub(gtg, 0, p, p, m)
		}
	}
	return genotypes, nil
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
	if partyID < auxiliaryPartyID || partyID > cohortBPartyID {
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

	firstChromosome := config.Chromosomes[0]
	var err error
	input.X, input.Y, err = readRows(
		expandInputPath(config.CovariateFile, ancestry, firstChromosome),
		expandInputPath(config.PhenotypeFile, ancestry, firstChromosome),
		len(config.CovariateColumns),
		len(config.PhenotypeColumns),
	)
	if err != nil {
		return nil, err
	}
	input.SampleCount, _ = input.X.Dims()

	for index, chromosome := range config.Chromosomes {
		genes, err := readGenes(
			expandInputPath(config.GenesFile, ancestry, chromosome),
			expandInputPath(config.VariantCountsFile, ancestry, chromosome),
		)
		if err != nil {
			return nil, err
		}

		input.Chromosomes[index] = ChromosomeInput{
			Chromosome: chromosome,
			GenotypeDirectory: expandInputPath(
				config.GenotypeDirectory,
				ancestry,
				chromosome,
			),
			PrivateGenotypeDirectory: expandInputPath(
				config.PrivateGenotypeDirectory,
				ancestry,
				chromosome,
			),
			Genes: genes,
		}
	}

	if err := readCohortCounts(input, partyID); err != nil {
		return nil, fmt.Errorf("read cohort allele counts (rerun prepare): %w", err)
	}
	return input, nil
}

func expandInputPath(path, ancestry string, chromosome int) string {
	path = strings.ReplaceAll(path, "{ancestry}", ancestry)
	return strings.ReplaceAll(path, "{chromosome}", strconv.Itoa(chromosome))
}
