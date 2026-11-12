package jointstat

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/hhcho/sfgwas/rvas/workflow"
)

// All inputs in this file are LOCAL. Private variant keys never enter the wire layout.
type geneInput struct {
	ID      string
	Public  []string
	Private []string
}
type variant struct {
	Key                  string
	AC, Called, CohortAC int
}
type prepared struct {
	Samples, CohortSize       int
	Genes                     []geneInput
	Public, Private           []variant
	PublicGenes, PrivateGenes [][]int
}

func expand(path, ancestry string, chromosome int) string {
	return strings.NewReplacer("{ancestry}", ancestry, "{chromosome}", strconv.Itoa(chromosome)).Replace(path)
}
func sitePath(config *workflow.Config, ancestry string, chromosome int) string {
	return filepath.Dir(expand(config.GenotypeDirectory, ancestry, chromosome))
}
func readGenes(site string, party int) ([]geneInput, error) {
	parent := filepath.Dir(site)
	raw, err := os.ReadFile(filepath.Join(parent, "genes.txt"))
	if err != nil {
		return nil, err
	}
	ids := strings.Fields(string(raw))
	raw, err = os.ReadFile(filepath.Join(parent, "block_sizes.txt"))
	if err != nil {
		return nil, err
	}
	sizes := strings.Fields(string(raw))
	if len(ids) == 0 || len(ids) != len(sizes) {
		return nil, fmt.Errorf("gene/block count mismatch at %s", parent)
	}
	genes := make([]geneInput, len(ids))
	for i, id := range ids {
		genes[i].ID = id
	}
	file, err := os.Open(filepath.Join(parent, "public_variants.tsv"))
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
	geneCol, keyCol := -1, -1
	for i, name := range header {
		switch name {
		case "gene_id":
			geneCol = i
		case "variant_key":
			keyCol = i
		}
	}
	if geneCol < 0 || keyCol < 0 {
		return nil, fmt.Errorf("public_variants.tsv requires gene_id and variant_key")
	}
	g := -1
	for {
		row, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if g < 0 || row[geneCol] != genes[g].ID {
			g++
			if g >= len(genes) || row[geneCol] != genes[g].ID {
				return nil, fmt.Errorf("public_variants.tsv gene order differs from genes.txt")
			}
		}
		if row[keyCol] != "" {
			genes[g].Public = append(genes[g].Public, row[keyCol])
		}
	}
	if g+1 != len(genes) {
		return nil, fmt.Errorf("public_variants.tsv gene count differs from genes.txt")
	}
	for i, size := range sizes {
		n, err := strconv.Atoi(size)
		if err != nil || n != len(genes[i].Public) {
			return nil, fmt.Errorf("gene %s: public variant count differs from block_sizes.txt", ids[i])
		}
	}
	if party == 2 {
		raw, err := os.ReadFile(filepath.Join(site, "private_variants.txt"))
		if err != nil {
			return nil, fmt.Errorf("read private variant keys (rerun prepare into a new run_dir for older inputs): %w", err)
		}
		lines := strings.Split(strings.TrimSuffix(string(raw), "\n"), "\n")
		if len(lines) != len(genes) {
			return nil, fmt.Errorf("private_variants.txt gene count differs from genes.txt")
		}
		for i, line := range lines {
			genes[i].Private = strings.Fields(line)
		}
	}
	return genes, nil
}
func countLines(path string, genes []geneInput, private bool) ([][]int, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	lines := strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
	if len(lines) != len(genes) {
		return nil, fmt.Errorf("%s: gene count mismatch", path)
	}
	result := make([][]int, len(genes))
	for i, line := range lines {
		keys := genes[i].Public
		if private {
			keys = genes[i].Private
		}
		fields := strings.Fields(line)
		if len(fields) != len(keys) {
			return nil, fmt.Errorf("%s: gene %d variant count mismatch", path, i)
		}
		for _, field := range fields {
			n, e := strconv.Atoi(field)
			if e != nil || n < 0 {
				return nil, fmt.Errorf("invalid AC in %s", path)
			}
			result[i] = append(result[i], n)
		}
	}
	return result, nil
}

type block struct {
	Dosage, Calls       []byte
	Columns, CallOffset int
}

func readBlock(site string, index, samples, columns int, private bool) (block, error) {
	return readBlockRows(site, index, samples, columns, private, 0, samples)
}

// Read just the current sample batch, including non-byte-aligned callability bits.
func readBlockRows(site string, index, samples, columns int, private bool, start, end int) (block, error) {
	folder, prefix := "geno", "public"
	if private {
		folder, prefix = "private", "private"
	}
	b := block{Columns: columns, CallOffset: start * columns % 8}
	read := func(path string, size, offset, length int) ([]byte, error) {
		f, err := os.Open(path)
		if err != nil {
			return nil, err
		}
		defer f.Close()
		info, err := f.Stat()
		if err != nil {
			return nil, err
		}
		if info.Size() != int64(size) {
			return nil, fmt.Errorf("gene %d: dosage/callability dimensions differ", index)
		}
		data := make([]byte, length)
		_, err = f.ReadAt(data, int64(offset))
		return data, err
	}
	var err error
	b.Dosage, err = read(filepath.Join(site, folder, fmt.Sprintf("block.%d.bin", index)), samples*columns, start*columns, (end-start)*columns)
	if err != nil {
		return b, err
	}
	b.Calls, err = read(filepath.Join(site, "joint_stat", fmt.Sprintf("%s.%d.calls", prefix, index)), (samples*columns+7)/8, start*columns/8, (end*columns+7)/8-start*columns/8)
	if err != nil {
		return b, err
	}
	for j, d := range b.Dosage {
		if d > 2 || (!b.called(j) && d != 0) {
			return b, fmt.Errorf("gene %d: invalid dosage or missing-call sidecar", index)
		}
	}
	return b, nil
}
func (b block) called(index int) bool {
	index += b.CallOffset
	return b.Calls[index/8]&(1<<uint(index%8)) != 0
}

// prepareInput uses the already selected analysis rows. It never selects new samples.
func prepareInput(site string, party int) (*prepared, error) {
	genes, err := readGenes(site, party)
	if err != nil {
		return nil, err
	}
	file, err := os.Open(filepath.Join(site, "pheno.txt"))
	if err != nil {
		return nil, err
	}
	defer file.Close()
	samples := 0
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) != "" {
			samples++
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	if samples == 0 {
		return nil, fmt.Errorf("empty phenotype input at %s", site)
	}
	raw, err := os.ReadFile(filepath.Join(site, "cohort_size.txt"))
	if err != nil {
		return nil, err
	}
	cohort, err := strconv.Atoi(strings.TrimSpace(string(raw)))
	if err != nil || cohort < samples {
		return nil, fmt.Errorf("invalid cohort size at %s", site)
	}
	counts, err := countLines(filepath.Join(site, "allele_counts.txt"), genes, false)
	if err != nil {
		return nil, err
	}
	privateCounts := make([][]int, len(genes))
	if party == 2 {
		privateCounts, err = countLines(filepath.Join(site, "private_allele_counts.txt"), genes, true)
		if err != nil {
			return nil, err
		}
	}
	p := &prepared{Samples: samples, CohortSize: cohort, Genes: genes, PublicGenes: make([][]int, len(genes)), PrivateGenes: make([][]int, len(genes))}
	indices := []map[string]int{{}, {}}
	seenGenes := map[string]bool{}
	for g, gene := range genes {
		if gene.ID == "" || seenGenes[gene.ID] {
			return nil, fmt.Errorf("empty/duplicate gene ID")
		}
		seenGenes[gene.ID] = true
		if party == 1 && len(gene.Private) != 0 {
			return nil, fmt.Errorf("party 1 cannot have private columns")
		}
		for side, keys := range [][]string{gene.Public, gene.Private} {
			if side == 1 && party == 1 {
				continue
			}
			b, e := readBlock(site, g, samples, len(keys), side == 1)
			if e != nil {
				return nil, e
			}
			localCounts := counts[g]
			target := &p.Public
			mapping := &p.PublicGenes[g]
			if side == 1 {
				localCounts = privateCounts[g]
				target = &p.Private
				mapping = &p.PrivateGenes[g]
			}
			seen := map[string]bool{}
			for col, key := range keys {
				// Python accepts both chr21 and 21. Normalize before deduplicating.
				key = strings.TrimPrefix(key, "chr")
				if seen[key] {
					return nil, fmt.Errorf("duplicate normalized key %s", key)
				}
				keys[col] = key
				seen[key] = true
				v := variant{Key: key, CohortAC: localCounts[col]}
				if v.CohortAC > 2*cohort {
					return nil, fmt.Errorf("AC exceeds cohort size")
				}
				for row := 0; row < samples; row++ {
					j := row*len(keys) + col
					v.AC += int(b.Dosage[j])
					if b.called(j) {
						v.Called++
					}
				}
				if v.AC > v.CohortAC {
					return nil, fmt.Errorf("analysis AC exceeds ALL-cohort AC for %s", key)
				}
				index, ok := indices[side][key]
				if ok {
					if (*target)[index] != v {
						return nil, fmt.Errorf("inconsistent duplicate variant %s", key)
					}
				} else {
					index = len(*target)
					indices[side][key] = index
					*target = append(*target, v)
				}
				*mapping = append(*mapping, index)
			}
		}
	}
	for key := range indices[1] {
		if _, ok := indices[0][key]; ok {
			return nil, fmt.Errorf("variant appears in both public and private inputs: %s", key)
		}
	}
	return p, nil
}
