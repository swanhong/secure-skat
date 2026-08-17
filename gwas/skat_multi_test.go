package gwas

import (
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"gonum.org/v1/gonum/mat"
)

func TestLocalMultiNullEquationsMatchIndependentPhenotypes(t *testing.T) {
	cov := mat.NewDense(5, 2, []float64{
		0.2, -0.3,
		1.1, 0.4,
		-0.7, 0.8,
		0.5, -1.2,
		1.4, 0.6,
	})
	pheno := mat.NewDense(5, 3, []float64{
		1.0, 3.1, -0.2,
		2.2, 2.4, 0.7,
		-0.5, 1.3, 1.1,
		0.8, -0.4, 2.0,
		1.7, 0.9, -1.3,
	})

	multi := localMultiNullEquations(cov, pheno, 0, 3)
	for phenotype := 0; phenotype < 3; phenotype++ {
		column := mat.NewDense(5, 1, mat.Col(nil, phenotype, pheno))
		single := localNullEquations(cov, column, 0)
		for row := 0; row < 3; row++ {
			if !approxEqual(multi.Xty0.At(row, phenotype), single.Xty0[row], 1e-12) {
				t.Fatalf("Xty[%d,%d]=%v want %v", row, phenotype, multi.Xty0.At(row, phenotype), single.Xty0[row])
			}
		}
		if !approxEqual(multi.Yty0[phenotype], single.Y0ty0, 1e-12) {
			t.Fatalf("yty[%d]=%v want %v", phenotype, multi.Yty0[phenotype], single.Y0ty0)
		}
	}
}

func TestFlattenGenePhenotypesUsesGeneMajorOrder(t *testing.T) {
	var rtype mpc_core.LElem128
	values := []mpc_core.RVec{
		mpc_core.IntToRVec(rtype, []int{10, 20}),
		mpc_core.IntToRVec(rtype, []int{11, 21}),
		mpc_core.IntToRVec(rtype, []int{12, 22}),
	}
	flattened := flattenGenePhenotypes(rtype, values)
	want := []uint64{10, 11, 12, 20, 21, 22}
	for i := range want {
		if flattened[i].Uint64() != want[i] {
			t.Fatalf("flattened[%d]=%d want %d", i, flattened[i].Uint64(), want[i])
		}
	}

	repeated := repeatGenes(rtype, mpc_core.IntToRVec(rtype, []int{7, 8}), 3)
	want = []uint64{7, 7, 7, 8, 8, 8}
	for i := range want {
		if repeated[i].Uint64() != want[i] {
			t.Fatalf("repeated[%d]=%d want %d", i, repeated[i].Uint64(), want[i])
		}
	}
}

func TestPackedMultiPhenotypeEqualsIndependentQ1(t *testing.T) {
	t.Setenv("SFGWAS_MODE", "skat_fed")
	const genes, variants, privateVariants = 2, 3, 2
	nByParty := []int{0, 12, 11}
	q := 3
	if value := os.Getenv("SFGWAS_TEST_Q"); value != "" {
		var err error
		q, err = strconv.Atoi(value)
		if err != nil || q < 1 {
			t.Fatalf("invalid SFGWAS_TEST_Q=%q", value)
		}
	}
	probes := 0
	if value := os.Getenv("SFGWAS_TEST_PROBES"); value != "" {
		var err error
		probes, err = strconv.Atoi(value)
		if err != nil {
			t.Fatal(err)
		}
	}

	writeMatrix := func(path string, rows, columns int, value func(int, int) float64) {
		var text strings.Builder
		for i := 0; i < rows; i++ {
			for j := 0; j < columns; j++ {
				if j > 0 {
					text.WriteByte('\t')
				}
				fmt.Fprintf(&text, "%.12g", value(i, j))
			}
			text.WriteByte('\n')
		}
		if err := os.WriteFile(path, []byte(text.String()), 0644); err != nil {
			t.Fatal(err)
		}
	}
	writeBlock := func(path string, rows, columns, party, gene int) {
		data := make([]byte, rows*columns)
		for i := 0; i < rows; i++ {
			for j := 0; j < columns; j++ {
				data[i*columns+j] = byte(int8((i + 2*j + party + gene) % 3))
			}
		}
		if err := os.WriteFile(path, data, 0644); err != nil {
			t.Fatal(err)
		}
	}

	prot := InitProtocolForTestWithConfig(t, func(config *Config, pid int) {
		root := t.TempDir()
		config.CkksParams = crypto.CKKSParamsPN14QP436S45
		config.NumInds = append([]int(nil), nByParty...)
		config.NumSnps = genes * variants
		config.NumCovs = 1
		config.NumPhenos = q
		config.GenoFileFormat = "blocks"
		config.GenoNumBlocks = genes
		config.PrivatePid = 2
		config.SkatPValueProbes = probes
		config.RotKeyPow2Only = true
		config.SkipQC = true
		config.SkipPCA = true
		config.GenoBlockSizeFile = filepath.Join(root, "block_sizes.txt")
		config.GeneIDFile = filepath.Join(root, "genes.txt")
		config.OutDir = filepath.Join(root, "out")
		config.CacheDir = filepath.Join(root, "cache")
		if err := os.WriteFile(config.GenoBlockSizeFile, []byte("3\n3\n"), 0644); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(config.GeneIDFile, []byte("gene0\ngene1\n"), 0644); err != nil {
			t.Fatal(err)
		}
		if pid == 0 {
			return
		}

		n := nByParty[pid]
		prefix := filepath.Join(root, "geno")
		for gene := 0; gene < genes; gene++ {
			writeBlock(fmt.Sprintf("%s.%d.bin", prefix, gene), n, variants, pid, gene)
		}
		config.GenoFilePrefix = prefix
		config.PhenoFile = filepath.Join(root, "pheno.txt")
		config.CovFile = filepath.Join(root, "cov.txt")
		config.SnpPosFile = filepath.Join(root, "pos.txt")
		writeMatrix(config.PhenoFile, n, q, func(i, phenotype int) float64 {
			return 0.2*float64(i+1) + 0.7*float64(phenotype) + float64((i+phenotype*2)%4)
		})
		writeMatrix(config.CovFile, n, 1, func(i, _ int) float64 {
			return float64((i+pid)%7) / 5
		})
		writeMatrix(config.SnpPosFile, genes*variants, 2, func(i, column int) float64 {
			if column == 0 {
				return 22
			}
			return float64(1000 + i)
		})
		if pid == 2 {
			privatePrefix := filepath.Join(root, "priv.%d.bin")
			for gene := 0; gene < genes; gene++ {
				writeBlock(fmt.Sprintf(privatePrefix, gene), n, privateVariants, pid+2, gene)
			}
			config.PrivateOnlyPrefix = privatePrefix
		}
	})
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	pid := prot.mpcObj[0].GetPid()
	var allPheno *mat.Dense
	if pid > 0 {
		allPheno = mat.DenseCopyOf(prot.pheno)
	}
	multiQ, multiBurden, multiSkatP := prot.runFederatedPrivate()
	independentQ := make([][]float64, q)
	independentBurden := make([][]float64, q)
	independentSkatP := make([][]float64, q)
	for phenotype := 0; phenotype < q; phenotype++ {
		prot.config.NumPhenos = 1
		if pid > 0 {
			rows, _ := allPheno.Dims()
			prot.pheno = mat.NewDense(rows, 1, mat.Col(nil, phenotype, allPheno))
		}
		independentQ[phenotype], independentBurden[phenotype], independentSkatP[phenotype] = prot.runFederatedPrivate()
	}
	if pid != 1 {
		return
	}
	for gene := 0; gene < genes; gene++ {
		for phenotype := 0; phenotype < q; phenotype++ {
			index := gene*q + phenotype
			if probes == 0 {
				if !approxEqual(multiQ[index], independentQ[phenotype][gene], 2e-2) {
					t.Fatalf("Q[%d,%d]=%v independent=%v", gene, phenotype, multiQ[index], independentQ[phenotype][gene])
				}
			} else if !approxEqual(multiSkatP[index], independentSkatP[phenotype][gene], 5e-2) {
				t.Fatalf("SKAT p[%d,%d]=%v independent=%v", gene, phenotype, multiSkatP[index], independentSkatP[phenotype][gene])
			}
			if !approxEqual(multiBurden[index], independentBurden[phenotype][gene], 2e-3) {
				t.Fatalf("burden p[%d,%d]=%v independent=%v", gene, phenotype, multiBurden[index], independentBurden[phenotype][gene])
			}
		}
	}
}
