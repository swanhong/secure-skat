package gwas

import (
	"fmt"
	"os"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"gonum.org/v1/gonum/mat"
)

// skatBlockNumSnps returns block b's SNP count, syncing pid 0 over the network.
func (ast *AssocTest) skatBlockNumSnps(block int) int {
	mpcObj := ast.general.mpcObj[0]
	pid := mpcObj.GetPid()
	hubPid := mpcObj.GetHubPid()
	if pid == 0 {
		return mpcObj.Network.ReceiveInt(hubPid)
	}

	var nsnpsBlock int
	blockSize := ast.general.genoBlockSizes[block]
	shift := uint64(0)
	for i := 0; i < block; i++ {
		shift += uint64(ast.general.genoBlockSizes[i])
	}
	if ast.general.IsPgen() {
		if ast.general.gwasParams.snpFilt == nil {
			nsnpsBlock = blockSize
		} else {
			nsnpsBlock = SumBool(ast.general.gwasParams.snpFilt[shift : shift+uint64(blockSize)])
		}
	} else {
		nsnpsBlock = int(ast.general.genoBlocks[block].NumColsToKeep())
	}

	if pid == hubPid {
		mpcObj.Network.SendInt(nsnpsBlock, 0)
	}
	return nsnpsBlock
}

// openBlockGenoStream returns this party's genotype stream for block b (the open "blocks"
// stream, or a pgen block materialized to a temp file). nil for an empty block.
func (ast *AssocTest) openBlockGenoStream(b int) *GenoFileStream {
	if !ast.general.IsPgen() {
		if b < len(ast.general.genoBlocks) {
			return ast.general.genoBlocks[b]
		}
		return nil
	}

	gp := ast.general.gwasParams
	blockSize := ast.general.genoBlockSizes[b]
	shift := 0
	for i := 0; i < b; i++ {
		shift += ast.general.genoBlockSizes[i]
	}
	var snpFilt []bool
	if gp.snpFilt == nil {
		snpFilt = OnesBool(blockSize)
	} else {
		snpFilt = gp.snpFilt[shift : shift+blockSize]
	}
	nsnps := SumBool(snpFilt)
	if nsnps == 0 {
		return nil
	}
	numInd := ast.skatNumInds()[ast.general.mpcObj[0].GetPid()]
	pgenFile := fmt.Sprintf(ast.general.config.GenoFilePrefix, b+1)
	tmp := ast.general.CachePath(fmt.Sprintf("lowrank_pgen_gfs.%d.tmp", b))
	FilterMatrixFilePgen(pgenFile, numInd, nsnps, ast.general.config.SampleKeepFile,
		ast.general.config.SnpIdsFile, shift, snpFilt, tmp, false)
	return NewGenoFileStream(tmp, uint64(numInd), uint64(nsnps), true)
}

// denseFromStream reads a genotype stream into a dense matrix (samples × variants), missing
// (negative) dosages → 0. nil/empty streams yield nil.
func denseFromStream(gfs *GenoFileStream) *mat.Dense {
	if gfs == nil {
		return nil
	}
	gfs.Reset()
	n, m := int(gfs.NumRowsToKeep()), int(gfs.NumColsToKeep())
	if n == 0 || m == 0 {
		for gfs.NextRow() != nil {
		}
		return nil
	}

	data := make([]float64, n*m)
	for i := 0; i < n; i++ {
		row := gfs.NextRow()
		if len(row) != m {
			panic(fmt.Sprintf("denseFromStream: row %d has %d columns, want %d", i, len(row), m))
		}
		for j, v := range row {
			if v > 0 {
				data[i*m+j] = float64(v)
			}
		}
	}
	if row := gfs.NextRow(); row != nil {
		panic(fmt.Sprintf("denseFromStream: got more than %d rows", n))
	}
	return mat.NewDense(n, m, data)
}

// orientGenotypeLocal recodes each locally major-coded column in place so that dosage always counts
// that cohort's local minor allele. The strict sum>n test is p_i>1/2; ties are left unchanged. All
// downstream contractions must use this same matrix, because the recode affects score, Gram, Burden,
// and public-private cross terms together.
func orientGenotypeLocal(G *mat.Dense) *mat.Dense {
	if G == nil {
		return nil
	}
	n, m := G.Dims()
	for j := 0; j < m; j++ {
		sum := 0.0
		for i := 0; i < n; i++ {
			sum += G.At(i, j)
		}
		if sum <= float64(n) {
			continue
		}
		for i := 0; i < n; i++ {
			G.Set(i, j, 2.0-G.At(i, j))
		}
	}
	return G
}

func orientedGenotypeLocalCopy(G *mat.Dense) *mat.Dense {
	if G == nil {
		return nil
	}
	return orientGenotypeLocal(mat.DenseCopyOf(G))
}

// loadDenseBlocks reads nGenes per-gene int8 genotype block files into dense matrices for the
// federated-private private side. path b = fmt.Sprintf(prefix, b), row-major n×m_b with m_b
// inferred from file size / n (int8 = 1 byte/cell).
func loadDenseBlocks(prefix string, nGenes, n int) ([]*mat.Dense, error) {
	if n <= 0 {
		return nil, fmt.Errorf("loadDenseBlocks: n=%d must be positive", n)
	}
	blocks := make([]*mat.Dense, nGenes)
	for b := 0; b < nGenes; b++ {
		path := fmt.Sprintf(prefix, b)
		fi, err := os.Stat(path)
		if err != nil {
			return nil, err
		}
		if fi.Size()%int64(n) != 0 {
			return nil, fmt.Errorf("loadDenseBlocks: %s size %d not divisible by n=%d", path, fi.Size(), n)
		}
		m := int(fi.Size()) / n
		if m == 0 {
			blocks[b] = nil // gene with no private variants; privateQRaw(nil) -> Enc(0)
			continue
		}
		blocks[b] = denseFromStream(NewGenoFileStream(path, uint64(n), uint64(m), false))
	}
	return blocks, nil
}

// readGenoBlockLocal streams this party's genotype block b into a dense matrix
// (samples × variants), missing (negative) dosages → 0. Empty on pid 0 or an empty block.
func (ast *AssocTest) readGenoBlockLocal(b int) *mat.Dense {
	return denseFromStream(ast.openBlockGenoStream(b))
}

// blockStat computes one block's raw statistics as 1-elem RVecs: qRawSS = Σw²s² and
// bLinSS = Σw·s (burden linear term, squared by the caller).
func (ast *AssocTest) blockStat(b, nsnps int, null skatNull, X *mat.Dense, y0 []float64, gl *geneLocal) (qRawSS, bLinSS, weightSS mpc_core.RVec) {
	mpcObj := ast.general.mpcObj[0]
	cps := ast.general.cps
	pid := mpcObj.GetPid()

	var GtX *mat.Dense
	var Gty0 []float64
	dosage := make([]float64, nsnps)
	if pid > 0 {
		lc := ast.localFor(b, nsnps, X, y0, gl).LocalContraction
		GtX, Gty0, dosage = lc.GtX, lc.Gty0, lc.DosageSum
	}

	// Score in secret shares — the Gᵀy₀ − GᵀX·β̂ cancellation is exact in fixed-point (β̂ from the
	// Cholesky null solve is accurate). All parties (incl. pid 0, with zero shares) take part.
	sCVec := mpcObj.SSToCVec(cps, ast.scoreSS(GtX, Gty0, null.betaSS, nsnps))
	weightEnc, weightSS := ast.blindWeightCKKS(dosage, nsnps)

	var qRaw, bLin crypto.CipherVector
	if pid > 0 && len(sCVec) > 0 {
		qRaw, bLin = ast.scoreCalculation(sCVec, weightEnc)
	}
	return ast.scalarCiphertextToShares(qRaw), ast.scalarCiphertextToShares(bLin), weightSS
}
