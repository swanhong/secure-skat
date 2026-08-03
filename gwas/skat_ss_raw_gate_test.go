package gwas

import (
	"fmt"
	"math"
	"os"
	"runtime"
	"strconv"
	"testing"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"golang.org/x/sys/unix"
	"gonum.org/v1/gonum/mat"
)

type ssRawGateResult struct {
	q, l             mpc_core.RVec
	variants, chunks int
}

// publicSSRawGate batches the public-list score, weight, Q/N, and L/sqrt(N)
// across gene boundaries. It is an experimental Step 3 backend, not production.
func (ast *AssocTest) publicSSRawGate(null skatNull, X *mat.Dense, y0 []float64, publicSizes []int, chunkSize int) ssRawGateResult {
	if chunkSize <= 0 {
		panic("SS raw chunk size must be positive")
	}

	mpcObj := ast.general.mpcObj[0]
	rtype := mpcObj.GetRType()
	pid := mpcObj.GetPid()
	c := len(null.betaSS)
	N := float64(ast.skatTotalNumInds())
	result := ssRawGateResult{
		q: mpc_core.InitRVec(rtype.Zero(), len(publicSizes)),
		l: mpc_core.InitRVec(rtype.Zero(), len(publicSizes)),
	}

	gtxData := make([]float64, chunkSize*c)
	gtyData := make([]float64, chunkSize)
	dosageData := make([]float64, chunkSize)
	used, outputGene, outputInGene := 0, 0, 0
	advanceEmptyGenes := func() {
		for outputGene < len(publicSizes) && publicSizes[outputGene] == 0 {
			outputGene++
		}
	}
	advanceEmptyGenes()

	flush := func(n int) {
		var gtx *mat.Dense
		var gty, dosage []float64
		if pid > 0 {
			gtx = mat.NewDense(n, c, gtxData[:n*c])
			gty = gtyData[:n]
			dosage = dosageData[:n]
		}

		score := ast.scoreSS(gtx, gty, null.betaSS, n)
		maf := mpc_core.InitRVec(rtype.Zero(), n)
		if pid > 0 {
			for i := range maf {
				maf[i] = rtype.FromFloat64(dosage[i]/(2*N), mpcObj.GetFracBits())
			}
		}
		base := ast.hubVec(1, n)
		base.Sub(maf)
		w2 := ast.ssSquare(base)
		w4 := ast.ssSquare(w2)
		w8 := ast.ssSquare(w4)
		w16 := ast.ssSquare(w8)
		weight := ast.ssMul(w16, w8)
		weight.MulScalar(rtype.FromFloat64(25, 0))

		x := ast.ssPMul(ast.ssMul(weight, score), 1/math.Sqrt(N))
		q := ast.ssSquare(x)
		result.variants += n
		result.chunks++

		for i := 0; i < n; i++ {
			result.q[outputGene] = result.q[outputGene].Add(q[i])
			result.l[outputGene] = result.l[outputGene].Add(x[i])
			outputInGene++
			if outputInGene == publicSizes[outputGene] {
				outputGene++
				outputInGene = 0
				advanceEmptyGenes()
			}
		}
	}

	for gene, m := range publicSizes {
		if m == 0 {
			continue
		}
		var local LocalContraction
		if pid > 0 {
			G := orientGenotypeLocal(ast.readGenoBlockLocal(gene))
			if G == nil {
				panic(fmt.Sprintf("gene %d genotype block is empty", gene))
			}
			_, got := G.Dims()
			if got != m {
				panic(fmt.Sprintf("gene %d has %d public variants, want %d", gene, got, m))
			}
			local = localGenotypeContract(G, X, y0)
		}

		for row := 0; row < m; {
			take := min(chunkSize-used, m-row)
			if pid > 0 {
				copy(gtyData[used:used+take], local.Gty0[row:row+take])
				copy(dosageData[used:used+take], local.DosageSum[row:row+take])
				for i := 0; i < take; i++ {
					for j := 0; j < c; j++ {
						gtxData[(used+i)*c+j] = local.GtX.At(row+i, j)
					}
				}
			}
			used += take
			row += take
			if used == chunkSize {
				flush(used)
				used = 0
			}
		}
	}
	if used > 0 {
		flush(used)
	}
	return result
}

func processPeakRSS() uint64 {
	var usage unix.Rusage
	if unix.Getrusage(unix.RUSAGE_SELF, &usage) != nil {
		return 0
	}
	rss := uint64(usage.Maxrss)
	if runtime.GOOS != "darwin" {
		rss *= 1024
	}
	return rss
}

// TestSKATPublicSSRawGate is retained only until the real-N Step 3 measurement is collected.
func TestSKATPublicSSRawGate(t *testing.T) {
	if os.Getenv("SS_RAW_GATE") != "1" {
		t.Skip("set SS_RAW_GATE=1 for the Step 3 measurement")
	}
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	chunkSize := 4096
	if value := os.Getenv("SS_RAW_CHUNK"); value != "" {
		parsed, err := strconv.Atoi(value)
		if err != nil || parsed <= 0 {
			t.Fatalf("invalid SS_RAW_CHUNK %q", value)
		}
		chunkSize = parsed
	}

	ast := prot.initSKAT()
	publicSizes := loadPublicGeneSizes(prot.config.GenoBlockSizeFile, prot.config.GenoNumBlocks)
	null, X, y0 := ast.nullSetup()
	runtime.GC()
	networks := prot.mpcObj.GetNetworks()
	commStart := networks.GetCommunicationStats()
	rssBefore := processPeakRSS()
	started := time.Now()
	result := ast.publicSSRawGate(null, X, y0, publicSizes, chunkSize)
	wall := time.Since(started)
	comm := networks.GetCommunicationStats().Sub(commStart)
	rssPeak := processPeakRSS()
	M := result.variants
	c := len(null.betaSS)
	fmt.Printf("[step3-ss] party=%d genes=%d variants=%d chunk=%d chunks=%d wall_ms=%d beaver_calls=%d beaver_mults=%d trunc_calls=%d trunc_elems=%d sent_bytes=%d received_bytes=%d sent_messages=%d received_messages=%d rss_before=%d rss_peak=%d fixed_integer_bits=%d\n",
		prot.mpcObj[0].GetPid(), len(publicSizes), M, chunkSize, result.chunks, wall.Milliseconds(),
		8*result.chunks, M*(c+7), 9*result.chunks, 9*M,
		comm.SentBytes, comm.ReceivedBytes, comm.SentMessages, comm.ReceivedMessages,
		rssBefore, rssPeak, prot.mpcObj[0].GetDataBits()-prot.mpcObj[0].GetFracBits()-1)
}
