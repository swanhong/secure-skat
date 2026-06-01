package gwas

import (
	"fmt"
	"math/bits"
	"runtime"
	"sync"
	"time"
	"unsafe"

	"go.dedis.ch/onet/v3/log"

	"github.com/hhcho/sfgwas/crypto"

	"math"

	"github.com/hhcho/sfgwas/mpc"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"gonum.org/v1/gonum/mat"
)

// Multiply Q (kp by nsnp) with X (nsnp by nind) with lazy normalization of X
// Compute Q * S * (X - m * 1^T) as (Q * S) * X - ((Q * S) * m) * 1^T
// S: diagonal matrix containing 1/stdev of each SNP, m: column vector containing mean of each SNP
func QXLazyNormStream(cps *crypto.CryptoParams, mpcObj *mpc.MPC, Q crypto.CipherMatrix, Xcachefile string, XMean, XStdInv crypto.CipherVector, numInd int) (out crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		return
	}
	slots := cps.GetSlots()

	start := time.Now()

	// Compute Q * S
	QS := make(crypto.CipherMatrix, len(Q))
	for i := range Q {
		QS[i] = crypto.CMult(cps, Q[i], XStdInv)
	}

	// Compute (Q * S) * X
	out = MatMult4StreamCompute(cps, QS, 5, Xcachefile)

	out = mpcObj.Network.BootstrapMatAll(cps, out) // TODO

	// Compute (Q * S) * m
	QSm := make(crypto.CipherVector, len(Q))
	for i := range Q {
		QSm[i] = crypto.InnerProd(cps, QS[i], XMean) // Already has output value in all slots
	}

	// Compute (Q * S) * X - ((Q * S) * m) * 1^T
	for i := range QS {
		cps.WithEvaluator(func(eval *ckks.Evaluator) error {
			for j := range out[i] {
				if err := eval.Sub(out[i][j], QSm[i], out[i][j]); err != nil {
					return err
				}
			}
			return nil
		})

		// Zero out empty slots at the end
		// TODO: Is there a better way?
		for j := range out[i] {
			var N int
			if j < len(out[i])-1 {
				N = slots
			} else {
				N = ((numInd - 1) % slots) + 1
			}
			out[i][j] = crypto.MaskTrunc(cps, out[i][j], N)
		}
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Matrix multiplication complete,", time.Since(start))

	return
}

// Multiply Q (kp by nind) with X^T (nind by nsnp) with lazy normalization of X
// Compute Q * (X^T - 1 * m^T) * S as ((Q * X^T) - ((Q * 1) * m^T)) * S
// S: diagonal matrix containing 1/stdev of each SNP, m: column vector containing mean of each SNP
// TODO: multiply with S AFTER aggregation across parties, that way bootstrap once for all
func QXtLazyNormStream(cps *crypto.CryptoParams, mpcObj *mpc.MPC, Q crypto.CipherMatrix, XTcachefile string, XMean, XStdInv crypto.CipherVector) (out crypto.CipherMatrix) {
	if mpcObj.GetPid() == 0 {
		return
	}

	start := time.Now()

	// Compute Q * X^T
	out = MatMult4StreamCompute(cps, Q, 5, XTcachefile)
	out = mpcObj.Network.BootstrapMatAll(cps, out) // TODO change to bootstrapping after aggregation

	// Compute (Q * X^T) - ((Q * 1) * m^T)
	for i := range out {
		rowSum := crypto.InnerSumAll(cps, Q[i])

		Q1m := crypto.CMultScalar(cps, XMean, rowSum)

		cps.WithEvaluator(func(eval *ckks.Evaluator) error {
			for j := range out[i] {
				if err := eval.Sub(out[i][j], Q1m[j], out[i][j]); err != nil {
					return err
				}
			}
			return nil
		})
	}

	// Compute ((Q * X^T) - ((Q * 1) * m^T)) * S
	for i := range out {
		out[i] = crypto.CMult(cps, out[i], XStdInv)
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Matrix multiplication complete,", time.Since(start))

	return
}

type matmulPlainInnerFn func(*crypto.CryptoParams, crypto.CipherVector, crypto.PlainMatrix, int) crypto.CipherVector
type matmulInnerFn func(*crypto.CryptoParams, crypto.CipherVector, crypto.CipherMatrix, int) crypto.CipherVector

func DCMatMulAAtB(cryptoParams *crypto.CryptoParams, mpcObj *mpc.MPC, A crypto.CipherMatrix, B crypto.CipherMatrix,
	nrows []int, ncol_out int, innerFn matmulInnerFn) crypto.CipherMatrix {

	slots := cryptoParams.GetSlots()
	pid := mpcObj.GetPid()

	out := crypto.CZeroMat(cryptoParams, ((nrows[pid]-1)/slots)+1, ncol_out)

	cTQloc := make(crypto.CipherVector, ncol_out)
	for c := range A {
		var wg sync.WaitGroup
		for j := range cTQloc {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				innerProd := innerFn(cryptoParams, A[c], B, j)
				cTQloc[j] = crypto.InnerSumAll(cryptoParams, innerProd)
			}(j)
		}
		wg.Wait()

		cTQ := mpcObj.Network.AggregateCVec(cryptoParams, cTQloc)

		for j := 0; j < ncol_out; j++ {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				ccTQ := crypto.CMult(cryptoParams, A[c], crypto.CipherVector{cTQ[j]})
				out[j] = crypto.CAdd(cryptoParams, out[j], ccTQ)
			}(j)
		}
		wg.Wait()
	}

	return out
}

func DCMatMulAAtBPlain(cryptoParams *crypto.CryptoParams, mpcObj *mpc.MPC, A crypto.CipherMatrix, B crypto.PlainMatrix,
	nrows []int, ncol_out int, innerFn matmulPlainInnerFn) crypto.CipherMatrix {

	slots := cryptoParams.GetSlots()
	pid := mpcObj.GetPid()

	// Align levels of A and B (if possible)
	A, _ = crypto.FlattenLevels(cryptoParams, A)

	out := crypto.CZeroMat(cryptoParams, int((nrows[pid]-1)/slots)+1, ncol_out)
	// Initialize out with correct scale for subsequent additions
	for i := range out {
		for j := range out[i] {
			out[i][j].Scale = cryptoParams.Params.DefaultScale()
		}
	}

	var wg sync.WaitGroup
	for c := range A {
		cTQloc := make(crypto.CipherVector, ncol_out)

		for j := range cTQloc {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				innerProd := innerFn(cryptoParams, A[c], B, j)
				cTQloc[j] = crypto.InnerSumAll(cryptoParams, innerProd)
			}(j)
		}
		wg.Wait()

		cTQ := mpcObj.Network.AggregateCVec(cryptoParams, cTQloc)
		// Align the aggregated cTQ to the same level as A[c]
		for col := range cTQ {
			if cTQ[col] == nil {
				continue
			}
			levelA := A[c][0].Level()
			if cTQ[col].Level() > levelA {
				cryptoParams.WithEvaluator(func(eval *ckks.Evaluator) error {
					eval.DropLevel(cTQ[col], cTQ[col].Level()-levelA)
					return nil
				})
			}
		}

		for j := 0; j < ncol_out; j++ {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				ctq := cTQ[j]
				// Drop A[c]'s level to match ctq if needed (per ciphertext)
				ac := make(crypto.CipherVector, len(A[c]))
				for k, ct := range A[c] {
					if ct == nil {
						continue
					}
					if ct.Level() > ctq.Level() {
						ctCopy := ct.CopyNew()
						cryptoParams.WithEvaluator(func(eval *ckks.Evaluator) error {
							eval.DropLevel(ctCopy, ctCopy.Level()-ctq.Level())
							return nil
						})
						ac[k] = ctCopy
					} else {
						ac[k] = ct
					}
				}

				ccTQ := crypto.CMult(cryptoParams, ac, crypto.CipherVector{ctq})
				out[j] = crypto.CAdd(cryptoParams, out[j], ccTQ)
			}(j)
		}
		wg.Wait()
	}

	return out
}

func DCMatMulAAtBPlainWithIntmd(cryptoParams *crypto.CryptoParams, mpcObj *mpc.MPC, A crypto.CipherMatrix, B crypto.PlainMatrix,
	nrows []int, ncol_out int, innerFn matmulPlainInnerFn) (crypto.CipherMatrix, crypto.CipherMatrix, error) {

	slots := cryptoParams.GetSlots()
	pid := mpcObj.GetPid()

	// Align levels of A and B (if possible)
	A, _ = crypto.FlattenLevels(cryptoParams, A)

	out := crypto.CZeroMat(cryptoParams, int((nrows[pid]-1)/slots)+1, ncol_out)
	// Initialize out with correct scale for subsequent additions
	for i := range out {
		for j := range out[i] {
			out[i][j].Scale = cryptoParams.Params.DefaultScale()
		}
	}

	var QTY crypto.CipherMatrix
	if pid > 0 {
		QTY = make(crypto.CipherMatrix, 1)
		QTY[0] = make(crypto.CipherVector, len(A))
	} else {
		QTY = make(crypto.CipherMatrix, 1)
		QTY[0] = nil
	}

	var wg sync.WaitGroup
	for c := range A {
		cTQloc := make(crypto.CipherVector, ncol_out)

		for j := range cTQloc {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				innerProd := innerFn(cryptoParams, A[c], B, j)
				cTQloc[j] = crypto.InnerSumAll(cryptoParams, innerProd)
			}(j)
		}
		wg.Wait()

		cTQ := mpcObj.Network.AggregateCVec(cryptoParams, cTQloc)

		if pid > 0 {
			QTY[0][c] = cTQ[0]
		}

		// Align the aggregated cTQ to the same level as A[c]
		for col := range cTQ {
			if cTQ[col] == nil {
				continue
			}
			levelA := A[c][0].Level()
			if cTQ[col].Level() > levelA {
				cryptoParams.WithEvaluator(func(eval *ckks.Evaluator) error {
					eval.DropLevel(cTQ[col], cTQ[col].Level()-levelA)
					return nil
				})
			}
		}

		for j := 0; j < ncol_out; j++ {
			wg.Add(1)
			go func(j int) {
				defer wg.Done()
				ctq := cTQ[j]
				// Drop A[c]'s level to match ctq if needed (per ciphertext)
				ac := make(crypto.CipherVector, len(A[c]))
				for k, ct := range A[c] {
					if ct == nil {
						continue
					}
					if ct.Level() > ctq.Level() {
						ctCopy := ct.CopyNew()
						cryptoParams.WithEvaluator(func(eval *ckks.Evaluator) error {
							eval.DropLevel(ctCopy, ctCopy.Level()-ctq.Level())
							return nil
						})
						ac[k] = ctCopy
					} else {
						ac[k] = ct
					}
				}

				ccTQ := crypto.CMult(cryptoParams, ac, crypto.CipherVector{ctq})
				out[j] = crypto.CAdd(cryptoParams, out[j], ccTQ)
			}(j)
		}
		wg.Wait()
	}

	return out, QTY, nil
}

type uint128 struct {
	hi uint64
	lo uint64
}

type CipherAccV2 struct {
	acc0 [][]uint128
	acc1 [][]uint128
}

type CipherVectorAccV2 struct {
	val []CipherAccV2
}

func NewCipherVectorAccV2(cryptoParams *crypto.CryptoParams, n int, level int) CipherVectorAccV2 {
	N := cryptoParams.Params.N()
	var out CipherVectorAccV2
	out.val = make([]CipherAccV2, n)
	for i := range out.val {
		out.val[i].acc0 = make([][]uint128, level)
		out.val[i].acc1 = make([][]uint128, level)
		for l := 0; l < level; l++ {
			out.val[i].acc0[l] = make([]uint128, N)
			out.val[i].acc1[l] = make([]uint128, N)
		}
	}

	return out
}

func MulCoeffsAndAdd128(a, b []uint64, c []uint128) {

	var hi, lo, carry uint64

	for j := 0; j < len(a); j = j + 8 {

		x := (*[8]uint64)(unsafe.Pointer(&a[j]))
		y := (*[8]uint64)(unsafe.Pointer(&b[j]))
		z := (*[8]uint128)(unsafe.Pointer(&c[j]))

		hi, lo = bits.Mul64(x[0], y[0])
		z[0].lo, carry = bits.Add64(z[0].lo, lo, 0)
		z[0].hi += hi + carry

		hi, lo = bits.Mul64(x[1], y[1])
		z[1].lo, carry = bits.Add64(z[1].lo, lo, 0)
		z[1].hi += hi + carry

		hi, lo = bits.Mul64(x[2], y[2])
		z[2].lo, carry = bits.Add64(z[2].lo, lo, 0)
		z[2].hi += hi + carry

		hi, lo = bits.Mul64(x[3], y[3])
		z[3].lo, carry = bits.Add64(z[3].lo, lo, 0)
		z[3].hi += hi + carry

		hi, lo = bits.Mul64(x[4], y[4])
		z[4].lo, carry = bits.Add64(z[4].lo, lo, 0)
		z[4].hi += hi + carry

		hi, lo = bits.Mul64(x[5], y[5])
		z[5].lo, carry = bits.Add64(z[5].lo, lo, 0)
		z[5].hi += hi + carry

		hi, lo = bits.Mul64(x[6], y[6])
		z[6].lo, carry = bits.Add64(z[6].lo, lo, 0)
		z[6].hi += hi + carry

		hi, lo = bits.Mul64(x[7], y[7])
		z[7].lo, carry = bits.Add64(z[7].lo, lo, 0)
		z[7].hi += hi + carry
	}
}

func ReduceAndAddUint128(in []uint128, out []uint64, qInv, q uint64) {

	var hhi uint64

	for j := 0; j < len(in); j = j + 8 {

		x := (*[8]uint128)(unsafe.Pointer(&in[j]))
		y := (*[8]uint64)(unsafe.Pointer(&out[j]))

		hhi, _ = bits.Mul64(x[0].lo*qInv, q)
		y[0] += x[0].hi - hhi + q

		hhi, _ = bits.Mul64(x[1].lo*qInv, q)
		y[1] += x[1].hi - hhi + q

		hhi, _ = bits.Mul64(x[2].lo*qInv, q)
		y[2] += x[2].hi - hhi + q

		hhi, _ = bits.Mul64(x[3].lo*qInv, q)
		y[3] += x[3].hi - hhi + q

		hhi, _ = bits.Mul64(x[4].lo*qInv, q)
		y[4] += x[4].hi - hhi + q

		hhi, _ = bits.Mul64(x[5].lo*qInv, q)
		y[5] += x[5].hi - hhi + q

		hhi, _ = bits.Mul64(x[6].lo*qInv, q)
		y[6] += x[6].hi - hhi + q

		hhi, _ = bits.Mul64(x[7].lo*qInv, q)
		y[7] += x[7].hi - hhi + q
	}
}

func ModularReduceV2(cryptoParams *crypto.CryptoParams, cva CipherVectorAccV2, outScale float64) crypto.CipherVector {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Q())
	level := len(cva.val[0].acc0)

	out := make(crypto.CipherVector, len(cva.val))
	for i := range out {
		ct := ckks.NewCiphertext(cryptoParams.Params, 1, level-1)
		ct.Scale = rlwe.NewScale(outScale)
		for l := 0; l < level; l++ {
			mredParams := ringQ.SubRings[l].MRedConstant
			qi := ringQ.SubRings[l].Modulus
			ReduceAndAddUint128(cva.val[i].acc0[l], ct.Value[0].Coeffs[l], mredParams, qi)
			ReduceAndAddUint128(cva.val[i].acc1[l], ct.Value[1].Coeffs[l], mredParams, qi)
		}
		ringQ.AtLevel(level-1).Reduce(ct.Value[0], ct.Value[0])
		ringQ.AtLevel(level-1).Reduce(ct.Value[1], ct.Value[1])
		out[i] = ct
	}
	return out
}

func CPMultAccWithoutMRedV2(X crypto.CipherVector, Y crypto.PlainVector, Acc CipherVectorAccV2) {
	n := len(Acc.val)
	for i := 0; i < n; i++ {
		// Broadcasting
		xi, yi := i, i
		if len(X) == 1 {
			xi = 0
		}
		if len(Y) == 1 {
			yi = 0
		}

		if X[xi] != nil && Y[yi] != nil {
			for l := 0; l < len(Acc.val[i].acc0); l++ {
				MulCoeffsAndAdd128(X[xi].Value[0].Coeffs[l], Y[yi].Value.Coeffs[l], Acc.val[i].acc0[l])
				MulCoeffsAndAdd128(X[xi].Value[1].Coeffs[l], Y[yi].Value.Coeffs[l], Acc.val[i].acc1[l])
			}
		}
	}
}

func ToMontgomeryForm(cryptoParams *crypto.CryptoParams, pt crypto.PlainVector) {
	N := cryptoParams.Params.N()
	ringQ, _ := ring.NewRing(N, cryptoParams.Params.Q())
	for i := range pt {
		if pt[i] != nil {
			MFormLvl(ringQ, pt[i].Level(), pt[i].Value, pt[i].Value)
		}
	}
}

func MFormLvl(r *ring.Ring, level int, p1, p2 ring.Poly) {
	for i := 0; i < level+1; i++ {
		qi := r.SubRings[i].Modulus
		bredParams := r.SubRings[i].BRedConstant
		p1tmp, p2tmp := p1.Coeffs[i], p2.Coeffs[i]
		for j := 0; j < r.N(); j = j + 8 {

			x := (*[8]uint64)(unsafe.Pointer(&p1tmp[j]))
			z := (*[8]uint64)(unsafe.Pointer(&p2tmp[j]))

			z[0] = MForm(x[0], qi, bredParams)
			z[1] = MForm(x[1], qi, bredParams)
			z[2] = MForm(x[2], qi, bredParams)
			z[3] = MForm(x[3], qi, bredParams)
			z[4] = MForm(x[4], qi, bredParams)
			z[5] = MForm(x[5], qi, bredParams)
			z[6] = MForm(x[6], qi, bredParams)
			z[7] = MForm(x[7], qi, bredParams)
		}
	}
}

func MForm(a, q uint64, u [2]uint64) (r uint64) {
	mhi, _ := bits.Mul64(a, u[1])
	r = -(a*u[0] + mhi) * q
	if r >= q {
		r -= q
	}
	return
}

type BlockI8 struct {
	data [][]int8
	r    int
	c    int
}

func NewBlockI8(r, c int) BlockI8 {
	return BlockI8{
		data: make([][]int8, r),
		r:    r,
		c:    c,
	}
}

func (b BlockI8) At(i, j int) float64 {
	return float64(b.data[i][j])
}

func (b BlockI8) Dims() (int, int) {
	return b.r, b.c
}

type BlockF64 mat.Dense

func (b BlockF64) At(i, j int) float64 {
	return b.At(i, j)
}

type Block interface {
	At(int, int) float64
	Dims() (int, int)
}

type BlockVector []Block
type BlockMatrix []BlockVector

func ToBlockMatrix(A *mat.Dense, d int) BlockMatrix {
	r, c := A.Dims()
	br, bc := int((r-1)/d)+1, int((c-1)/d)+1

	out := make(BlockMatrix, br)
	for bi := range out {
		out[bi] = make(BlockVector, bc)
		for bj := range out[bi] {
			i1, i2 := bi*d, Min((bi+1)*d, r)
			j1, j2 := bj*d, Min((bj+1)*d, c)
			out[bi][bj] = A.Slice(i1, i2, j1, j2)
		}
	}

	return out
}

// Return if a diagonal vector exists without extracting elements
func GetDiagBool(X Block, dim int, index int) bool {
	r, c := X.Dims()
	index = Mod(index, dim) // range [0, dim-1]
	return (dim+1-r) <= index || index <= c-1
}

// index 0 is the main diagonal
// max size of Block is dim by dim and index ranges from 0 to dim-1 (mod dim)
// If given diagonal does not overlap with X (matrix might be smaller), returns false
func GetDiag(dst []float64, X Block, dim int, index int) bool {
	r, c := X.Dims()

	index = Mod(index, dim) // range [0, dim-1]

	if (dim+1-r) <= index || index <= c-1 {

		if dst == nil {
			dst = make([]float64, dim)
		} else if len(dst) < c {
			panic("destination array is not large enough")
		}

		i := Mod(-index, dim)
		for j := 0; j < len(dst); j++ {
			if i < r && j < c {
				dst[j] = X.At(i, j)
			} else {
				dst[j] = 0
			}

			i = Mod(i+1, dim)
		}

		return true
	}

	return false
}

func convertToComplex128WithRot(v []float64, nrot int) []complex128 {
	res := make([]complex128, len(v))
	for i, el := range v {
		res[Mod(i+nrot, len(res))] = complex(el, 0)
	}
	return res
}

// Return if a diagonal vector exists without extracting/encoding the vectors
func EncodeDiagBool(X BlockVector, index int, slots int) bool {
	for i := range X {
		if GetDiagBool(X[i], slots, index) {
			return true
		}
	}
	return false
}

// index specifies which diagonal to extract
// applies right-rotation by nrot positions before encoding
func EncodeDiag(cryptoParams *crypto.CryptoParams, X BlockVector, index int, nrot int, level int) (crypto.PlainVector, bool) {
	slots := cryptoParams.GetSlots()

	buf := make([]float64, slots)
	out := make(crypto.PlainVector, len(X))
	anyFlag := false

	for i := range X {
		success := GetDiag(buf, X[i], slots, index)
		if success {
			anyFlag = true
			plaintext := ckks.NewPlaintext(cryptoParams.Params, level)
			plaintext.Scale = cryptoParams.Params.DefaultScale()
			cryptoParams.WithEncoder(func(encoder *ckks.Encoder) error {
				return encoder.Encode(convertToComplex128WithRot(buf, nrot), plaintext)
			})
			out[i] = plaintext
		} else {
			out[i] = nil
		}
	}

	return out, anyFlag
}

func EncodeDiagWithEncoder(cryptoParams *crypto.CryptoParams, X BlockVector, index int, nrot int, level int, enc *ckks.Encoder) (crypto.PlainVector, bool) {
	slots := cryptoParams.GetSlots()

	buf := make([]float64, slots)
	out := make(crypto.PlainVector, len(X))
	anyFlag := false

	for i := range X {
		success := GetDiag(buf, X[i], slots, index)
		if success {
			anyFlag = true
			plaintext := ckks.NewPlaintext(cryptoParams.Params, level)
			plaintext.Scale = cryptoParams.Params.DefaultScale()
			if err := enc.Encode(convertToComplex128WithRot(buf, nrot), plaintext); err != nil {
				panic(err)
			}
			out[i] = plaintext
		} else {
			out[i] = nil
		}
	}

	return out, anyFlag
}

type PlainMatrixDiagCache [][]crypto.PlainVector

func MatMult4StreamPreprocess(cryptoParams *crypto.CryptoParams, gfs *GenoFileStream, maxLevel int, cacheFilePrefix string) {
	gfs.Reset() // Reset to beginning of file just in case

	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))

	m_ct := ((gfs.NumCols() - 1) / uint64(slots)) + 1
	numBlockRows := ((gfs.NumRows() - 1) / uint64(slots)) + 1
	nproc := runtime.GOMAXPROCS(0)

	log.LLvl1(time.Now().Format(time.RFC3339), "MatMult4StreamPreprocess:", "input", gfs.NumRows(), gfs.NumCols(), "numBlockRows", numBlockRows, "numBlockCols", m_ct)

	for bi := 0; bi < int(numBlockRows); bi++ {

		dcs, flag := NewDiagCacheStream(cryptoParams, cacheFilePrefix, bi, true)
		if flag {
			continue
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "gathering submatrix")

		BSlice := make([]BlockI8, m_ct)
		nr := Min((bi+1)*slots, int(gfs.NumRows())) - bi*slots
		for ri := 0; ri < nr; ri++ {

			// Read one row from file
			row := gfs.NextRow()

			// Add slice to each block matrix
			for bj := range BSlice {
				j1 := bj * slots
				j2 := Min((bj+1)*slots, int(gfs.NumCols()))
				nc := j2 - j1
				if ri == 0 {
					BSlice[bj] = NewBlockI8(nr, nc)
				}
				BSlice[bj].data[ri] = row[j1:j2]
			}
		}

		blockVec := make(BlockVector, m_ct)
		for bj := range blockVec {
			blockVec[bj] = Block(BSlice[bj])
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "finding active diagonals")

		// Pre-collect active baby/giant indices
		babyTable := make([]bool, d)
		giantTable := make([]bool, d)
		shiftTable := make([]bool, slots)
		for shift := 0; shift < slots; shift++ {
			if EncodeDiagBool(blockVec, -shift, slots) {
				baby, giant := shift%d, shift/d
				babyTable[baby] = true
				giantTable[giant] = true
				shiftTable[shift] = true
			}
		}

		dcs.SetIndexTables(babyTable, giantTable)

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "extracting and caching diagonals")

		type dataItem struct {
			plainVec crypto.PlainVector
			shift    int
		}

		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}

		diagChannel := make(chan dataItem, 16)

		// Job feeder
		go func() {
			for shift, flag := range shiftTable {
				if flag {
					jobChannels[shift%nproc] <- shift
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Data writer
		var writer sync.WaitGroup
		writer.Add(1)
		go func() {
			defer writer.Done()
			for item := range diagChannel {
				dcs.WriteDiag(item.plainVec, uint32(item.shift))
			}
		}()

		// Data encoders
		var encoderGroup sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			encoderGroup.Add(1)
			go func(thread int) {
				defer encoderGroup.Done()

				enc := ckks.NewEncoder(cryptoParams.Params, cryptoParams.GetPrec())

				for shift := range jobChannels[thread] {
					_, giant := shift%d, shift/d

					plainVec, _ := EncodeDiagWithEncoder(cryptoParams, blockVec, -shift, d*giant, maxLevel, enc)

					ToMontgomeryForm(cryptoParams, plainVec)

					diagChannel <- dataItem{plainVec, shift}
				}
			}(thread)
		}

		encoderGroup.Wait()
		close(diagChannel)

		writer.Wait()
		dcs.Close()
	}

	return
}

func MatMult4StreamCompute(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, maxLevel int, cacheFilePrefix string) crypto.CipherMatrix {
	s := len(A)
	outScale := A[0][0].Scale.Float64() * cryptoParams.Params.DefaultScale().Float64()
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	nproc := runtime.GOMAXPROCS(0)

	numBlockRows := len(A[0])

	log.LLvl1(time.Now().Format(time.RFC3339), "MatMult4StreamCompute")
	if A[0][0].Level() > maxLevel {
		log.LLvl1(time.Now().Format(time.RFC3339), "Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}

	accCache := make([][]CipherVectorAccV2, s)
	accCacheMux := make([][]sync.Mutex, s)
	for i := range accCache {
		accCache[i] = make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups
		accCacheMux[i] = make([]sync.Mutex, d)
	}

	rotCache := make(crypto.CipherMatrix, s)
	for i := range rotCache {
		rotCache[i] = make(crypto.CipherVector, d)
	}

	var m_ct int

	for bi := 0; bi < numBlockRows; bi++ {

		dcs, _ := NewDiagCacheStream(cryptoParams, cacheFilePrefix, bi, false)

		babyTable, giantTable := dcs.GetIndexTables()

		m_ct = int(dcs.vectorLen)

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "generating rotation cache")

		// Dispatcher
		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 32)
		}

		go func() {
			for baby, flag := range babyTable {
				if flag {
					jobChannels[baby%nproc] <- baby
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		Aslice := make(crypto.CipherVector, len(A))
		for i := range A {
			Aslice[i] = A[i][bi]
		}
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, rlwe.NewMemEvaluationKeySet(cryptoParams.Rlk, cryptoParams.RotKs...))

				for baby := range jobChannels[thread] {
					for i := range A {
						rotCache[i][baby] = crypto.RotateRightWithEvaluator(cryptoParams, Aslice[i], -baby, eva)
					}
				}
			}(thread)
		}
		workerGroup.Wait()

		for giant, flag := range giantTable {
			if flag {
				for i := range A {
					if accCache[i][giant].val == nil {
						accCache[i][giant] = NewCipherVectorAccV2(cryptoParams, int(dcs.vectorLen), maxLevel)
					}
				}
			}
		}

		var wg sync.WaitGroup

		type dataItem struct {
			plainVec crypto.PlainVector
			shift    int
		}

		diagChannels := make([]chan dataItem, nproc)
		for i := range diagChannels {
			diagChannels[i] = make(chan dataItem, 8)
		}

		// Data feeder
		go func() {
			for plainVec, shift := dcs.ReadDiag(); plainVec != nil; plainVec, shift = dcs.ReadDiag() {
				diagChannels[shift%nproc] <- dataItem{plainVec, shift}
			}
			for _, c := range diagChannels {
				close(c)
			}
		}()

		// Data processors
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()
				for item := range diagChannels[thread] {
					plainVec, shift := item.plainVec, item.shift
					baby, giant := shift%d, shift/d
					for i := range A {
						accCacheMux[i][giant].Lock()
						CPMultAccWithoutMRedV2(crypto.CipherVector{rotCache[i][baby]}, plainVec, accCache[i][giant])
						accCacheMux[i][giant].Unlock()
					}
				}
			}(thread)
		}
		wg.Wait()
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Postprocessing accumulators")

	out := crypto.CZeroMat(cryptoParams, m_ct, s)
	for i := range out {
		aggChannel := make(chan crypto.CipherVector, 16)

		jobChannels := make([]chan int, nproc)
		for j := range jobChannels {
			jobChannels[j] = make(chan int, 32)
		}

		go func() {
			for l := range accCache[i] {
				if accCache[i][l].val != nil {
					jobChannels[l%nproc] <- l
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		var wg sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, rlwe.NewMemEvaluationKeySet(cryptoParams.Rlk, cryptoParams.RotKs...))

				for l := range jobChannels[thread] {
					cv := ModularReduceV2(cryptoParams, accCache[i][l], outScale)

					if l > 0 { // Giant step alignment
						for j := range cv {
							cv[j] = crypto.RotateRightWithEvaluator(cryptoParams, cv[j], -l*d, eva)
						}
					}

					aggChannel <- cv
				}
			}(thread)
		}

		var aggGroup sync.WaitGroup
		aggGroup.Add(1)
		go func() {
			defer aggGroup.Done()

			eva := ckks.NewEvaluator(cryptoParams.Params, rlwe.NewMemEvaluationKeySet(cryptoParams.Rlk, cryptoParams.RotKs...))

			for cv := range aggChannel {
				for j := range cv {
					eva.Add(out[i][j], cv[j], out[i][j])
				}
			}
		}()

		wg.Wait()
		close(aggChannel)
		aggGroup.Wait()
	}

	return out
}

func MatMult4Stream(cryptoParams *crypto.CryptoParams, A crypto.CipherMatrix, gfs *GenoFileStream, maxLevel int, computeSquaredSum, square bool, nproc int) (crypto.CipherMatrix, []float64, []float64) {
	gfs.Reset() // Reset to beginning of file just in case

	nrow, ncol := gfs.NumRowsToKeep(), gfs.NumColsToKeep()
	if nproc <= 0 { // If nproc is non-positive, use all cores
		nproc = runtime.GOMAXPROCS(0)
	}

	s := len(A)
	outScale := A[0][0].Scale.Float64() * cryptoParams.Params.DefaultScale().Float64()
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))

	//blockB := ToBlockMatrix(B, slots)
	//fmt.Println("blockB dims:", len(blockB), len(blockB[0]))
	m_ct := ((ncol - 1) / uint64(slots)) + 1
	numBlockRows := ((nrow - 1) / uint64(slots)) + 1

	if A[0][0].Level() > maxLevel {
		fmt.Println("Dropping level. Input:", A[0][0].Level())
		A = crypto.DropLevel(cryptoParams, A, maxLevel)
	}
	fmt.Println("A level:", A[0][0].Level())

	accCache := make([][]CipherVectorAccV2, s)
	accCacheMux := make([][]sync.Mutex, s)
	for i := range accCache {
		accCache[i] = make([]CipherVectorAccV2, d) // Cache each of the sqrt(slots) groups, initialize later on-the-fly
		accCacheMux[i] = make([]sync.Mutex, d)
	}

	rotCache := make(crypto.CipherMatrix, s)
	for i := range rotCache {
		rotCache[i] = make(crypto.CipherVector, d)
	}

	var sqSum, sum []float64
	if computeSquaredSum {
		sqSum = make([]float64, ncol)
		sum = make([]float64, ncol)
	}

	for bi := 0; bi < int(numBlockRows); bi++ {

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "gathering submatrix")

		BSlice := make([]BlockI8, m_ct)
		nr := Min((bi+1)*slots, int(nrow)) - bi*slots
		for ri := 0; ri < nr; ri++ {

			// Read one row from file
			row := gfs.NextRow()

			// Replace missing with zeros
			for rj := range row {
				if row[rj] < 0 {
					row[rj] = 0
				}

				if computeSquaredSum {
					sqSum[rj] += float64(row[rj] * row[rj])
					sum[rj] += float64(row[rj])
				}
				if square { // optionnally square the values
					row[rj] = row[rj] * row[rj]
				}
			}

			// Add slice to each block matrix
			for bj := range BSlice {
				j1 := bj * slots
				j2 := Min((bj+1)*slots, int(ncol))
				nc := j2 - j1
				if ri == 0 {
					BSlice[bj] = NewBlockI8(nr, nc)
				}
				BSlice[bj].data[ri] = row[j1:j2]
			}
		}

		blockVec := make(BlockVector, m_ct)
		for bj := range blockVec {
			blockVec[bj] = Block(BSlice[bj])
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "finding active diagonals")

		// Pre-collect active baby/giant indices
		babyTable := make([]bool, d)
		giantTable := make([]bool, d)
		shiftTable := make([]bool, slots)
		for shift := 0; shift < slots; shift++ {
			if EncodeDiagBool(blockVec, -shift, slots) {
				baby, giant := shift%d, shift/d
				babyTable[baby] = true
				giantTable[giant] = true
				shiftTable[shift] = true
			}
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "generating rotation cache")

		log.LLvl1(time.Now().Format(time.RFC3339), "Num procs", nproc)

		// Dispatcher
		jobChannels := make([]chan int, nproc)
		for i := range jobChannels {
			jobChannels[i] = make(chan int, 64)
		}
		go func() {
			index := 0
			for baby, flag := range babyTable {
				if flag {
					jobChannels[index%nproc] <- baby
					index++
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		// Workers
		var workerGroup sync.WaitGroup
		Aslice := make(crypto.CipherVector, len(A))
		for i := range A {
			Aslice[i] = A[i][bi]
		}
		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, rlwe.NewMemEvaluationKeySet(cryptoParams.Rlk, cryptoParams.RotKs...))

				for baby := range jobChannels[thread] {
					for i := range A {
						rotCache[i][baby] = crypto.RotateRightWithEvaluator(cryptoParams, Aslice[i], -baby, eva)
					}
				}
			}(thread)
		}
		workerGroup.Wait()

		for giant, flag := range giantTable {
			if flag {
				for i := range A {
					if accCache[i][giant].val == nil {
						accCache[i][giant] = NewCipherVectorAccV2(cryptoParams, int(m_ct), maxLevel)
					}
				}
			}
		}

		log.LLvl1(time.Now().Format(time.RFC3339), "Block row", bi+1, "/", numBlockRows, "extracting and multiplying diagonals")

		// Extract and encode diagonal vectors
		shiftChannels := make([]chan int, nproc)
		for i := range shiftChannels {
			shiftChannels[i] = make(chan int, 128)
		}

		go func() {
			index := 0
			for shift, flag := range shiftTable {
				if flag {
					if (index+1)%1000 == 0 {
						log.LLvl1(index + 1)
					}
					shiftChannels[index%nproc] <- shift
					index++
				}
			}
			for _, c := range shiftChannels {
				close(c)
			}
		}()

		for thread := 0; thread < nproc; thread++ {
			workerGroup.Add(1)
			go func(thread int) {
				defer workerGroup.Done()

				enc := ckks.NewEncoder(cryptoParams.Params, cryptoParams.GetPrec())

				for shift := range shiftChannels[thread] {
					baby, giant := shift%d, shift/d

					plainVec, _ := EncodeDiagWithEncoder(cryptoParams, blockVec, -shift, d*giant, maxLevel, enc)

					ToMontgomeryForm(cryptoParams, plainVec)

					for i := range A {
						accCacheMux[i][giant].Lock()
						CPMultAccWithoutMRedV2(crypto.CipherVector{rotCache[i][baby]}, plainVec, accCache[i][giant])
						accCacheMux[i][giant].Unlock()
					}
				}
			}(thread)
		}
		workerGroup.Wait()
	}

	log.LLvl1(time.Now().Format(time.RFC3339), "Postprocessing accumulators")

	out := crypto.CZeroMat(cryptoParams, int(m_ct), s)
	for i := range out {
		jobChannels := make([]chan int, nproc)
		for j := range jobChannels {
			jobChannels[j] = make(chan int, 32)
		}

		go func() {
			for l := range accCache[i] {
				if accCache[i][l].val != nil {
					jobChannels[l%nproc] <- l
				}
			}
			for _, c := range jobChannels {
				close(c)
			}
		}()

		aggChannel := make(chan crypto.CipherVector, 8)

		var wg sync.WaitGroup
		for thread := 0; thread < nproc; thread++ {
			wg.Add(1)
			go func(thread int) {
				defer wg.Done()

				eva := ckks.NewEvaluator(cryptoParams.Params, rlwe.NewMemEvaluationKeySet(cryptoParams.Rlk, cryptoParams.RotKs...))

				for l := range jobChannels[thread] {
					cv := ModularReduceV2(cryptoParams, accCache[i][l], outScale)

					if l > 0 { // Giant step alignment
						for j := range cv {
							cv[j] = crypto.RotateRightWithEvaluator(cryptoParams, cv[j], -l*d, eva)
						}
					}

					aggChannel <- cv
				}
			}(thread)
		}

		var aggGroup sync.WaitGroup
		aggGroup.Add(1)
		go func() {
			defer aggGroup.Done()

			eva := ckks.NewEvaluator(cryptoParams.Params, rlwe.NewMemEvaluationKeySet(cryptoParams.Rlk, cryptoParams.RotKs...))

			for cv := range aggChannel {
				for j := range cv {
					eva.Add(out[i][j], cv[j], out[i][j])
				}
			}
		}()

		wg.Wait()
		close(aggChannel)
		aggGroup.Wait()
		// Rescale out[i] back to DefaultScale after accumulation.
		out[i] = crypto.CRescale(cryptoParams, out[i])
	}

	return out, sum, sqSum
}

func MatMult4TransformB(cryptoParams *crypto.CryptoParams, B *mat.Dense) PlainMatrixDiagCache {
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	blockB := ToBlockMatrix(B, slots)

	cache := make(PlainMatrixDiagCache, len(blockB))

	for bi := range blockB {

		cache[bi] = make([]crypto.PlainVector, slots)

		for shift := 0; shift < slots; shift++ {
			giant := int(shift / d)
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, d*giant, cryptoParams.Params.MaxLevel())
			if !flag {
				cache[bi][shift] = nil
			} else {
				ToMontgomeryForm(cryptoParams, plainVec)
				cache[bi][shift] = plainVec
			}
		}
	}

	return cache
}

// Caches a banded matrix B with specified upper and lower bandwidth.
// Assumes that B is at most a (slots x slots) matrix, for `slots` the number of slots in a ciphertext.
func MatMult4TransformBandedB(cryptoParams *crypto.CryptoParams, B *mat.Dense, upperBandwidth, lowerBandwidth int) PlainMatrixDiagCache {
	slots := cryptoParams.GetSlots()
	d := int(math.Ceil(math.Sqrt(float64(slots))))
	blockB := ToBlockMatrix(B, slots)

	cache := make(PlainMatrixDiagCache, len(blockB))

	// The band is the diagonals with shift indices in [0, upperBandwidth + 1) and [slots - lowerBandwidth, slots)
	topDiagBound := upperBandwidth + 1 // The +1 is because we also want to include the main diagonal, which is index 0.
	leftDiagBound := slots - lowerBandwidth
	// log.LLvl1("topDiagBound:", topDiagBound, "; leftDiagBound", leftDiagBound)

	for bi := range blockB {

		cache[bi] = make([]crypto.PlainVector, slots)

		for shift := 0; shift < slots; shift++ {
			if topDiagBound <= shift && shift < leftDiagBound {
				// The diagonal this shift corresponds falls outside of the band.
				continue
			}

			giant := int(shift / d)
			plainVec, flag := EncodeDiag(cryptoParams, blockB[bi], -shift, d*giant, cryptoParams.Params.MaxLevel())
			if !flag {
				cache[bi][shift] = nil
			} else {
				ToMontgomeryForm(cryptoParams, plainVec)
				cache[bi][shift] = plainVec
			}
		}
	}

	return cache
}

func CMultMatInnerProd(cryptoParams *crypto.CryptoParams, M, N crypto.CipherMatrix, numThreads int) crypto.CipherMatrix {
	var mutex sync.Mutex
	log.LLvl1("Matrix multiplication, result size ", len(M), "x", len(N))
	result, _, _, err := crypto.InitEncryptedMatrix(cryptoParams, len(M), len(N))
	if err != nil {
		log.Fatal(err)
	}
	vparallelize := int(math.Ceil(float64(len(M)) / float64(numThreads)))
	wg := sync.WaitGroup{}
	for MRow := 0; MRow < len(M); MRow += vparallelize {
		wg.Add(1)
		go func(iT int) {
			defer wg.Done()
			for k := 0; k < vparallelize && (k+iT < len(M)); k++ {
				for NCol := 0; NCol < len(N); NCol++ {
					multi := crypto.InnerProd(cryptoParams, M[k+iT], crypto.CopyEncryptedVector(N[NCol]))
					multiMasked := crypto.Mask(cryptoParams, multi, NCol, false)
					mutex.Lock()
					result[k+iT] = crypto.CAdd(cryptoParams, crypto.CipherVector{multiMasked}, result[k+iT])
					mutex.Unlock()
				}
			}
			log.LLvl1("processed row ", iT, "to", iT+vparallelize, "of", len(M), "rows")
		}(MRow)
	}
	wg.Wait()

	return result
}

// Matrix multiplication using inner products, with the result being a vector
func CMultMatInnerProdVector(cryptoParams *crypto.CryptoParams, M crypto.CipherMatrix, N crypto.CipherVector, MCols, numThreads int) crypto.CipherVector {
	var mutex sync.Mutex
	log.LLvl1("Multiply Matrix with column vector: result size will be a vector of length", len(M))
	result := crypto.CZeros(cryptoParams, int(math.Ceil(float64(len(M))/float64(cryptoParams.GetSlots()))))

	// row col approach
	vparallelize := int(math.Ceil(float64(len(M)) / float64(numThreads)))
	wg := sync.WaitGroup{}
	maskClear := make([]float64, cryptoParams.GetSlots()*len(M[0]))
	for i := 0; i < MCols; i++ {
		maskClear[i] = 1
	}
	log.LLvl1("NUmber of values kept in mask", MCols)
	maskEncoded, _ := crypto.EncodeFloatVector(cryptoParams, maskClear)
	N = crypto.CPMult(cryptoParams, N, maskEncoded)

	for MRow := 0; MRow < len(M); MRow += vparallelize {
		wg.Add(1)
		go func(iT int) {
			defer wg.Done()
			for k := 0; k < vparallelize && (k+iT < len(M)); k++ {
				// assume len(M) smaller than number of slots
				if len(M) > cryptoParams.GetSlots() {
					log.LLvl1("ERROR ! Number of rows in M is larger than the number of slots")
				}
				Mcurrent := crypto.CPMult(cryptoParams, M[k+iT], maskEncoded)
				multi := crypto.InnerProd(cryptoParams, Mcurrent, crypto.CopyEncryptedVector(N))
				multiMasked := crypto.Mask(cryptoParams, multi, k+iT, false)
				mutex.Lock()
				result = crypto.CAdd(cryptoParams, crypto.CipherVector{multiMasked}, result)
				mutex.Unlock()
			}
			log.LLvl1("processed row ", iT, "to", iT+vparallelize, "of", len(M), "rows")
		}(MRow)
	}
	wg.Wait()

	return result
}

// Matrix multiplication between column encrypted matrices
func CMultMatColTimesColToCol(cryptoParams *crypto.CryptoParams, M, N crypto.CipherMatrix, numRowsM,
	numThreads int) crypto.CipherMatrix {

	// numRowsM = numInds, numColsN = len(C)
	log.LLvl1("result size ", numRowsM, "x", len(N))
	result, _, _, err := crypto.InitEncryptedMatrix(cryptoParams, len(N), numRowsM)
	if err != nil {
		log.Fatal(err)
	}

	vparallelize := int(math.Ceil(float64(len(M)) / float64(numThreads)))
	mutex := sync.Mutex{}
	wg := sync.WaitGroup{}

	for MCol := 0; MCol < len(M); MCol += vparallelize {
		wg.Add(1)
		go func(iT int) {
			defer wg.Done()
			for k := 0; k < vparallelize && (k+iT < len(M)); k++ {
				for Ncol := 0; Ncol < len(N); Ncol++ {
					elemRep := crypto.CMask(cryptoParams, N[Ncol], k+iT, false)
					elemRepCiph := crypto.InnerSumAll(cryptoParams, elemRep)
					elemRepNew := make(crypto.CipherVector, len(M[k+iT]))
					for j := range elemRepNew {
						elemRepNew[j] = elemRepCiph.CopyNew()
					}
					multi := crypto.CMult(cryptoParams, elemRepNew, M[k+iT])
					mutex.Lock()
					result[Ncol] = crypto.CAdd(cryptoParams, multi, result[Ncol])
					mutex.Unlock()
				}
			}
			log.LLvl1("processed col ", iT, "to", iT+vparallelize, "of", len(M), "cols")
		}(MCol)

	}
	wg.Wait()

	return result
}

// Matrix multiplication between column and row encrypted matrices
func CMultMatColTimesRowToCol(cryptoParams *crypto.CryptoParams, M, N crypto.CipherMatrix, numRowsM, numColsN,
	numThreads int) crypto.CipherMatrix {

	// numRowsM = numInds, numColsN = len(C)
	log.LLvl1("Matrix multiplication, result size ", numRowsM, "x", numColsN)
	result, _, _, err := crypto.InitEncryptedMatrix(cryptoParams, numColsN, numRowsM)
	if err != nil {
		log.Fatal(err)
	}
	vparallelize := int(math.Ceil(float64(len(M)) / float64(numThreads)))
	mutex := sync.Mutex{}
	wg := sync.WaitGroup{}

	for MCol := 0; MCol < len(M); MCol += vparallelize {
		wg.Add(1)
		go func(iT int) {
			defer wg.Done()
			for k := 0; k < vparallelize && (k+iT < len(M)); k++ {
				for Ncol := 0; Ncol < len(N); Ncol++ {
					elemRep := crypto.CMask(cryptoParams, N[k+iT], Ncol, false)
					elemRepCiph := crypto.InnerSumAll(cryptoParams, elemRep)
					elemRepNew := make(crypto.CipherVector, len(M[k+iT]))
					for j := range elemRepNew {
						elemRepNew[j] = elemRepCiph.CopyNew()
					}
					multi := crypto.CMult(cryptoParams, elemRepNew, M[k+iT])
					mutex.Lock()
					result[Ncol] = crypto.CAdd(cryptoParams, multi, result[Ncol])
					mutex.Unlock()
				}
			}
			log.LLvl1("Processed col ", iT, "to", iT+vparallelize, "of", len(M), "rows")
		}(MCol)
	}
	wg.Wait()

	return result
}
