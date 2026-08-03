package gwas

import (
	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func packFlatVariants(genes [][]float64) []float64 {
	n := 0
	for _, gene := range genes {
		n += len(gene)
	}
	packed := make([]float64, 0, n)
	for _, gene := range genes {
		packed = append(packed, gene...)
	}
	return packed
}

func scatterFlatVariants(packed []float64, sizes []int) [][]float64 {
	genes := make([][]float64, len(sizes))
	offset := 0
	for gene, size := range sizes {
		genes[gene] = append([]float64(nil), packed[offset:offset+size]...)
		offset += size
	}
	if offset != len(packed) {
		panic("flat variant size mismatch")
	}
	return genes
}

func windowGroupCount(window GeneBatchWindow, rhs int) int {
	if rhs <= 0 {
		panic("window RHS count must be positive")
	}
	return (rhs + window.H - 1) / window.H
}

// packWindowShares maps tile-local [variant][RHS] shares to group-major CKKS slots.
func packWindowShares(rtype mpc_core.RElem, bucket GeneBatchBucket, window GeneBatchWindow, tileValues []mpc_core.RMat, rhs int) mpc_core.RVec {
	if len(tileValues) != len(window.Tiles) {
		panic("window tile count mismatch")
	}
	slots := bucket.P * bucket.L
	packed := mpc_core.InitRVec(rtype.Zero(), windowGroupCount(window, rhs)*slots)
	for tileIndex, tile := range window.Tiles {
		values := tileValues[tileIndex]
		if len(values) != tile.Variants {
			panic("window variant count mismatch")
		}
		for variant, row := range values {
			if len(row) != rhs {
				panic("window RHS width mismatch")
			}
			for column, value := range row {
				group := column / window.H
				lane := tile.LaneBase + column%window.H
				packed[group*slots+variant*bucket.L+lane] = value
			}
		}
	}
	return packed
}

func scatterWindowShares(rtype mpc_core.RElem, bucket GeneBatchBucket, window GeneBatchWindow, packed mpc_core.RVec, rhs int) []mpc_core.RMat {
	slots := bucket.P * bucket.L
	if len(packed) != windowGroupCount(window, rhs)*slots {
		panic("packed window size mismatch")
	}
	values := make([]mpc_core.RMat, len(window.Tiles))
	for tileIndex, tile := range window.Tiles {
		values[tileIndex] = mpc_core.InitRMat(rtype.Zero(), tile.Variants, rhs)
		for variant := 0; variant < tile.Variants; variant++ {
			for column := 0; column < rhs; column++ {
				group := column / window.H
				lane := tile.LaneBase + column%window.H
				values[tileIndex][variant][column] = packed[group*slots+variant*bucket.L+lane]
			}
		}
	}
	return values
}

func windowActiveMask(bucket GeneBatchBucket, window GeneBatchWindow, rhs int) []float64 {
	slots := bucket.P * bucket.L
	mask := make([]float64, windowGroupCount(window, rhs)*slots)
	for _, tile := range window.Tiles {
		for variant := 0; variant < tile.Variants; variant++ {
			for column := 0; column < rhs; column++ {
				group := column / window.H
				lane := tile.LaneBase + column%window.H
				mask[group*slots+variant*bucket.L+lane] = 1
			}
		}
	}
	return mask
}

func (ast *AssocTest) applyPackedMask(values crypto.CipherVector, mask []float64) crypto.CipherVector {
	if len(mask) != len(values)*ast.general.cps.GetSlots() {
		panic("packed mask size mismatch")
	}
	if ast.general.mpcObj[0].GetPid() == 0 {
		return make(crypto.CipherVector, len(values))
	}
	encoded, _ := crypto.EncodeFloatVector(ast.general.cps, mask)
	return crypto.CPMult(ast.general.cps, values, encoded)
}

func (ast *AssocTest) maskWindowActive(bucket GeneBatchBucket, window GeneBatchWindow, values crypto.CipherVector, rhs int) crypto.CipherVector {
	return ast.applyPackedMask(values, windowActiveMask(bucket, window, rhs))
}

// sumWindowRows masks padding and sums P variant rows per lane.
func (ast *AssocTest) sumWindowRows(bucket GeneBatchBucket, window GeneBatchWindow, values crypto.CipherVector, rhs int) crypto.CipherVector {
	if ast.general.mpcObj[0].GetPid() == 0 {
		return make(crypto.CipherVector, len(values))
	}
	values = ast.maskWindowActive(bucket, window, values, rhs)
	out := make(crypto.CipherVector, len(values))
	if err := ast.general.cps.WithEvaluator(func(eval *ckks.Evaluator) error {
		for i, ct := range values {
			out[i] = ct.CopyNew()
			if err := eval.InnerSum(ct, bucket.L, bucket.P, out[i]); err != nil {
				return err
			}
		}
		return nil
	}); err != nil {
		panic(err)
	}
	return out
}

func (ast *AssocTest) windowSharesToCiphertexts(bucket GeneBatchBucket, window GeneBatchWindow, tileValues []mpc_core.RMat, rhs int) crypto.CipherVector {
	mpcObj := ast.general.mpcObj[0]
	values := mpcObj.SSToCVec(ast.general.cps, packWindowShares(mpcObj.GetRType(), bucket, window, tileValues, rhs))
	if mpcObj.GetPid() == 0 {
		return make(crypto.CipherVector, windowGroupCount(window, rhs))
	}
	return values
}

func (ast *AssocTest) windowCiphertextsToShares(bucket GeneBatchBucket, window GeneBatchWindow, values crypto.CipherVector, rhs int) []mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	if len(values) != windowGroupCount(window, rhs) {
		panic("window ciphertext count mismatch")
	}
	slots := bucket.P * bucket.L
	packed := mpcObj.CVecToSS(ast.general.cps, mpcObj.GetRType(), values, -1, len(values), len(values)*slots)
	return scatterWindowShares(mpcObj.GetRType(), bucket, window, packed, rhs)
}

// windowSumsToShares converts row-0 lane sums and scatters one value per tile/RHS.
func (ast *AssocTest) windowSumsToShares(bucket GeneBatchBucket, window GeneBatchWindow, values crypto.CipherVector, rhs int) []mpc_core.RVec {
	mpcObj := ast.general.mpcObj[0]
	if len(values) != windowGroupCount(window, rhs) {
		panic("window ciphertext count mismatch")
	}
	cm := make(crypto.CipherMatrix, len(values))
	for group, ct := range values {
		cm[group] = crypto.CipherVector{ct}
	}
	packed := mpcObj.CMatToSS(ast.general.cps, mpcObj.GetRType(), cm, -1, len(values), 1, bucket.L)
	out := make([]mpc_core.RVec, len(window.Tiles))
	for tileIndex, tile := range window.Tiles {
		out[tileIndex] = mpc_core.InitRVec(mpcObj.GetRType().Zero(), rhs)
		if tile.Variants == 0 {
			continue
		}
		for column := 0; column < rhs; column++ {
			group := column / window.H
			lane := tile.LaneBase + column%window.H
			out[tileIndex][column] = packed[group][lane]
		}
	}
	return out
}
