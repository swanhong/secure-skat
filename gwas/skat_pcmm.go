package gwas

import (
	"fmt"
	"time"

	mpc_core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/lintrans"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"go.dedis.ch/onet/v3/log"
)

const pcmmInputLevel = 5

// pcmmWindowBytes is the encoded transform plus the raw diagonals alive during encoding.
func pcmmWindowBytes(params ckks.Parameters, bucket GeneBatchBucket) uint64 {
	levels := pcmmInputLevel + 1 + params.MaxLevelP() + 1
	encoded := uint64(bucket.P) * uint64(levels) * uint64(params.N()) * 8
	raw := uint64(bucket.P) * uint64(params.MaxSlots()) * 8
	return encoded + raw
}

func pcmmParameters(params ckks.Parameters, bucket GeneBatchBucket) lintrans.Parameters {
	diagonals := make([]int, bucket.P)
	for d := range diagonals {
		diagonals[d] = bucket.L * d
	}
	return lintrans.Parameters{
		DiagonalsIndexList:        diagonals,
		LevelQ:                    pcmmInputLevel,
		LevelP:                    params.MaxLevelP(),
		Scale:                     params.GetOptimalScalingFactor(params.DefaultScale(), params.DefaultScale(), pcmmInputLevel),
		LogDimensions:             params.LogMaxDimensions(),
		LogBabyStepGiantStepRatio: 0,
	}
}

func pcmmGaloisElements(params ckks.Parameters, bucket GeneBatchBucket) []uint64 {
	return lintrans.GaloisElements(params, pcmmParameters(params, bucket))
}

func windowGramDiagonals(bucket GeneBatchBucket, window GeneBatchWindow, local []windowLocalContraction) lintrans.Diagonals[float64] {
	if len(local) != len(window.Tiles) {
		panic("window local contraction count mismatch")
	}
	diagonals := make(lintrans.Diagonals[float64], bucket.P)
	for d := 0; d < bucket.P; d++ {
		values := make([]float64, bucket.P*bucket.L)
		for tileIndex, tile := range window.Tiles {
			gamma := local[tileIndex].Gamma
			if gamma == nil {
				continue
			}
			for row := 0; row < tile.Variants; row++ {
				col := (row + d) % bucket.P
				if col >= tile.Variants {
					continue
				}
				for lane := 0; lane < window.H; lane++ {
					values[row*bucket.L+tile.LaneBase+lane] = gamma.At(row, col)
				}
			}
		}
		diagonals[bucket.L*d] = values
	}
	return diagonals
}

func (ast *AssocTest) encodeWindowGram(bucket GeneBatchBucket, window GeneBatchWindow, local []windowLocalContraction) lintrans.LinearTransformation {
	transform := lintrans.NewTransformation(ast.general.cps.Params, pcmmParameters(ast.general.cps.Params, bucket))
	if err := ast.general.cps.WithEncoder(func(encoder *ckks.Encoder) error {
		return lintrans.Encode(encoder, windowGramDiagonals(bucket, window, local), transform)
	}); err != nil {
		panic(err)
	}
	return transform
}

// packedGramAction evaluates each gene's local Gamma on the same packed ciphertext.
func (ast *AssocTest) packedGramAction(bucket GeneBatchBucket, window GeneBatchWindow, local []windowLocalContraction, values []mpc_core.RMat, rhs int) []mpc_core.RMat {
	mpcObj := ast.general.mpcObj[0]
	started := time.Now()
	defer func() {
		if mpcObj.GetPid() == mpcObj.GetHubPid() {
			log.LLvl1(fmt.Sprintf("[skat_fed] PCMM P=%d genes=%d rhs=%d groups=%d elapsed=%v",
				bucket.P, len(window.Tiles), rhs, windowGroupCount(window, rhs), time.Since(started).Round(time.Millisecond)))
		}
	}()
	input := ast.windowSharesToCiphertexts(bucket, window, values, rhs)
	output := make(crypto.CipherVector, len(input))
	if mpcObj.GetPid() > 0 {
		input = crypto.DropLevel(ast.general.cps, crypto.CipherMatrix{input}, pcmmInputLevel)[0]
		transform := ast.encodeWindowGram(bucket, window, local)
		if err := ast.general.cps.WithEvaluator(func(evaluator *ckks.Evaluator) error {
			linear := lintrans.NewEvaluator(evaluator)
			for group, ct := range input {
				var err error
				if output[group], err = linear.EvaluateNew(ct, transform); err != nil {
					return err
				}
				if err = evaluator.Rescale(output[group], output[group]); err != nil {
					return err
				}
			}
			return nil
		}); err != nil {
			panic(err)
		}
		// Hide each party's deterministic local transform before aggregation.
		zero := crypto.CZeros(ast.general.cps, len(output))
		zero = crypto.DropLevel(ast.general.cps, crypto.CipherMatrix{zero}, output[0].Level())[0]
		output = crypto.CAdd(ast.general.cps, output, zero)
		output = mpcObj.Network.AggregateCVec(ast.general.cps, output)
	}
	return ast.windowCiphertextsToShares(bucket, window, output, rhs)
}
