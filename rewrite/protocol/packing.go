package protocol

func PackGeneBatch[T any](
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	values [][]T,
) []T {
	nu := cryptoParams.Slots / batch.W
	h := nu / len(batch.GeneIndices)
	packed := make([]T, cryptoParams.Slots)

	for pos, geneIndex := range batch.GeneIndices {
		laneBase := pos * h
		variantCount := dataParams.Genes[geneIndex].VariantCount
		for variantIndex := 0; variantIndex < variantCount; variantIndex++ {
			slot := variantIndex*nu + laneBase
			packed[slot] = values[pos][variantIndex]
		}
	}
	return packed
}

func UnpackGeneBatch[T any](
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
	packed []T,
) [][]T {
	nu := cryptoParams.Slots / batch.W
	h := nu / len(batch.GeneIndices)
	values := make([][]T, len(batch.GeneIndices))

	for pos, geneIndex := range batch.GeneIndices {
		laneBase := pos * h
		variantCount := dataParams.Genes[geneIndex].VariantCount
		values[pos] = make([]T, variantCount)
		for variantIndex := 0; variantIndex < variantCount; variantIndex++ {
			slot := variantIndex*nu + laneBase
			values[pos][variantIndex] = packed[slot]
		}
	}
	return values
}

func ActiveMask(
	dataParams DataParams,
	cryptoParams CryptoParams,
	batch GeneBatch,
) []float64 {
	activeValues := make([][]float64, len(batch.GeneIndices))

	for pos, geneIndex := range batch.GeneIndices {
		variantCount := dataParams.Genes[geneIndex].VariantCount
		activeValues[pos] = make([]float64, variantCount)

		for variantIndex := range activeValues[pos] {
			activeValues[pos][variantIndex] = 1.0
		}
	}

	return PackGeneBatch(dataParams, cryptoParams, batch, activeValues)
}
