package protocol

import "fmt"

type Gene struct {
	GeneID       string
	VariantCount int
}

type DataParams struct {
	Genes          []Gene
	NA             int
	NB             int
	N              int
	C              int
	PhenotypeCount int
}

type GeneBatch struct {
	W           int
	GeneIndices []int
}

type CryptoParams struct {
	Slots   int
	R       int
	Batches []GeneBatch
}

func BuildCryptoParams(dataParams DataParams, probeCount int, slots int) (CryptoParams, error) {
	if slots < 1 {
		return CryptoParams{}, fmt.Errorf("slots must be a positive integer")
	}

	if probeCount < 1 {
		return CryptoParams{}, fmt.Errorf("probeCount must be a positive integer")
	}

	genesByWidth := make(map[int][]int)
	maxWidth := 1

	for geneIndex, gene := range dataParams.Genes {
		width := 1
		for width < gene.VariantCount {
			width *= 2
		}
		if width > slots {
			return CryptoParams{}, fmt.Errorf("gene %s has %d variants, which exceeds the number of slots %d", gene.GeneID, gene.VariantCount, slots)
		}
		genesByWidth[width] = append(genesByWidth[width], geneIndex)
		if width > maxWidth {
			maxWidth = width
		}
	}

	cryptoParams := CryptoParams{
		Slots:   slots,
		R:       probeCount,
		Batches: []GeneBatch{},
	}

	for width := 1; width <= maxWidth; width *= 2 {
		geneIndices := genesByWidth[width]
		if len(geneIndices) == 0 {
			continue
		}
		genesPerBatch := slots / width
		for start := 0; start < len(geneIndices); start += genesPerBatch {
			end := min(start+genesPerBatch, len(geneIndices))
			batchIndices := append([]int(nil), geneIndices[start:end]...)
			cryptoParams.Batches = append(cryptoParams.Batches, GeneBatch{
				W:           width,
				GeneIndices: batchIndices,
			})
		}
	}

	return cryptoParams, nil
}
