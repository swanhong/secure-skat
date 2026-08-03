package gwas

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/json"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/hhcho/sfgwas/mpc"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"go.dedis.ch/onet/v3/log"
)

const (
	geneBatchRaw        = "raw"
	geneBatchExact      = "exact"
	geneBatchHutchinson = "hutchinson"
)

type GeneBatchManifest struct {
	Slots          int               `json:"slots"`
	Probes         int               `json:"probes"`
	Buckets        []GeneBatchBucket `json:"buckets"`
	GaloisElements []uint64          `json:"galois_elements"`
}

type GeneBatchBucket struct {
	Mode           string            `json:"mode"`
	P              int               `json:"p"`
	L              int               `json:"lanes"`
	TransformBytes uint64            `json:"-"`
	Windows        []GeneBatchWindow `json:"windows"`
}

type GeneBatchWindow struct {
	H     int             `json:"rhs_per_gene"`
	Tiles []GeneBatchTile `json:"tiles"`
}

type GeneBatchTile struct {
	Gene     int    `json:"gene"`
	GeneID   string `json:"gene_id"`
	Variants int    `json:"variants"`
	LaneBase int    `json:"lane_base"`
}

func loadPublicGeneSizes(path string, n int) []int {
	fields := readPublicFields(path, n)
	sizes := make([]int, n)
	for i := range sizes {
		var err error
		sizes[i], err = strconv.Atoi(fields[i])
		if err != nil || sizes[i] < 0 {
			panic(fmt.Sprintf("invalid gene size in %s at line %d", path, i+1))
		}
	}
	return sizes
}

func readPublicFields(path string, n int) []string {
	b, err := os.ReadFile(path)
	if err != nil {
		panic(err)
	}
	fields := strings.Fields(string(b))
	if len(fields) != n {
		panic(fmt.Sprintf("%s has %d entries, want %d", path, len(fields), n))
	}
	return fields
}

func genePaddingSize(m int) int {
	p := 1
	for p < m {
		p <<= 1
	}
	return p
}

func buildGeneBatchManifest(params ckks.Parameters, geneIDs []string, publicSizes []int, probes int) GeneBatchManifest {
	if len(geneIDs) != len(publicSizes) {
		panic("gene ID and public-size counts differ")
	}
	if probes < 0 {
		panic("SKAT probe count cannot be negative")
	}

	slots := params.MaxSlots()
	paddings := make([]int, len(geneIDs))
	modes := make([]string, len(geneIDs))
	seenIDs := make(map[string]bool, len(geneIDs))
	for gene, m := range publicSizes {
		if m < 0 || geneIDs[gene] == "" || seenIDs[geneIDs[gene]] {
			panic(fmt.Sprintf("invalid public gene %d", gene))
		}
		seenIDs[geneIDs[gene]] = true
		paddings[gene] = genePaddingSize(m)
		if paddings[gene] > slots {
			panic(fmt.Sprintf("gene %s requires P=%d, exceeds %d CKKS slots", geneIDs[gene], paddings[gene], slots))
		}
		modes[gene] = geneBatchRaw
		if probes > 0 && probes >= m {
			modes[gene] = geneBatchExact
		} else if probes > 0 {
			modes[gene] = geneBatchHutchinson
		}
	}

	manifest := GeneBatchManifest{Slots: slots, Probes: probes}
	galois := make(map[uint64]bool)
	for _, mode := range []string{geneBatchRaw, geneBatchExact, geneBatchHutchinson} {
		for p := 1; p <= slots; p <<= 1 {
			var tiles []GeneBatchTile
			for gene, m := range publicSizes {
				if modes[gene] == mode && paddings[gene] == p {
					tiles = append(tiles, GeneBatchTile{Gene: gene, GeneID: geneIDs[gene], Variants: m})
				}
			}
			if len(tiles) == 0 {
				continue
			}

			lanes := slots / p
			bucket := GeneBatchBucket{Mode: mode, P: p, L: lanes}
			if mode == geneBatchHutchinson && params.MaxLevel() >= pcmmInputLevel {
				bucket.TransformBytes = pcmmWindowBytes(params, bucket)
			}
			for start := 0; start < len(tiles); start += lanes {
				window := GeneBatchWindow{Tiles: tiles[start:min(start+lanes, len(tiles))]}
				window.H = lanes / len(window.Tiles)
				for gene := range window.Tiles {
					window.Tiles[gene].LaneBase = gene * window.H
				}
				bucket.Windows = append(bucket.Windows, window)
			}
			manifest.Buckets = append(manifest.Buckets, bucket)
			for _, galEl := range params.GaloisElementsForInnerSum(lanes, p) {
				if galEl != 1 {
					galois[galEl] = true
				}
			}
			if mode == geneBatchHutchinson && params.MaxLevel() >= pcmmInputLevel {
				for _, galEl := range pcmmGaloisElements(params, bucket) {
					if galEl != 1 {
						galois[galEl] = true
					}
				}
			}
		}
	}
	for galEl := range galois {
		manifest.GaloisElements = append(manifest.GaloisElements, galEl)
	}
	sort.Slice(manifest.GaloisElements, func(i, j int) bool {
		return manifest.GaloisElements[i] < manifest.GaloisElements[j]
	})
	validateGeneBatchManifest(manifest, len(geneIDs))
	return manifest
}

func validateGeneBatchManifest(manifest GeneBatchManifest, nGenes int) {
	seen := make([]bool, nGenes)
	for _, bucket := range manifest.Buckets {
		for _, window := range bucket.Windows {
			if len(window.Tiles) == 0 || len(window.Tiles) > bucket.L || window.H != bucket.L/len(window.Tiles) {
				panic("invalid batch window")
			}
			for gene, tile := range window.Tiles {
				if tile.Gene < 0 || tile.Gene >= nGenes || seen[tile.Gene] {
					panic("duplicate or invalid gene in batch manifest")
				}
				if tile.LaneBase != gene*window.H || tile.LaneBase+window.H > bucket.L {
					panic("lane collision in batch manifest")
				}
				seen[tile.Gene] = true
			}
		}
	}
	for gene, ok := range seen {
		if !ok {
			panic(fmt.Sprintf("gene %d missing from batch manifest", gene))
		}
	}
}

func (manifest GeneBatchManifest) hash() [sha256.Size]byte {
	b, err := json.Marshal(manifest)
	if err != nil {
		panic(err)
	}
	return sha256.Sum256(b)
}

func syncGeneBatchManifest(mpcObj *mpc.MPC, hash [sha256.Size]byte) {
	net := mpcObj.Network
	pid, hub := net.GetPid(), mpcObj.GetHubPid()
	words := make([]uint64, sha256.Size/8)
	for i := range words {
		words[i] = binary.LittleEndian.Uint64(hash[i*8 : (i+1)*8])
	}

	if pid != hub {
		net.SendIntVector(words, hub)
		if net.ReceiveInt(hub) != 1 {
			panic("gene batch manifest mismatch")
		}
		return
	}

	match := true
	for party := 0; party < net.GetNParty(); party++ {
		if party == hub {
			continue
		}
		got := net.ReceiveIntVector(len(words), party)
		for i := range words {
			match = match && got[i] == words[i]
		}
	}
	verdict := 0
	if match {
		verdict = 1
	}
	for party := 0; party < net.GetNParty(); party++ {
		if party != hub {
			net.SendInt(verdict, party)
		}
	}
	if !match {
		panic("gene batch manifest mismatch")
	}
}

func (ast *AssocTest) prepareGeneBatchManifest() (GeneBatchManifest, []int) {
	mpcObj := ast.general.mpcObj[0]
	n, pid, hub := ast.general.config.GenoNumBlocks, mpcObj.GetPid(), mpcObj.GetHubPid()
	publicSizes := loadPublicGeneSizes(ast.general.config.GenoBlockSizeFile, n)
	geneIDs := readPublicFields(ast.general.config.GeneIDFile, n)
	manifest := buildGeneBatchManifest(ast.general.cps.Params, geneIDs, publicSizes, ast.general.config.SkatPValueProbes)
	hash := manifest.hash()
	syncGeneBatchManifest(mpcObj, hash)
	if pid == hub {
		logGeneBatchManifest(manifest, hash)
	}
	return manifest, publicSizes
}

func logGeneBatchManifest(manifest GeneBatchManifest, hash [sha256.Size]byte) {
	mHistogram := make(map[int]int)
	for _, bucket := range manifest.Buckets {
		genes, padding := 0, 0
		for _, window := range bucket.Windows {
			for _, tile := range window.Tiles {
				genes++
				padding += bucket.P - tile.Variants
				mHistogram[tile.Variants]++
			}
		}
		log.LLvl1(fmt.Sprintf("[skat_fed] manifest mode=%s P=%d genes=%d lanes=%d windows=%d padding=%d transform_mib=%.1f",
			bucket.Mode, bucket.P, genes, bucket.L, len(bucket.Windows), padding, float64(bucket.TransformBytes)/(1<<20)))
	}
	sizes := make([]int, 0, len(mHistogram))
	for size := range mHistogram {
		sizes = append(sizes, size)
	}
	sort.Ints(sizes)
	parts := make([]string, len(sizes))
	for i, size := range sizes {
		parts[i] = fmt.Sprintf("%d:%d", size, mHistogram[size])
	}
	log.LLvl1(fmt.Sprintf("[skat_fed] manifest hash=%x m_pub={%s} galois_elements=%d",
		hash[:8], strings.Join(parts, ","), len(manifest.GaloisElements)))
}
