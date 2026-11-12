package jointstat

import (
	"crypto/rand"
	"encoding/csv"
	"fmt"
	"math"
	"net"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rvas/workflow"
)

type fixtureVariant struct {
	Key     string
	Calls   [2][]int
	Cohort  [2]int
	Private bool
}

func fixture() ([]fixtureVariant, [][]int) {
	variants := make([]fixtureVariant, 6)
	for j := range variants {
		variants[j].Key = fmt.Sprintf("21:%d:A:G", j+1)
		for s := 0; s < 2; s++ {
			variants[j].Calls[s] = make([]int, 24)
		}
	}
	// Shared carriers, overlapping within a gene; B must not count variant carriers twice.
	for s := 0; s < 2; s++ {
		for i := 0; i < 20; i++ {
			variants[0].Calls[s][i] = 1
		}
		variants[0].Cohort[s] = 20
	}
	for i := 0; i < 8; i++ {
		variants[1].Calls[0][i] = 1
	}
	variants[1].Cohort = [2]int{8, 0}
	for i := 0; i < 12; i++ {
		variants[2].Calls[1][i] = 1
	}
	variants[2].Cohort = [2]int{0, 12}
	// True B-private variant appears in two genes; C must count it once.
	variants[3].Private = true
	for i := 0; i < 24; i++ {
		variants[3].Calls[1][i] = 1
	}
	variants[3].Cohort = [2]int{0, 24}
	// ALT-major with missing call: missing must not become two REF alleles.
	for s := 0; s < 2; s++ {
		for i := 0; i < 24; i++ {
			variants[4].Calls[s][i] = 2
		}
		variants[4].Calls[s][0] = 1
		variants[4].Calls[s][1] = -9
		variants[4].Cohort[s] = 398
	}
	// Exactly on the pooled MAF=0.01 boundary: AC=8 out of 800.
	variants[5].Calls[0][0] = 2
	variants[5].Calls[1][0] = 1
	variants[5].Cohort = [2]int{4, 4}
	return variants, [][]int{{0, 1}, {0, 2, 3}, {3}, {}, {4, 5}}
}
func writeFixture(t *testing.T, site string, pid int, variants []fixtureVariant, groups [][]int) {
	t.Helper()
	var genes []geneInput
	var publicCounts, privateCounts []string
	for g, indices := range groups {
		gene := geneInput{ID: fmt.Sprintf("g%d", g)}
		for side := 0; side < 2; side++ {
			if side == 1 && pid == 1 {
				continue
			}
			var selected []fixtureVariant
			var counts []string
			for _, j := range indices {
				v := variants[j]
				if v.Private != (side == 1) {
					continue
				}
				selected = append(selected, v)
				counts = append(counts, fmt.Sprint(v.Cohort[pid-1]))
				if side == 0 {
					gene.Public = append(gene.Public, v.Key)
				} else {
					gene.Private = append(gene.Private, v.Key)
				}
			}
			raw := make([]byte, 24*len(selected))
			calls := make([]byte, (len(raw)+7)/8)
			for r := 0; r < 24; r++ {
				for col, v := range selected {
					at := r*len(selected) + col
					d := v.Calls[pid-1][r]
					if d >= 0 {
						raw[at] = byte(d)
						calls[at/8] |= 1 << uint(at%8)
					}
				}
			}
			folder, prefix := "geno", "public"
			if side == 1 {
				folder, prefix = "private", "private"
			}
			for name, data := range map[string][]byte{filepath.Join(folder, fmt.Sprintf("block.%d.bin", g)): raw, filepath.Join("joint_stat", fmt.Sprintf("%s.%d.calls", prefix, g)): calls} {
				path := filepath.Join(site, name)
				if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
					t.Fatal(err)
				}
				if err := os.WriteFile(path, data, 0600); err != nil {
					t.Fatal(err)
				}
			}
			if side == 0 {
				publicCounts = append(publicCounts, strings.Join(counts, " "))
			} else {
				privateCounts = append(privateCounts, strings.Join(counts, " "))
			}
		}
		genes = append(genes, gene)
	}
	var geneNames, sizes []string
	for _, g := range genes {
		geneNames = append(geneNames, g.ID)
		sizes = append(sizes, fmt.Sprint(len(g.Public)))
	}
	for name, text := range map[string]string{"cov.txt": strings.Repeat("1\n", 24), "pheno.txt": strings.Repeat("0\n", 24), "../genes.txt": strings.Join(geneNames, "\n") + "\n", "../block_sizes.txt": strings.Join(sizes, "\n") + "\n"} {
		if err := os.WriteFile(filepath.Join(site, name), []byte(text), 0600); err != nil {
			t.Fatal(err)
		}
	}
	var publicRows [][]string
	var privateKeys []string
	for i, gene := range genes {
		keys := gene.Public
		if len(keys) == 0 {
			keys = []string{""}
		}
		for _, key := range keys {
			position := ""
			if key != "" {
				position = strings.SplitN(key, ":", 3)[1]
			}
			publicRows = append(publicRows, []string{strings.SplitN(variants[0].Key, ":", 2)[0], gene.ID, gene.ID, fmt.Sprint(i), key, position})
		}
		privateKeys = append(privateKeys, strings.Join(gene.Private, " "))
	}
	if err := writeTSV(filepath.Join(filepath.Dir(site), "public_variants.tsv"), []string{"chromosome", "gene_id", "gene_symbol", "order_index", "variant_key", "position"}, publicRows); err != nil {
		t.Fatal(err)
	}
	if pid == 2 {
		if err := os.WriteFile(filepath.Join(site, "private_variants.txt"), []byte(strings.Join(privateKeys, "\n")+"\n"), 0600); err != nil {
			t.Fatal(err)
		}
	}
	for name, data := range map[string]string{"cohort_size.txt": "200\n", "allele_counts.txt": strings.Join(publicCounts, "\n") + "\n", "private_allele_counts.txt": strings.Join(privateCounts, "\n") + "\n"} {
		if err := os.WriteFile(filepath.Join(site, name), []byte(data), 0600); err != nil {
			t.Fatal(err)
		}
	}
}

func TestPreparedInputLayout(t *testing.T) {
	variants, groups := fixture()
	for _, name := range []string{"valid", "gene_order", "column_count", "private_rows", "missing_private_keys"} {
		t.Run(name, func(t *testing.T) {
			site := filepath.Join(t.TempDir(), "B")
			writeFixture(t, site, 2, variants, groups)
			var path, text string
			switch name {
			case "gene_order":
				path, text = filepath.Join(site, "..", "genes.txt"), "g1\ng0\ng2\ng3\ng4\n"
			case "column_count":
				path, text = filepath.Join(site, "..", "block_sizes.txt"), "1\n2\n0\n0\n2\n"
			case "private_rows":
				path, text = filepath.Join(site, "private_variants.txt"), "\n"
			case "missing_private_keys":
				if err := os.Remove(filepath.Join(site, "private_variants.txt")); err != nil {
					t.Fatal(err)
				}
			}
			if path != "" {
				if err := os.WriteFile(path, []byte(text), 0600); err != nil {
					t.Fatal(err)
				}
			}
			p, err := prepareInput(site, 2)
			if name != "valid" {
				if err == nil {
					t.Fatal("invalid prepared layout accepted")
				}
				return
			}
			if err != nil {
				t.Fatal(err)
			}
			if len(p.Genes[0].Private) != 0 || len(p.Genes[1].Private) != 1 || len(p.Genes[3].Private) != 0 || len(p.Genes[4].Private) != 0 {
				t.Fatal("empty private rows lost gene alignment")
			}
		})
	}
}

// Independent oracle works on raw participant calls, not on MPC helper functions.
func oracle(variants []fixtureVariant, groups [][]int, maf float64, thresholds, bins []int) []row {
	keep := make([]bool, len(variants))
	flip := make([]bool, len(variants))
	ac := make([][2]int, len(variants))
	for j, v := range variants {
		cohort := v.Cohort[0] + v.Cohort[1]
		keep[j] = float64(min(cohort, 800-cohort))/800 <= maf
		total := 0
		for s := 0; s < 2; s++ {
			for _, x := range v.Calls[s] {
				if x > 0 {
					total += x
				}
			}
		}
		flip[j] = total > 48
		for s := 0; s < 2; s++ {
			if s == 0 && v.Private {
				continue
			}
			for _, x := range v.Calls[s] {
				if x < 0 || !keep[j] {
					continue
				}
				if flip[j] {
					x = 2 - x
				}
				ac[j][s] += x
			}
		}
	}
	counts := make([][3]int, len(groups))
	mac := make([][3]int, len(groups))
	for g, indices := range groups {
		for s := 0; s < 2; s++ {
			for _, j := range indices {
				mac[g][s] += ac[j][s]
			}
			for sample := 0; sample < 24; sample++ {
				has := false
				for _, j := range indices {
					v := variants[j]
					if !keep[j] || (s == 0 && v.Private) {
						continue
					}
					x := v.Calls[s][sample]
					if x < 0 {
						continue
					}
					if flip[j] {
						x = 2 - x
					}
					has = has || x > 0
				}
				if has {
					counts[g][s]++
				}
			}
		}
		counts[g][2] = counts[g][0] + counts[g][1]
		mac[g][2] = mac[g][0] + mac[g][1]
	}
	var out []row
	for _, threshold := range thresholds {
		n := [4]int{}
		for _, c := range counts {
			if c[0] >= threshold {
				n[0]++
			}
			if c[1] >= threshold {
				n[1]++
			}
			if c[0] >= threshold || c[1] >= threshold {
				n[2]++
			}
			if c[2] >= threshold {
				n[3]++
			}
		}
		for i, v := range n {
			out = append(out, row{Panel: "A", Series: []string{"MVP", "AoU", "local_any", "joint"}[i], Lower: fmt.Sprint(threshold), Upper: "NA", Count: v})
		}
	}
	for s := 0; s < 3; s++ {
		for b, lo := range bins {
			hi := math.MaxInt
			upper := "Inf"
			if b+1 < len(bins) {
				hi = bins[b+1]
				upper = fmt.Sprint(hi)
			}
			n := 0
			for g, c := range counts {
				if c[0] >= 20 && c[1] >= 20 && mac[g][s] >= lo && mac[g][s] < hi {
					n++
				}
			}
			out = append(out, row{Panel: "B", Series: []string{"MVP", "AoU", "joint"}[s], Lower: fmt.Sprint(lo), Upper: upper, Count: n})
		}
	}
	n := [3]int{}
	for _, c := range ac {
		switch {
		case c[0] > 0 && c[1] > 0:
			n[0]++
		case c[0] > 0:
			n[1]++
		case c[1] > 0:
			n[2]++
		}
	}
	for i, v := range n {
		out = append(out, row{Panel: "C", Series: []string{"shared", "MVP_only", "AoU_only"}[i], Lower: "NA", Upper: "NA", Count: v})
	}
	return out
}

func TestThreePartyWorkflow(t *testing.T) {
	t.Run("public_and_private", func(t *testing.T) { testWorkflow(t, false, false, 1048576, false, 1) })
	t.Run("only_private_and_empty_genes", func(t *testing.T) { testWorkflow(t, true, false, 1048576, false, 1) })
}

func TestCarrierBatches(t *testing.T) {
	testWorkflow(t, false, false, 87, false, 1)
}

func TestTieOrientation(t *testing.T) {
	testWorkflow(t, false, true, 87, false, 1)
}

func TestReusedAncestrySession(t *testing.T) {
	testWorkflow(t, false, false, 1048576, true, 1)
}

func TestCarrierMessageScaling(t *testing.T) {
	small := testWorkflow(t, false, false, 1048576, false, 1)
	large := testWorkflow(t, false, false, 1048576, false, 20)
	if small != large {
		t.Fatalf("messages grew with genes in one batch: %d -> %d", small, large)
	}
	t.Logf("5 and 100 genes: %d sent messages per chromosome across three MAFs", small)
}

func testWorkflow(t *testing.T, onlyPrivate, ties bool, batchCells int, multiple bool, geneCopies int) uint64 {
	root := t.TempDir()
	variants, groups := fixture()
	if onlyPrivate {
		groups = [][]int{{3}, {3}, {}}
		for j := range variants {
			if j != 3 {
				variants[j].Calls = [2][]int{make([]int, 24), make([]int, 24)}
			}
		}
	}
	if ties {
		// Joint AC=48=N with missing calls: A has 21 ALT carriers but 18 REF carriers.
		for i := 0; i < 24; i++ {
			v := 1
			if i >= 18 {
				v = 2
			}
			if i >= 21 {
				v = -9
			}
			variants[4].Calls[0][i], variants[4].Calls[1][i] = v, 1
			variants[3].Calls[1][i] = 2 // B-private AC also exactly N: keep ALT, with no observed REF.
		}
		variants[4].Cohort = [2]int{24, 24}
		variants[3].Cohort = [2]int{0, 48}
	}
	baseGroups := groups
	for i := 1; i < geneCopies; i++ {
		groups = append(groups, baseGroups...)
	}
	ancestries, chromosomes := []string{"EUR"}, []int{21}
	if multiple {
		ancestries, chromosomes = []string{"EUR", "AFR"}, []int{21, 22}
	}
	for _, anc := range ancestries {
		keys := filepath.Join(root, "keys", anc)
		if err := os.MkdirAll(keys, 0700); err != nil {
			t.Fatal(err)
		}
		for _, name := range []string{"global", "0_1", "0_2", "1_2"} {
			key := make([]byte, 32)
			rand.Read(key)
			if err := os.WriteFile(filepath.Join(keys, "shared_key_"+name+".bin"), key, 0600); err != nil {
				t.Fatal(err)
			}
		}
	}
	var listeners []net.Listener
	var ports []string
	for i := 0; i < 3; i++ {
		l, err := net.Listen("tcp", "127.0.0.1:0")
		if err != nil {
			t.Fatal(err)
		}
		listeners = append(listeners, l)
		ports = append(ports, fmt.Sprint(l.Addr().(*net.TCPAddr).Port))
	}
	for _, l := range listeners {
		l.Close()
	}
	servers := map[string]mpc.Server{"party0": {IpAddr: "127.0.0.1", Ports: map[string]string{"party1": ports[0], "party2": ports[1]}}, "party1": {IpAddr: "127.0.0.1", Ports: map[string]string{"party2": ports[2]}}, "party2": {IpAddr: "127.0.0.1", Ports: map[string]string{}}}
	opt := options{Thresholds: []int{20, 30, 50}, Bins: []int{20, 40, 80}, BatchCells: batchCells}
	configs := make([]*workflow.Config, 3)
	for pid := 0; pid < 3; pid++ {
		dir := filepath.Join(root, fmt.Sprintf("party%d", pid))
		side := []string{"helper", "A", "B"}[pid]
		if pid > 0 {
			for _, anc := range ancestries {
				for _, chr := range chromosomes {
					copyVariants := append([]fixtureVariant(nil), variants...)
					for i := range copyVariants {
						copyVariants[i].Key = fmt.Sprintf("%d:%d:A:G", chr, i+1)
					}
					writeFixture(t, filepath.Join(dir, "prepared", anc, fmt.Sprintf("chr%d", chr), side), pid, copyVariants, groups)
				}
			}
		}
		configs[pid] = &workflow.Config{RunDir: dir, SharedKeysPath: filepath.Join(root, "keys"), Ancestries: ancestries, Chromosomes: chromosomes, GlobalMaxMAF: []float64{0.01, 0.1, 0.5}, GenotypeDirectory: filepath.Join(dir, "prepared", "{ancestry}", "chr{chromosome}", side, "geno"), BindingIP: "127.0.0.1", Servers: servers, DataBits: 70, FractionalBits: 40}
		if multiple {
			// Shipped 1KG uses 60/30; AoU/MVP fixtures above use 70/40.
			configs[pid].DataBits, configs[pid].FractionalBits = 60, 30
		}
	}
	errors := make([]error, 3)
	var wg sync.WaitGroup
	for attempt := 0; attempt < 2; attempt++ {
		for pid := 0; pid < 3; pid++ {
			wg.Add(1)
			go func(pid int) { defer wg.Done(); errors[pid] = run(configs[pid], pid, opt) }(pid)
			time.Sleep(50 * time.Millisecond)
		}
		wg.Wait()
		for pid, err := range errors {
			if err != nil {
				t.Fatalf("attempt%d party%d: %v", attempt, pid, err)
			}
		}
	}
	a, err := os.ReadFile(filepath.Join(configs[1].RunDir, "resource_stats", "resource_stats.tsv"))
	if err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile(filepath.Join(configs[2].RunDir, "resource_stats", "resource_stats.tsv"))
	if err != nil {
		t.Fatal(err)
	}
	if string(a) != string(b) {
		t.Fatal("site TSVs differ")
	}
	var want []row
	for _, anc := range ancestries {
		for _, chr := range chromosomes {
			for _, maf := range configs[1].GlobalMaxMAF {
				for _, r := range oracle(variants, groups, maf, opt.Thresholds, opt.Bins) {
					r.Ancestry, r.Chromosome, r.MAF, r.Annotation = anc, chr, maf, "HC_pLoF"
					want = append(want, r)
				}
			}
		}
	}
	expected := filepath.Join(root, "oracle.tsv")
	writeResults(expected, want)
	oracleData, _ := os.ReadFile(expected)
	if string(a) != string(oracleData) {
		t.Fatalf("secure result differs from independent participant-level oracle\nGOT:\n%s\nWANT:\n%s", a, oracleData)
	}
	if _, err = os.Stat(filepath.Join(configs[0].RunDir, "resource_stats", "resource_stats.tsv")); !os.IsNotExist(err) {
		t.Fatal("helper received output")
	}
	for pid := 0; pid < 3; pid++ {
		if _, err = os.Stat(filepath.Join(configs[pid].RunDir, "metrics")); !os.IsNotExist(err) {
			t.Fatal("joint-stat wrote main metrics")
		}
	}
	for pid := 0; pid < 3; pid++ {
		if _, err := os.Stat(filepath.Join(configs[pid].RunDir, "cache", "joint_stat")); !os.IsNotExist(err) {
			t.Fatal("unnecessary count cache created")
		}
	}
	if err := os.WriteFile(filepath.Join(sitePath(configs[1], "EUR", 21), "joint_stat", "public.0.calls"), []byte("invalid size"), 0600); err != nil {
		t.Fatal(err)
	}
	if _, err := prepareInput(sitePath(configs[1], "EUR", 21), 1); err == nil {
		t.Fatal("invalid callability dimensions accepted")
	}
	file, err := os.Open(filepath.Join(configs[1].RunDir, "resource_stats", "metrics_party1.tsv"))
	if err != nil {
		t.Fatal(err)
	}
	defer file.Close()
	reader := csv.NewReader(file)
	reader.Comma = '\t'
	records, err := reader.ReadAll()
	if err != nil {
		t.Fatal(err)
	}
	column := -1
	for i, name := range records[0] {
		if name == "sent_messages" {
			column = i
		}
	}
	if column < 0 {
		t.Fatal("missing message metric")
	}
	var messages uint64
	for _, record := range records[1:] {
		if record[0] == "secure_aggregate" {
			n, err := strconv.ParseUint(record[column], 10, 64)
			if err != nil {
				t.Fatal(err)
			}
			messages += n
		}
	}
	return messages

}
