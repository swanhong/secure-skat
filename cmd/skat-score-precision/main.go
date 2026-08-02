// Command skat-score-precision screens the numerical stability of replacing the
// public SS score with across-gene packed CKKS. It runs entirely in one process
// over an existing fed_prep output directory: no MPC, networking, collective
// decryption, PLINK, or raw PGEN access is used.
//
// The plaintext-only mode is a necessary cancellation screen, not a CKKS
// acceptance test. The default mode additionally executes the actual Lattigo
// plaintext-times-ciphertext score circuit and compares it with a residual-first
// stable plaintext reference.
package main

import (
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"
)

type options struct {
	Root       string
	CKKS       string
	FracBits   int
	MaxGenes   int
	CSV        string
	Top        int
	PlainOnly  bool
	PublicOnly bool
	Threads    int
	Progress   int
}

func defaultRoot() string {
	if v := os.Getenv("FED_OUT"); v != "" {
		return v
	}
	h, err := os.UserHomeDir()
	if err != nil {
		return "fed_prep_out"
	}
	return filepath.Join(h, "fed_prep_out")
}

func expandHome(path string) string {
	if path == "~" || strings.HasPrefix(path, "~/") {
		if h, err := os.UserHomeDir(); err == nil {
			if path == "~" {
				return h
			}
			return filepath.Join(h, path[2:])
		}
	}
	return path
}

func parseOptions() options {
	var o options
	flag.StringVar(&o.Root, "out", defaultRoot(), "fed_prep output directory (default $FED_OUT or ~/fed_prep_out)")
	flag.StringVar(&o.CKKS, "ckks", "auto", "CKKS parameters: auto, PN13QP218, PN14QP438, PN15QP880, or PN16QP1761")
	flag.IntVar(&o.FracBits, "frac-bits", -1, "quantize beta to this many fractional bits; -1 reads config")
	flag.IntVar(&o.MaxGenes, "max-genes", 0, "evenly sample this many genes; 0 means all")
	flag.StringVar(&o.CSV, "csv", "auto", "per-gene aggregate CSV path; auto writes under --out")
	flag.IntVar(&o.Top, "top", 15, "number of worst aggregate genes to print")
	flag.BoolVar(&o.PlainOnly, "plain-only", false, "run only the fast cancellation screen, without CKKS")
	flag.BoolVar(&o.PublicOnly, "public-only", false, "check only the public-list score proposed for SS-to-CKKS migration")
	flag.IntVar(&o.Threads, "threads", 1, "local CKKS object-pool size (1 is sufficient for this sequential diagnostic)")
	flag.IntVar(&o.Progress, "progress-every", 50, "print progress every this many selected genes; 0 disables")
	flag.Parse()
	o.Root = expandHome(o.Root)
	if o.MaxGenes < 0 {
		fatalf("--max-genes must be >=0")
	}
	if o.Threads <= 0 {
		fatalf("--threads must be >0")
	}
	return o
}

func fatalf(format string, args ...interface{}) {
	fmt.Fprintf(os.Stderr, "error: "+format+"\n", args...)
	os.Exit(2)
}

func geneName(m manifest, idx int) string {
	if len(m.GeneNames) == m.NGenes && m.GeneNames[idx] != "" {
		return m.GeneNames[idx]
	}
	return fmt.Sprintf("gene_%d", idx)
}

func main() {
	o := parseOptions()
	started := time.Now()
	m, err := loadManifest(o.Root)
	if err != nil {
		fatalf("load manifest: %v", err)
	}
	cfg, err := loadGlobalConfig(o.Root)
	if err != nil {
		fatalf("load config: %v", err)
	}
	if o.CKKS == "auto" {
		o.CKKS = cfg.CKKSParams
	}
	if o.FracBits < 0 {
		o.FracBits = cfg.MPCFracBits
	}
	if o.FracBits < 0 || o.FracBits > 52 {
		fatalf("--frac-bits=%d is outside the safe float64 diagnostic range [0,52]", o.FracBits)
	}
	if o.CSV == "auto" {
		tag := strings.ToLower(o.CKKS)
		if o.PlainOnly {
			tag = "plain"
		}
		if o.PublicOnly {
			tag += "_public"
		}
		o.CSV = filepath.Join(o.Root, "score_precision_"+tag+".csv")
	} else {
		o.CSV = expandHome(o.CSV)
	}

	a, err := loadCohort(o.Root, "A", cfg.BinaryPheno, cfg.NumCovs)
	if err != nil {
		fatalf("load cohort A: %v", err)
	}
	b, err := loadCohort(o.Root, "B", cfg.BinaryPheno, cfg.NumCovs)
	if err != nil {
		fatalf("load cohort B: %v", err)
	}
	if a.C != b.C {
		fatalf("cohort covariate widths differ: A=%d B=%d", a.C, b.C)
	}
	beta, err := fitNull(&a, &b, o.FracBits)
	if err != nil {
		fatalf("fit null model: %v", err)
	}
	selected := chooseGenes(m.NGenes, o.MaxGenes)
	selectedSet := make(map[int]bool, len(selected))
	for _, g := range selected {
		selectedSet[g] = true
	}

	publicMetrics := make([]categoryMetrics, m.NGenes)
	privateMetrics := make([]categoryMetrics, m.NGenes)
	var engine *heEngine
	var publicChunk, privateChunk *scoreChunker
	if !o.PlainOnly {
		engine, err = newHEEngine(o.CKKS, beta, uint(cfg.MPCFieldSize), o.Threads)
		if err != nil {
			fatalf("initialize CKKS: %v", err)
		}
		publicChunk = newScoreChunker(engine, publicMetrics, true)
		if !o.PublicOnly {
			privateChunk = newScoreChunker(engine, privateMetrics, false)
		}
	}

	fmt.Println("=== SKAT packed-score precision diagnostic ===")
	fmt.Println("  input=existing Workbench fed_prep output (path kept out of summary)")
	fmt.Printf("  cohorts: A=%d B=%d total=%d; covariates including intercept=%d\n", a.N, b.N, a.N+b.N, a.C)
	fmt.Printf("  genes=%d/%d (even deterministic sample when limited); beta fractional bits=%d\n", len(selected), m.NGenes, o.FracBits)
	if o.PlainOnly {
		fmt.Println("  mode=plain cancellation screen only (cannot by itself accept CKKS)")
	} else {
		fmt.Printf("  mode=single-process packed CKKS score circuit; params=%s slots=%d threads=%d\n", o.CKKS, engine.Slots, o.Threads)
	}
	if o.PublicOnly {
		fmt.Println("  scope=public-list score only (B-private score skipped)")
	}

	totalN := a.N + b.N
	for pos, g := range selected {
		pubA, err := contractBlock(filepath.Join(o.Root, "A", fmt.Sprintf("geno.%d.bin", g)), &a, m.PublicM[g])
		if err != nil {
			fatalf("gene %d A public block: %v", g, err)
		}
		pubB, err := contractBlock(filepath.Join(o.Root, "B", fmt.Sprintf("geno.%d.bin", g)), &b, m.PublicM[g])
		if err != nil {
			fatalf("gene %d B public block: %v", g, err)
		}
		pubRecords, pubKappa := recordsForGene(g, pubA, pubB, beta, totalN, o.FracBits)
		for _, r := range pubRecords {
			publicMetrics[g].addPlain(r)
			if publicChunk != nil {
				if err := publicChunk.Add(r); err != nil {
					fatalf("public CKKS tile: %v", err)
				}
			}
		}
		publicMetrics[g].setKappas(pubKappa)

		if !o.PublicOnly {
			privB, err := contractBlock(filepath.Join(o.Root, "B", fmt.Sprintf("priv.%d.bin", g)), &b, m.PrivateM[g])
			if err != nil {
				fatalf("gene %d B private block: %v", g, err)
			}
			privRecords, privKappa := recordsForGene(g, privB, emptyContraction(a.C, m.PrivateM[g]), beta, totalN, o.FracBits)
			for _, r := range privRecords {
				privateMetrics[g].addPlain(r)
				if privateChunk != nil {
					if err := privateChunk.Add(r); err != nil {
						fatalf("private CKKS tile: %v", err)
					}
				}
			}
			privateMetrics[g].setKappas(privKappa)
		}

		if o.Progress > 0 && ((pos+1)%o.Progress == 0 || pos+1 == len(selected)) {
			fmt.Printf("  processed %d/%d genes (%.1fs)\n", pos+1, len(selected), time.Since(started).Seconds())
		}
	}
	if publicChunk != nil {
		if err := publicChunk.Flush(); err != nil {
			fatalf("flush public CKKS tile: %v", err)
		}
		if privateChunk != nil {
			if err := privateChunk.Flush(); err != nil {
				fatalf("flush private CKKS tile: %v", err)
			}
		}
	}

	rows := make([]geneResult, 0, len(selected))
	for g := 0; g < m.NGenes; g++ {
		if !selectedSet[g] {
			continue
		}
		rows = append(rows, geneResult{Index: g, Name: geneName(m, g), Public: publicMetrics[g], Private: privateMetrics[g]})
	}
	printReport(rows, !o.PlainOnly, o.Top)
	if o.PlainOnly {
		fmt.Println("\nSCREEN_STATUS=PLAINTEXT_REJECTION_SCREEN_COMPLETE")
	} else {
		invalidGenes := 0
		for _, g := range rows {
			metrics := g.Public
			if !o.PublicOnly {
				metrics = g.Total()
			}
			if metrics.HEInvalid > 0 || (metrics.Variants > 0 && !metrics.HasHE) {
				invalidGenes++
			}
		}
		if invalidGenes > 0 {
			fmt.Printf("\nSCREEN_STATUS=INVALID_NONFINITE genes=%d\n", invalidGenes)
		} else {
			fmt.Println("\nSCREEN_STATUS=DIRECT_CKKS_REJECTION_SCREEN_COMPLETE")
		}
	}
	if err := os.MkdirAll(filepath.Dir(o.CSV), 0o755); err != nil {
		fatalf("create CSV directory: %v", err)
	}
	if err := writeCSV(o.CSV, rows, !o.PlainOnly); err != nil {
		fatalf("write CSV: %v", err)
	}
	fmt.Printf("\nCSV written under --out: %s\n", filepath.Base(o.CSV))
	fmt.Printf("TOTAL: %.2fs\n", time.Since(started).Seconds())
	if o.PlainOnly {
		fmt.Println("Interpretation: this mode can reject a poorly conditioned CKKS score, but passing it is not sufficient; run again without --plain-only for the actual local CKKS differential.")
	} else {
		fmt.Println("Interpretation: this is a direct-Enc(beta) score-circuit differential and a rejection screen, not a production acceptance gate. It does not simulate MPC fixed-point null-solve error, collective SS<->CKKS conversion, encrypted weights, moments, or the final p-value circuit.")
	}
}
