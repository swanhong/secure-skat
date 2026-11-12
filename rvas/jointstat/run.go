package jointstat

import (
	"crypto/hmac"
	"crypto/rand"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"reflect"
	"runtime"
	"strconv"
	"strings"
	"time"

	core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rvas/workflow"
)

const Version = "1"
const annotationClass = "HC_pLoF"

type options struct {
	Thresholds, Bins []int
	BatchCells       int
}
type measurement struct {
	Stage, Ancestry                string
	Chromosome                     int
	MAF                            float64
	Seconds                        float64
	Sent, Received                 uint64
	SentMessages, ReceivedMessages uint64
}

func integers(text string, minimum int) ([]int, error) {
	var values []int
	for _, s := range strings.Split(text, ",") {
		v, err := strconv.Atoi(s)
		if err != nil || v < minimum || v > 1000000000 || (len(values) > 0 && v <= values[len(values)-1]) {
			return nil, fmt.Errorf("expected increasing integers >=%d: %s", minimum, text)
		}
		values = append(values, v)
	}
	return values, nil
}

func Main(args []string) error {
	flags := flag.NewFlagSet("joint-stat", flag.ContinueOnError)
	dir := flags.String("config", "", "party configuration directory")
	pid := flags.Int("party", -1, "0=helper, 1=MVP, 2=AoU")
	ts := flags.String("carrier-thresholds", "20,30,50,100,200,500,1000", "public carrier thresholds")
	bs := flags.String("cmac-bins", "20,40,80,160,320,640,1280,2560,5120", "lower bin edges; final bin is unbounded")
	cells := flags.Int("batch-cells", 1048576, "target input field elements per carrier batch (at least one row)")
	if err := flags.Parse(args); err != nil {
		if err == flag.ErrHelp {
			return nil
		}
		return err
	}
	if *dir == "" || *pid < 0 || *pid > 2 || flags.NArg() != 0 || *cells < 1 {
		return fmt.Errorf("config, party 0/1/2 and positive batch-cells are required; no positional arguments")
	}
	thresholds, err := integers(*ts, 20)
	if err != nil {
		return err
	}
	bins, err := integers(*bs, 20)
	if err != nil {
		return err
	}
	if bins[0] != 20 {
		return fmt.Errorf("first cMAC edge must be 20")
	}
	config, err := workflow.LoadPartyConfig(*dir, *pid)
	if err != nil {
		return err
	}
	if config.RunDir == "" || config.SharedKeysPath == "" || len(config.Ancestries) == 0 || len(config.Chromosomes) == 0 || len(config.GlobalMaxMAF) == 0 {
		return fmt.Errorf("run_dir, shared_keys_path, ancestries, chromosomes and global_max_maf are required")
	}
	if config.DataBits < 64 || config.DataBits > 120 || config.FractionalBits < 1 || config.FractionalBits > 40 || config.DataBits-config.FractionalBits < 24 {
		return fmt.Errorf("joint-stat requires data_bits 64..120, fractional_bits 1..40 and >=24 integer bits")
	}
	if config.LocalNumThreads > 0 {
		runtime.GOMAXPROCS(config.LocalNumThreads)
	}
	return run(config, *pid, options{Thresholds: thresholds, Bins: bins, BatchCells: *cells})
}

type hello struct {
	Spec, LayoutDigest, Nonce string
	Samples, Cohort           int
	Shape                     layout
}
type agreement struct {
	Error, Nonce string
	Dims         dimensions
}

func sendJSON(net *mpc.Network, to int, value interface{}) {
	raw, err := json.Marshal(value)
	if err != nil {
		panic(err)
	}
	v := make([]uint64, len(raw))
	for i, b := range raw {
		v[i] = uint64(b)
	}
	net.SendInt(len(v), to)
	net.SendIntVector(v, to)
}
func receiveJSON(net *mpc.Network, from int, value interface{}) error {
	n := net.ReceiveInt(from)
	if n < 0 || n > 64<<20 {
		return fmt.Errorf("oversized public handshake")
	}
	v := net.ReceiveIntVector(n, from)
	raw := make([]byte, n)
	for i, b := range v {
		if b > 255 {
			return fmt.Errorf("invalid handshake byte")
		}
		raw[i] = byte(b)
	}
	return json.Unmarshal(raw, value)
}
func hashJSON(v interface{}) string {
	raw, _ := json.Marshal(v)
	h := sha256.Sum256(raw)
	return hex.EncodeToString(h[:])
}

func agree(net *mpc.Network, pid int, spec string, p *prepared) (agreement, error) {
	nonce := make([]byte, 32)
	if _, err := rand.Read(nonce); err != nil {
		return agreement{}, err
	}
	h := hello{Spec: spec, Nonce: hex.EncodeToString(nonce)}
	if pid > 0 {
		l := publicLayout(p)
		h.LayoutDigest = hashJSON(l)
		h.Samples = p.Samples
		h.Cohort = p.CohortSize
		// Helper learns dimensions and duplicate-column routing, never private keys
		// or even public variant/gene identities from this program.
		h.Shape = layout{Genes: make([]string, len(l.Genes)), Keys: make([]string, len(l.Keys)), Mapping: l.Mapping}
	}
	var a agreement
	if pid == 0 {
		peers := []hello{h, {}, {}}
		for s := 1; s <= 2; s++ {
			if err := receiveJSON(net, s, &peers[s]); err != nil {
				return a, err
			}
		}
		if peers[1].Spec != spec || peers[2].Spec != spec || peers[1].LayoutDigest != peers[2].LayoutDigest || !reflect.DeepEqual(peers[1].Shape, peers[2].Shape) {
			a.Error = "public settings or variant/gene alignment differ across parties"
		}
		for s := 1; s <= 2; s++ {
			if peers[s].Samples <= 0 || peers[s].Cohort < peers[s].Samples || peers[s].Cohort > 10000000 {
				a.Error = "invalid public sample dimensions"
			}
		}
		a.Nonce = hashJSON([]string{peers[0].Nonce, peers[1].Nonce, peers[2].Nonce})
		a.Dims = dimensions{peers[1].Samples, peers[2].Samples, peers[1].Cohort + peers[2].Cohort, peers[1].Shape}
		for s := 1; s <= 2; s++ {
			sendJSON(net, s, a)
		}
	} else {
		sendJSON(net, 0, h)
		if err := receiveJSON(net, 0, &a); err != nil {
			return a, err
		}
	}
	if a.Error != "" {
		return a, fmt.Errorf("%s", a.Error)
	}
	return a, nil
}

// Derive a fresh, domain-separated stream for every session. Reusing SKAT's
// initial PRG stream or restarting a joint-stat run must never reuse masks.
func freshPRG(base, nonce string, pid int, networks []*mpc.Network) error {
	tmp, err := os.MkdirTemp("", "joint-stat-prg-")
	if err != nil {
		return err
	}
	defer os.RemoveAll(tmp)
	names := [][]string{{"global", "0_1", "0_2"}, {"global", "0_1", "1_2"}, {"global", "0_2", "1_2"}}[pid]
	for _, name := range names {
		filename := "shared_key_" + name + ".bin"
		key, err := os.ReadFile(filepath.Join(base, filename))
		if err != nil {
			return err
		}
		if len(key) != 32 {
			return fmt.Errorf("invalid PRG key")
		}
		h := hmac.New(sha256.New, key)
		h.Write([]byte("joint-stat/v" + Version + "/" + nonce + "/" + name))
		if err = os.WriteFile(filepath.Join(tmp, filename), h.Sum(nil), 0600); err != nil {
			return err
		}
	}
	mpc.InitializeParallelPRG(tmp, networks, pid, 3)
	return nil
}

func run(config *workflow.Config, pid int, opt options) error {
	totalStart := time.Now()
	out := filepath.Join(config.RunDir, "resource_stats")
	if err := os.MkdirAll(out, 0700); err != nil {
		return err
	}
	marker := filepath.Join(out, fmt.Sprintf("_SUCCESS_party%d", pid))
	if err := os.Remove(marker); err != nil && !os.IsNotExist(err) {
		return err
	}
	var results []row
	var metrics []measurement
	add := func(stage, anc string, chr int, start time.Time, before mpc.CommunicationStats, networks []*mpc.Network) {
		var delta mpc.CommunicationStats
		if len(networks) > 0 {
			delta = mpc.ParallelNetworks(networks).GetCommunicationStats().Sub(before)
		}
		m := measurement{Stage: stage, Ancestry: anc, Chromosome: chr, Seconds: time.Since(start).Seconds(), Sent: delta.SentBytes, Received: delta.ReceivedBytes, SentMessages: delta.SentMessages, ReceivedMessages: delta.ReceivedMessages}
		metrics = append(metrics, m)
		fmt.Printf("[joint-stat party%d %s chr%d] %s: %.3fs, sent %d messages, received %d messages\n", pid, anc, chr, stage, m.Seconds, m.SentMessages, m.ReceivedMessages)
	}
	for _, anc := range config.Ancestries {
		err := func() error {
			var networks []*mpc.Network
			var obj *mpc.MPC
			defer func() {
				if len(networks) > 0 {
					networks[0].CloseAll()
				}
			}()
			for _, chr := range config.Chromosomes {
				start := time.Now()
				var p *prepared
				site := ""
				if pid > 0 {
					site = sitePath(config, anc, chr)
					var err error
					p, err = prepareInput(site, pid)
					if err != nil {
						return err
					}
				}
				add("local_prepare", anc, chr, start, mpc.CommunicationStats{}, nil)
				start = time.Now()
				keyPath := filepath.Join(config.SharedKeysPath, anc)
				if len(networks) == 0 {
					networks = mpc.InitCommunication(config.BindingIP, config.Servers, pid, 3, 1, keyPath)
					networks[0].EnableLogging()
				}
				before := mpc.ParallelNetworks(networks).GetCommunicationStats()
				spec := hashJSON(struct {
					Version, Ancestry      string
					Chromosome             int
					Anc                    []string
					Chr                    []int
					MAF                    []float64
					Options                options
					Data, Frac, Phenotypes int
				}{Version, anc, chr, config.Ancestries, config.Chromosomes, config.GlobalMaxMAF, opt, config.DataBits, config.FractionalBits, len(config.PhenotypeColumns)})
				// Each chromosome has its own public layout, even on a reused connection.
				a, err := agree(networks[0], pid, spec, p)
				if err != nil {
					return err
				}
				if obj == nil {
					if err = freshPRG(keyPath, a.Nonce, pid, networks); err != nil {
						return err
					}
					obj = mpc.InitParallelMPCEnv(networks, core.LElem256Zero, config.DataBits, config.FractionalBits)[0]
					obj.SetHubPid(1)
				}
				add("setup", anc, chr, start, before, networks)
				start = time.Now()
				before = mpc.ParallelNetworks(networks).GetCommunicationStats()
				calculations, err := compute(obj, p, site, a.Dims, config.GlobalMaxMAF, opt.Thresholds, opt.Bins, opt.BatchCells)
				if err != nil {
					return err
				}
				add("secure_aggregate", anc, chr, start, before, networks)
				start = time.Now()
				before = mpc.ParallelNetworks(networks).GetCommunicationStats()
				for m, c := range calculations {
					// RevealSymVec reconstructs to parties 1 and 2 only. Party 0 returns zeros.
					plain := obj.RevealSymVec(c.Shares)
					if pid == 0 {
						continue
					}
					counts := plain.ToInt()
					for i, r := range c.Labels {
						if counts[i] < 0 || counts[i] > 1000000000 {
							return fmt.Errorf("invalid released count")
						}
						r.Ancestry, r.Chromosome, r.MAF, r.Annotation, r.Count = anc, chr, config.GlobalMaxMAF[m], annotationClass, counts[i]
						results = append(results, r)
					}
				}
				add("release", anc, chr, start, before, networks)
			}
			return nil
		}()
		if err != nil {
			return err
		}
	}
	start := time.Now()
	if pid > 0 {
		if err := writeResults(filepath.Join(out, "resource_stats.tsv"), results); err != nil {
			return err
		}
		metadata := map[string]interface{}{"schema_version": Version, "annotation": annotationClass, "carrier_thresholds": opt.Thresholds, "cmac_bins": opt.Bins, "batch_cells": opt.BatchCells, "ancestries": config.Ancestries, "chromosomes": config.Chromosomes, "maf_thresholds": config.GlobalMaxMAF}
		if err := writeJSON(filepath.Join(out, "metadata.json"), metadata); err != nil {
			return err
		}
	}
	add("write", "", 0, start, mpc.CommunicationStats{}, nil)
	total := measurement{Stage: "figure_total", Seconds: time.Since(totalStart).Seconds()}
	for _, m := range metrics {
		total.Sent += m.Sent
		total.Received += m.Received
		total.SentMessages += m.SentMessages
		total.ReceivedMessages += m.ReceivedMessages
	}
	metrics = append(metrics, total)
	if err := writeMetrics(filepath.Join(out, fmt.Sprintf("metrics_party%d.tsv", pid)), metrics); err != nil {
		return err
	}
	return os.WriteFile(marker, []byte("complete\n"), 0600)
}
