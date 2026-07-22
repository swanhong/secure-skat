package gwas

import (
	"fmt"
	"math"
	"sort"
	"time"

	"github.com/hhcho/sfgwas/mpc"
)

// fedStageDef names non-overlapping secure-SKAT leaves. The math-oriented names
// describe the secure action using a local contraction; they do not imply that
// G^T G, G^T X, or G^T y is transmitted as a plaintext matrix.
type fedStageDef struct {
	name   string
	parent string
}

var fedStageDefs = []fedStageDef{
	{"load_private_only", "secure_run"},
	{"assoc_init", "secure_run"},
	{"null_local_xtx_xty_yty", "null_model"},
	{"null_aggregate_xtx", "null_model"},
	{"null_aggregate_xty", "null_model"},
	{"null_aggregate_yty", "null_model"},
	{"null_solve", "null_model"},
	{"null_beta_pack", "null_model"},
	{"null_rss", "null_model"},
	{"null_other", "null_model"},
	{"pre_block_setup", "secure_compute"},
	{"gene_shape_sync", "genes"},
	{"gene_local_public_gtg_gtx_gty", "genes"},
	{"gene_public_score_gtx_gty", "genes"},
	{"gene_public_weight_dosage_maf", "genes"},
	{"gene_public_stat_share", "genes"},
	{"gene_local_private_gtg_gtx_gty", "genes"},
	{"gene_private_stat_share", "genes"},
	{"gene_burden_public_gtg", "genes"},
	{"gene_burden_public_gtx", "genes"},
	{"gene_burden_private_cross", "genes"},
	{"gene_burden_projection", "genes"},
	{"gene_moments_setup_gtx", "genes"},
	{"gene_moments_public_trace_gtg_gtx", "genes"},
	{"gene_moments_private_cross", "genes"},
	{"gene_other", "genes"},
	{"finalize_scale", "post_block_finalize"},
	{"finalize_burden_pvalue", "post_block_finalize"},
	{"finalize_skat_pvalue", "post_block_finalize"},
	{"finalize_output_pack", "post_block_finalize"},
	{"finalize_other", "post_block_finalize"},
	{"decrypt_outputs", "secure_run"},
	{"run_other", "secure_run"},
}

type fedMetricMark struct {
	at   time.Time
	comm mpc.CommunicationStats
}

type fedStageValue struct {
	duration time.Duration
	comm     mpc.CommunicationStats
	count    int
}

type fedRunMetrics struct {
	networks   mpc.ParallelNetworks
	mode       string
	pid        int
	runStarted time.Time
	runStart   mpc.CommunicationStats
	setupComm  mpc.CommunicationStats
	initComm   mpc.CommunicationStats
	values     map[string]*fedStageValue
}

func newFedRunMetrics(networks mpc.ParallelNetworks, mode string, setupComm mpc.CommunicationStats) *fedRunMetrics {
	pid := -1
	if len(networks) > 0 {
		pid = networks[0].GetPid()
	}
	preRun := networks.GetCommunicationStats()
	m := &fedRunMetrics{
		networks:   networks,
		mode:       mode,
		pid:        pid,
		runStarted: time.Now(),
		runStart:   preRun,
		setupComm:  setupComm,
		initComm:   preRun.Sub(setupComm),
		values:     make(map[string]*fedStageValue, len(fedStageDefs)),
	}
	for _, def := range fedStageDefs {
		m.values[def.name] = &fedStageValue{}
	}
	return m
}

func communicationAdd(a, b mpc.CommunicationStats) mpc.CommunicationStats {
	return mpc.CommunicationStats{
		SentBytes:        a.SentBytes + b.SentBytes,
		ReceivedBytes:    a.ReceivedBytes + b.ReceivedBytes,
		SentMessages:     a.SentMessages + b.SentMessages,
		ReceivedMessages: a.ReceivedMessages + b.ReceivedMessages,
	}
}

func (m *fedRunMetrics) mark() fedMetricMark {
	if m == nil {
		return fedMetricMark{}
	}
	return fedMetricMark{at: time.Now(), comm: m.networks.GetCommunicationStats()}
}

func (m *fedRunMetrics) end(stage string, mark fedMetricMark) {
	if m == nil {
		return
	}
	v, ok := m.values[stage]
	if !ok {
		panic("unknown secure-SKAT metric stage: " + stage)
	}
	v.duration += time.Since(mark.at)
	v.comm = communicationAdd(v.comm, m.networks.GetCommunicationStats().Sub(mark.comm))
	v.count++
}

func (m *fedRunMetrics) addDuration(stage string, duration time.Duration) {
	m.addDurationCount(stage, duration, 1)
}

func (m *fedRunMetrics) addDurationCount(stage string, duration time.Duration, count int) {
	if m == nil {
		return
	}
	v, ok := m.values[stage]
	if !ok {
		panic("unknown secure-SKAT metric stage: " + stage)
	}
	v.duration += duration
	v.count += count
}

func (m *fedRunMetrics) stageDuration(stage string) time.Duration {
	if m == nil || m.values[stage] == nil {
		return 0
	}
	return m.values[stage].duration
}

func (m *fedRunMetrics) parentLeafDuration(parent, except string) time.Duration {
	if m == nil {
		return 0
	}
	var total time.Duration
	for _, def := range fedStageDefs {
		if def.parent == parent && def.name != except {
			total += m.values[def.name].duration
		}
	}
	return total
}

func (m *fedRunMetrics) finishRun() time.Duration {
	if m == nil {
		return 0
	}
	runDuration := time.Since(m.runStarted)
	runTotal := m.networks.GetCommunicationStats().Sub(m.runStart)
	var classified mpc.CommunicationStats
	for _, def := range fedStageDefs {
		if def.name != "run_other" {
			classified = communicationAdd(classified, m.values[def.name].comm)
		}
	}
	m.values["run_other"].comm = runTotal.Sub(classified)
	classifiedDuration := time.Duration(0)
	for _, def := range fedStageDefs {
		if def.name != "run_other" {
			classifiedDuration += m.values[def.name].duration
		}
	}
	m.values["run_other"].duration = nonNegativeDuration(runDuration - classifiedDuration)
	m.values["run_other"].count = 1
	return runDuration
}

func printTimingRecord(mode string, pid int, stage, parent, kind string, duration time.Duration,
	count int, minMs, q1Ms, meanMs, q3Ms, maxMs int64) {
	status := "done"
	if count == 0 {
		status = "skipped"
	}
	fmt.Printf("[timing] scope=secure mode=%s party=%d stage=%s parent=%s kind=%s status=%s milliseconds=%.3f count=%d",
		mode, pid, stage, parent, kind, status, float64(duration)/float64(time.Millisecond), count)
	if kind == "distribution" {
		fmt.Printf(" min_ms=%d q1_ms=%d mean_ms=%d q3_ms=%d max_ms=%d", minMs, q1Ms, meanMs, q3Ms, maxMs)
	}
	fmt.Println()
}

func printCommunicationStage(mode string, pid int, stage, parent, kind string, stats mpc.CommunicationStats) {
	fmt.Printf("[communication-stage] scope=skat_fed_total mode=%s party=%d stage=%s parent=%s kind=%s sent_bytes=%d received_bytes=%d sent_messages=%d received_messages=%d\n",
		mode, pid, stage, parent, kind, stats.SentBytes, stats.ReceivedBytes, stats.SentMessages, stats.ReceivedMessages)
}

func (m *fedRunMetrics) printRecords() {
	if m == nil {
		return
	}
	// These two leaves partition all traffic before runFederatedPrivate. Local
	// file loading normally makes initialization_other exactly zero bytes.
	printCommunicationStage(m.mode, m.pid, "collective_key_setup", "protocol_total", "leaf", m.setupComm)
	printCommunicationStage(m.mode, m.pid, "initialization_other", "protocol_total", "leaf", m.initComm)
	for _, def := range fedStageDefs {
		v := m.values[def.name]
		printTimingRecord(m.mode, m.pid, def.name, def.parent, "leaf", v.duration, v.count, 0, 0, 0, 0, 0)
		printCommunicationStage(m.mode, m.pid, def.name, def.parent, "leaf", v.comm)
	}
}

func (ast *AssocTest) metricMark() fedMetricMark {
	if ast == nil {
		return fedMetricMark{}
	}
	return ast.fedMetrics.mark()
}

func (ast *AssocTest) metricEnd(stage string, mark fedMetricMark) {
	if ast != nil {
		ast.fedMetrics.end(stage, mark)
	}
}

func nonNegativeDuration(d time.Duration) time.Duration {
	if d < 0 {
		return 0
	}
	return d
}

func blockTimingDistribution(secs []float64) (count int, minMs, q1Ms, meanMs, q3Ms, maxMs int64) {
	if len(secs) == 0 {
		return
	}
	sorted := append([]float64(nil), secs...)
	sort.Float64s(sorted)
	var sum float64
	for _, value := range sorted {
		sum += value
	}
	q := func(p float64) float64 { return sorted[int(p*float64(len(sorted)-1)+0.5)] }
	toMs := func(seconds float64) int64 { return int64(math.Round(seconds * 1000)) }
	return len(sorted), toMs(sorted[0]), toMs(q(0.25)), toMs(sum / float64(len(sorted))),
		toMs(q(0.75)), toMs(sorted[len(sorted)-1])
}

// printFedTimingRecords emits the same timing tree as machine-readable records.
// Inclusive rows must not be summed; leaf rows are the additive breakdown.
func printFedTimingRecords(mode string, pid int, secureRun time.Duration) {
	st := mpc.SetupTiming
	cryptoSetup := st.PubKey + st.RelinKey + st.RotKey
	initOther := nonNegativeDuration(fedTimings.initTotal - cryptoSetup)
	postBlock := nonNegativeDuration(fedTimings.total - fedTimings.nullTotal - fedTimings.preBlock - fedTimings.blocks)
	protocolTotal := fedTimings.initTotal + secureRun

	printTimingRecord(mode, pid, "protocol_init", "protocol_total", "inclusive", fedTimings.initTotal, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "collective_key_setup", "protocol_init", "inclusive", cryptoSetup, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "pubkey_gen", "collective_key_setup", "leaf", st.PubKey, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "relin_key_gen", "collective_key_setup", "leaf", st.RelinKey, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "rotation_key_gen", "collective_key_setup", "leaf", st.RotKey, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "initialization_other", "protocol_init", "leaf", initOther, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "secure_run", "protocol_total", "inclusive", secureRun, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "secure_compute", "secure_run", "inclusive", fedTimings.total, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "null_model", "secure_compute", "inclusive", fedTimings.nullTotal, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "genes", "secure_compute", "inclusive", fedTimings.blocks, len(fedTimings.blockSecs), 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "post_block_finalize", "secure_compute", "inclusive", postBlock, 1, 0, 0, 0, 0, 0)
	printTimingRecord(mode, pid, "protocol_total", "none", "inclusive", protocolTotal, 1, 0, 0, 0, 0, 0)
	count, minMs, q1Ms, meanMs, q3Ms, maxMs := blockTimingDistribution(fedTimings.blockSecs)
	printTimingRecord(mode, pid, "gene_duration_distribution", "genes", "distribution", time.Duration(meanMs)*time.Millisecond,
		count, minMs, q1Ms, meanMs, q3Ms, maxMs)
}
