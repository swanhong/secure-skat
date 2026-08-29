package workflow

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/hhcho/sfgwas/mpc"
)

type stageDefinition struct {
	parent          string
	measurementKind string
}

var stageDefinitions = map[string]stageDefinition{
	"input_loading":            {measurementKind: "leaf"},
	"network_init":             {parent: "session_setup", measurementKind: "leaf"},
	"sample_count_exchange":    {parent: "session_setup", measurementKind: "leaf"},
	"collective_setup":         {parent: "session_setup", measurementKind: "inclusive"},
	"pubkey_gen":               {parent: "collective_setup", measurementKind: "leaf"},
	"relin_key_gen":            {parent: "collective_setup", measurementKind: "leaf"},
	"rotkey_gen":               {parent: "collective_setup", measurementKind: "leaf"},
	"null_model":               {parent: "session_setup", measurementKind: "inclusive"},
	"null_local_equations":     {parent: "null_model", measurementKind: "leaf"},
	"null_aggregate_shares":    {parent: "null_model", measurementKind: "leaf"},
	"null_factor_solve":        {parent: "null_model", measurementKind: "leaf"},
	"null_rss":                 {parent: "null_model", measurementKind: "leaf"},
	"chromosome_total":         {measurementKind: "inclusive"},
	"compute_weights":          {parent: "chromosome_total", measurementKind: "leaf"},
	"packed_statistics":        {parent: "chromosome_total", measurementKind: "inclusive"},
	"beta_packing":             {parent: "packed_statistics", measurementKind: "leaf"},
	"batch_preparation":        {parent: "packed_statistics", measurementKind: "leaf"},
	"public_ql":                {parent: "packed_statistics", measurementKind: "leaf"},
	"private_ql":               {parent: "packed_statistics", measurementKind: "leaf"},
	"kernel_inputs":            {parent: "packed_statistics", measurementKind: "leaf"},
	"first_gtg_action":         {parent: "packed_statistics", measurementKind: "leaf"},
	"second_gtg_action":        {parent: "packed_statistics", measurementKind: "leaf"},
	"private_trace_correction": {parent: "packed_statistics", measurementKind: "leaf"},
	"burden_variance":          {parent: "packed_statistics", measurementKind: "leaf"},
	"finalize":                 {parent: "chromosome_total", measurementKind: "inclusive"},
	"alpha_score_assembly":     {parent: "finalize", measurementKind: "leaf"},
	"burden_statistic":         {parent: "finalize", measurementKind: "leaf"},
	"skat_wilson_hilferty":     {parent: "finalize", measurementKind: "leaf"},
	"release":                  {parent: "chromosome_total", measurementKind: "inclusive"},
	"shares_to_ciphertext":     {parent: "release", measurementKind: "leaf"},
	"collective_decrypt":       {parent: "release", measurementKind: "leaf"},
	"local_pvalues":            {parent: "release", measurementKind: "leaf"},
	"write_results":            {measurementKind: "leaf"},
}

type metricEvent struct {
	stage         string
	chromosome    int
	duration      time.Duration
	communication mpc.CommunicationStats
	count         int
}

type processMetric struct {
	process      string
	duration     time.Duration
	peakRSSBytes uint64
}

type metricRecorder struct {
	process     string
	ancestry    string
	printEvents bool
	startedAt   time.Time
	events      []metricEvent
}

func newMetricRecorder(
	process string,
	ancestry string,
	printEvents bool,
) *metricRecorder {
	return &metricRecorder{
		process:     process,
		ancestry:    ancestry,
		printEvents: printEvents,
		startedAt:   time.Now(),
	}
}

func (recorder *metricRecorder) start(
	stage string,
	chromosome int,
	networks mpc.ParallelNetworks,
) func() {
	return recorder.startEvent(stage, chromosome, networks, true)
}

func (recorder *metricRecorder) observe(
	chromosome int,
	networks mpc.ParallelNetworks,
) func(stage string) func() {
	return func(stage string) func() {
		return recorder.startEvent(stage, chromosome, networks, false)
	}
}

func (recorder *metricRecorder) startEvent(
	stage string,
	chromosome int,
	networks mpc.ParallelNetworks,
	printEvent bool,
) func() {
	startedAt := time.Now()
	var startCommunication mpc.CommunicationStats
	if len(networks) > 0 {
		startCommunication = networks.GetCommunicationStats()
	}

	return func() {
		communication := mpc.CommunicationStats{}
		if len(networks) > 0 {
			communication = networks.GetCommunicationStats().Sub(
				startCommunication,
			)
		}

		event := metricEvent{
			stage:         stage,
			chromosome:    chromosome,
			duration:      time.Since(startedAt),
			communication: communication,
			count:         1,
		}
		recorder.addEvent(event)

		if !recorder.printEvents || !printEvent {
			return
		}
		if event.chromosome == 0 {
			fmt.Printf(
				"[%s %s] %s: %.3fs\n",
				recorder.process,
				recorder.ancestry,
				event.stage,
				event.duration.Seconds(),
			)
		} else {
			fmt.Printf(
				"[%s %s] chr%d %s: %.3fs\n",
				recorder.process,
				recorder.ancestry,
				event.chromosome,
				event.stage,
				event.duration.Seconds(),
			)
		}
	}
}

func (recorder *metricRecorder) addDuration(
	stage string,
	chromosome int,
	duration time.Duration,
) {
	recorder.addEvent(metricEvent{
		stage:      stage,
		chromosome: chromosome,
		duration:   duration,
		count:      1,
	})
}

func (recorder *metricRecorder) addEvent(event metricEvent) {
	for index := range recorder.events {
		current := &recorder.events[index]
		if current.stage != event.stage || current.chromosome != event.chromosome {
			continue
		}

		current.duration += event.duration
		current.communication.SentBytes += event.communication.SentBytes
		current.communication.ReceivedBytes += event.communication.ReceivedBytes
		current.communication.SentMessages += event.communication.SentMessages
		current.communication.ReceivedMessages += event.communication.ReceivedMessages
		current.count += event.count
		return
	}
	recorder.events = append(recorder.events, event)
}

func (recorder *metricRecorder) writeCSV(path string) error {
	rows := make([][]string, 0, len(recorder.events))
	for _, event := range recorder.events {
		definition := stageDefinitions[event.stage]
		chromosome := ""
		if event.chromosome != 0 {
			chromosome = strconv.Itoa(event.chromosome)
		}
		communication := event.communication
		rows = append(rows, []string{
			recorder.process,
			recorder.ancestry,
			event.stage,
			definition.parent,
			definition.measurementKind,
			chromosome,
			strconv.Itoa(event.count),
			strconv.FormatFloat(
				event.duration.Seconds(),
				'f',
				9,
				64,
			),
			strconv.FormatUint(communication.SentBytes, 10),
			strconv.FormatUint(communication.ReceivedBytes, 10),
			strconv.FormatUint(communication.SentMessages, 10),
			strconv.FormatUint(communication.ReceivedMessages, 10),
		})
	}

	return writeCSV(path, []string{
		"process",
		"ancestry",
		"stage",
		"parent_stage",
		"measurement_kind",
		"chromosome",
		"count",
		"duration_seconds",
		"sent_bytes",
		"received_bytes",
		"sent_message_count",
		"received_message_count",
	}, rows)
}

func (recorder *metricRecorder) printTimeTree() {
	if !recorder.printEvents {
		return
	}
	fmt.Print(recorder.timeTree(time.Since(recorder.startedAt)))
}

func (recorder *metricRecorder) timeTree(total time.Duration) string {
	stageDuration := func(stage string, chromosome int) time.Duration {
		for _, event := range recorder.events {
			if event.stage == stage && event.chromosome == chromosome {
				return event.duration
			}
		}
		return 0
	}
	format := func(duration time.Duration) time.Duration {
		return duration.Round(time.Millisecond)
	}

	setup := stageDuration("network_init", 0) +
		stageDuration("sample_count_exchange", 0) +
		stageDuration("collective_setup", 0) +
		stageDuration("null_model", 0)
	inputLoading := stageDuration("input_loading", 0)
	writeResults := stageDuration("write_results", 0)
	classified := inputLoading + setup + writeResults

	var tree strings.Builder
	fmt.Fprintf(
		&tree,
		"\n[%s %s] timing tree:\n",
		recorder.process,
		recorder.ancestry,
	)
	fmt.Fprintf(&tree, "  ├─ input loading          %v\n", format(inputLoading))
	fmt.Fprintf(&tree, "  ├─ session setup          %v\n", format(setup))
	fmt.Fprintf(&tree, "  │  ├─ network init         %v\n", format(stageDuration("network_init", 0)))
	fmt.Fprintf(&tree, "  │  ├─ sample count exchange %v\n", format(stageDuration("sample_count_exchange", 0)))
	fmt.Fprintf(&tree, "  │  ├─ collective setup     %v\n", format(stageDuration("collective_setup", 0)))
	fmt.Fprintf(&tree, "  │  │  ├─ PubKeyGen          %v\n", format(stageDuration("pubkey_gen", 0)))
	fmt.Fprintf(&tree, "  │  │  ├─ RelinKeyGen        %v\n", format(stageDuration("relin_key_gen", 0)))
	fmt.Fprintf(&tree, "  │  │  └─ RotKeyGen          %v\n", format(stageDuration("rotkey_gen", 0)))
	fmt.Fprintf(&tree, "  │  └─ null model           %v\n", format(stageDuration("null_model", 0)))
	fmt.Fprintf(&tree, "  │     ├─ local XTX/XTY/YTY  %v\n", format(stageDuration("null_local_equations", 0)))
	fmt.Fprintf(&tree, "  │     ├─ aggregate shares   %v\n", format(stageDuration("null_aggregate_shares", 0)))
	fmt.Fprintf(&tree, "  │     ├─ factor and solve   %v\n", format(stageDuration("null_factor_solve", 0)))
	fmt.Fprintf(&tree, "  │     └─ RSS                %v\n", format(stageDuration("null_rss", 0)))

	for _, event := range recorder.events {
		if event.stage != "chromosome_total" {
			continue
		}
		chromosome := event.chromosome
		classified += event.duration
		fmt.Fprintf(&tree, "  ├─ chr%-19d %v\n", chromosome, format(event.duration))
		fmt.Fprintf(&tree, "  │  ├─ weights              %v\n", format(stageDuration("compute_weights", chromosome)))
		fmt.Fprintf(&tree, "  │  ├─ packed statistics    %v\n", format(stageDuration("packed_statistics", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ beta packing       %v\n", format(stageDuration("beta_packing", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ batch preparation  %v\n", format(stageDuration("batch_preparation", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ public Q/L         %v\n", format(stageDuration("public_ql", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ private Q/L        %v\n", format(stageDuration("private_ql", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ kernel inputs      %v\n", format(stageDuration("kernel_inputs", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ first GtG action   %v\n", format(stageDuration("first_gtg_action", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ second GtG action  %v\n", format(stageDuration("second_gtg_action", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ private trace correction %v\n", format(stageDuration("private_trace_correction", chromosome)))
		fmt.Fprintf(&tree, "  │  │  └─ burden variance    %v\n", format(stageDuration("burden_variance", chromosome)))
		fmt.Fprintf(&tree, "  │  ├─ finalize             %v\n", format(stageDuration("finalize", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ alpha and score assembly %v\n", format(stageDuration("alpha_score_assembly", chromosome)))
		fmt.Fprintf(&tree, "  │  │  ├─ burden statistic   %v\n", format(stageDuration("burden_statistic", chromosome)))
		fmt.Fprintf(&tree, "  │  │  └─ SKAT Wilson-Hilferty %v\n", format(stageDuration("skat_wilson_hilferty", chromosome)))
		fmt.Fprintf(&tree, "  │  └─ release              %v\n", format(stageDuration("release", chromosome)))
		fmt.Fprintf(&tree, "  │     ├─ shares to ciphertext %v\n", format(stageDuration("shares_to_ciphertext", chromosome)))
		fmt.Fprintf(&tree, "  │     ├─ collective decrypt %v\n", format(stageDuration("collective_decrypt", chromosome)))
		fmt.Fprintf(&tree, "  │     └─ local p-values     %v\n", format(stageDuration("local_pvalues", chromosome)))
	}

	other := total - classified
	if other < 0 {
		other = 0
	}
	fmt.Fprintf(&tree, "  ├─ write results          %v\n", format(writeResults))
	fmt.Fprintf(&tree, "  ├─ other overhead         %v\n", format(other))
	fmt.Fprintf(&tree, "  └─ TOTAL                  %v\n", format(total))
	return tree.String()
}

func writeProcessMetricsCSV(path string, metrics []processMetric) error {
	rows := make([][]string, 0, len(metrics))
	for _, metric := range metrics {
		rows = append(rows, []string{
			metric.process,
			strconv.FormatFloat(metric.duration.Seconds(), 'f', 9, 64),
			strconv.FormatUint(metric.peakRSSBytes, 10),
		})
	}

	return writeCSV(path, []string{
		"process",
		"duration_seconds",
		"peak_rss_bytes",
	}, rows)
}

func writeCSV(path string, header []string, rows [][]string) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return fmt.Errorf("create CSV directory: %w", err)
	}

	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create CSV file: %w", err)
	}

	writer := csv.NewWriter(file)
	writer.WriteAll(append([][]string{header}, rows...))
	if err := writer.Error(); err != nil {
		file.Close()
		return fmt.Errorf("write CSV file: %w", err)
	}
	if err := file.Close(); err != nil {
		return fmt.Errorf("close CSV file: %w", err)
	}
	return nil
}

func peakRSSBytes(usage *syscall.Rusage) (uint64, error) {
	peakRSS := uint64(usage.Maxrss)
	switch runtime.GOOS {
	case "darwin":
		return peakRSS, nil
	case "linux":
		return peakRSS * 1024, nil
	default:
		return 0, fmt.Errorf("peak RSS is not supported on %s", runtime.GOOS)
	}
}

func currentProcessPeakRSSBytes() (uint64, error) {
	var usage syscall.Rusage
	if err := syscall.Getrusage(syscall.RUSAGE_SELF, &usage); err != nil {
		return 0, fmt.Errorf("read process resource usage: %w", err)
	}
	return peakRSSBytes(&usage)
}

func childProcessPeakRSSBytes(state *os.ProcessState) (uint64, error) {
	if state == nil {
		return 0, fmt.Errorf("child process state is unavailable")
	}
	usage, ok := state.SysUsage().(*syscall.Rusage)
	if !ok {
		return 0, fmt.Errorf("child process resource usage is unavailable")
	}
	return peakRSSBytes(usage)
}
