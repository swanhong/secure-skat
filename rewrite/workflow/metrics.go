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

type metricEvent struct {
	stage         string
	chromosome    int
	duration      time.Duration
	communication mpc.CommunicationStats
}

type processMetric struct {
	process      string
	duration     time.Duration
	peakRSSBytes uint64
}

type metricRecorder struct {
	process     string
	printEvents bool
	startedAt   time.Time
	events      []metricEvent
}

func newMetricRecorder(process string, printEvents bool) *metricRecorder {
	return &metricRecorder{
		process:     process,
		printEvents: printEvents,
		startedAt:   time.Now(),
	}
}

func (recorder *metricRecorder) start(
	stage string,
	chromosome int,
	networks mpc.ParallelNetworks,
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
		}
		recorder.events = append(recorder.events, event)

		if !recorder.printEvents {
			return
		}
		if event.chromosome == 0 {
			fmt.Printf(
				"[%s] %s: %.3fs\n",
				recorder.process,
				event.stage,
				event.duration.Seconds(),
			)
		} else {
			fmt.Printf(
				"[%s] chr%d %s: %.3fs\n",
				recorder.process,
				event.chromosome,
				event.stage,
				event.duration.Seconds(),
			)
		}
	}
}

func (recorder *metricRecorder) writeCSV(path string) error {
	rows := make([][]string, 0, len(recorder.events))
	for _, event := range recorder.events {
		chromosome := ""
		if event.chromosome != 0 {
			chromosome = strconv.Itoa(event.chromosome)
		}
		communication := event.communication
		rows = append(rows, []string{
			recorder.process,
			event.stage,
			chromosome,
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
		"stage",
		"chromosome",
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
		stageDuration("setup_null", 0)
	writeResults := stageDuration("write_results", 0)
	classified := setup + writeResults

	var tree strings.Builder
	fmt.Fprintf(&tree, "\n[%s] timing tree:\n", recorder.process)
	fmt.Fprintf(&tree, "  ├─ session setup          %v\n", format(setup))
	fmt.Fprintf(&tree, "  │  ├─ network init         %v\n", format(stageDuration("network_init", 0)))
	fmt.Fprintf(&tree, "  │  ├─ sample count exchange %v\n", format(stageDuration("sample_count_exchange", 0)))
	fmt.Fprintf(&tree, "  │  ├─ collective setup     %v\n", format(stageDuration("collective_setup", 0)))
	fmt.Fprintf(&tree, "  │  └─ setup null           %v\n", format(stageDuration("setup_null", 0)))

	for _, event := range recorder.events {
		if event.stage != "chromosome_total" {
			continue
		}
		chromosome := event.chromosome
		classified += event.duration
		fmt.Fprintf(&tree, "  ├─ chr%-19d %v\n", chromosome, format(event.duration))
		fmt.Fprintf(&tree, "  │  ├─ compute weights      %v\n", format(stageDuration("compute_weights", chromosome)))
		fmt.Fprintf(&tree, "  │  ├─ packed statistics    %v\n", format(stageDuration("packed_statistics", chromosome)))
		fmt.Fprintf(&tree, "  │  ├─ finalize             %v\n", format(stageDuration("finalize", chromosome)))
		fmt.Fprintf(&tree, "  │  └─ release              %v\n", format(stageDuration("release", chromosome)))
	}

	fmt.Fprintf(&tree, "  ├─ write results          %v\n", format(writeResults))
	fmt.Fprintf(&tree, "  ├─ other overhead         %v\n", format(total-classified))
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
