package workflow

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"github.com/hhcho/sfgwas/mpc"
)

type metricEvent struct {
	stage         string
	chromosome    int
	duration      time.Duration
	communication mpc.CommunicationStats
}

type metricRecorder struct {
	process     string
	printEvents bool
	events      []metricEvent
}

func newMetricRecorder(process string, printEvents bool) *metricRecorder {
	return &metricRecorder{
		process:     process,
		printEvents: printEvents,
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
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return fmt.Errorf("create metrics directory: %w", err)
	}

	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("create metrics file: %w", err)
	}

	writer := csv.NewWriter(file)
	if err := writer.Write([]string{
		"process",
		"stage",
		"chromosome",
		"duration_seconds",
		"sent_bytes",
		"received_bytes",
		"sent_message_count",
		"received_message_count",
	}); err != nil {
		file.Close()
		return fmt.Errorf("write metrics header: %w", err)
	}

	for _, event := range recorder.events {
		chromosome := ""
		if event.chromosome != 0 {
			chromosome = strconv.Itoa(event.chromosome)
		}
		communication := event.communication

		if err := writer.Write([]string{
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
		}); err != nil {
			file.Close()
			return fmt.Errorf("write metrics row: %w", err)
		}
	}

	writer.Flush()
	if err := writer.Error(); err != nil {
		file.Close()
		return fmt.Errorf("flush metrics file: %w", err)
	}
	if err := file.Close(); err != nil {
		return fmt.Errorf("close metrics file: %w", err)
	}
	return nil
}
