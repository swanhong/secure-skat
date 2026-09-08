package workflow

import (
	"context"
	"crypto/rand"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"time"

	"github.com/aead/chacha20/chacha"
)

var partyKeyNames = [][]string{
	{"shared_key_global.bin", "shared_key_0_1.bin", "shared_key_0_2.bin"},
	{"shared_key_global.bin", "shared_key_0_1.bin", "shared_key_1_2.bin"},
	{"shared_key_global.bin", "shared_key_0_2.bin", "shared_key_1_2.bin"},
}

func runKeygenCommand(args []string) error {
	flags := flag.NewFlagSet("secure-rvas keygen", flag.ContinueOnError)
	configDirectory := flags.String("config", "config/1kg", "configuration directory")
	if err := flags.Parse(args); err != nil {
		return fmt.Errorf("parse keygen arguments: %w", err)
	}
	if flags.NArg() != 0 {
		return fmt.Errorf("unexpected keygen arguments: %v", flags.Args())
	}
	return GenerateSharedPRGKeys(*configDirectory)
}

func GenerateSharedPRGKeys(configDirectory string) error {
	configs := make([]*Config, partyCount)
	for partyID := range configs {
		config, err := LoadPartyConfig(configDirectory, partyID)
		if err != nil {
			return err
		}
		if err := validatePartyConfig(config, partyID); err != nil {
			return fmt.Errorf("validate party %d config: %w", partyID, err)
		}
		configs[partyID] = config
	}

	for _, ancestry := range configs[0].Ancestries {
		keys := make(map[string][]byte, 4)
		for _, name := range []string{
			"shared_key_global.bin",
			"shared_key_0_1.bin",
			"shared_key_0_2.bin",
			"shared_key_1_2.bin",
		} {
			keys[name] = make([]byte, chacha.KeySize)
			if _, err := rand.Read(keys[name]); err != nil {
				return err
			}
		}

		for partyID, config := range configs {
			directory := filepath.Join(config.SharedKeysPath, ancestry)
			if err := os.MkdirAll(directory, 0o700); err != nil {
				return err
			}
			for _, name := range partyKeyNames[partyID] {
				if err := os.WriteFile(filepath.Join(directory, name), keys[name], 0o600); err != nil {
					return err
				}
			}
		}
	}
	return nil
}

func runPartyCommand(args []string) error {
	flags := flag.NewFlagSet("secure-rvas party", flag.ContinueOnError)
	configDirectory := flags.String("config", "config/1kg", "configuration directory")
	partyID := flags.Int("party", -1, "party ID")
	if err := flags.Parse(args); err != nil {
		return fmt.Errorf("parse party arguments: %w", err)
	}
	if flags.NArg() != 0 {
		return fmt.Errorf("unexpected party arguments: %v", flags.Args())
	}
	if *partyID < 0 || *partyID >= partyCount {
		return fmt.Errorf("party ID must be between 0 and %d", partyCount-1)
	}

	config, err := LoadPartyConfig(*configDirectory, *partyID)
	if err != nil {
		return err
	}
	if err := validatePartyConfig(config, *partyID); err != nil {
		return fmt.Errorf("validate config: %w", err)
	}
	runtime.GOMAXPROCS(config.LocalNumThreads)
	return runParty(config, *partyID)
}

func runPartyProcesses(executable, configDirectory string) ([]processMetric, error) {
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	type partyResult struct {
		partyID int
		metric  processMetric
		err     error
	}
	results := make(chan partyResult, partyCount)
	for partyID := 0; partyID < partyCount; partyID++ {
		go func(partyID int) {
			command := exec.CommandContext(
				ctx,
				executable,
				"party",
				"--config", configDirectory,
				"--party", strconv.Itoa(partyID),
			)
			command.Stdout = os.Stdout
			command.Stderr = os.Stderr

			startedAt := time.Now()
			runErr := command.Run()
			peakRSS, peakRSSErr := childProcessPeakRSSBytes(command.ProcessState)
			result := partyResult{
				partyID: partyID,
				metric: processMetric{
					process:      fmt.Sprintf("party%d", partyID),
					duration:     time.Since(startedAt),
					peakRSSBytes: peakRSS,
				},
			}
			if runErr != nil {
				result.err = fmt.Errorf("party %d failed: %w", partyID, runErr)
			} else if peakRSSErr != nil {
				result.err = fmt.Errorf("party %d peak RSS: %w", partyID, peakRSSErr)
			}
			results <- result
			if result.err != nil {
				cancel()
			}
		}(partyID)
	}

	metrics := make([]processMetric, partyCount)
	var firstError error
	for party := 0; party < partyCount; party++ {
		result := <-results
		metrics[result.partyID] = result.metric
		if result.err != nil && firstError == nil {
			firstError = result.err
		}
	}
	return metrics, firstError
}

func Run(configDirectory string) error {
	startedAt := time.Now()
	config, err := LoadPartyConfig(configDirectory, cohortAPartyID)
	if err != nil {
		return err
	}
	successPath := filepath.Join(config.RunDir, "secure", "_SUCCESS")
	if err := os.Remove(successPath); err != nil && !os.IsNotExist(err) {
		return err
	}

	executable, err := os.Executable()
	if err != nil {
		return err
	}
	partyMetrics, err := runPartyProcesses(executable, configDirectory)
	if err != nil {
		return err
	}

	parentPeakRSS, err := currentProcessPeakRSSBytes()
	if err != nil {
		return err
	}
	parentMetric := processMetric{
		process:      "parent",
		duration:     time.Since(startedAt),
		peakRSSBytes: parentPeakRSS,
	}
	processMetrics := append([]processMetric{parentMetric}, partyMetrics...)
	if err := writeProcessMetricsCSV(
		filepath.Join(config.RunDir, "metrics", "process_summary.csv"),
		processMetrics,
	); err != nil {
		return err
	}

	fmt.Printf(
		"[parent] secure_run_total: %.3fs, peak_rss: %d bytes\n",
		parentMetric.duration.Seconds(),
		parentMetric.peakRSSBytes,
	)
	cohortAMetric := partyMetrics[cohortAPartyID]
	fmt.Printf(
		"[party1] process_total: %.3fs, peak_rss: %d bytes\n",
		cohortAMetric.duration.Seconds(),
		cohortAMetric.peakRSSBytes,
	)

	if err := os.MkdirAll(filepath.Dir(successPath), 0o755); err != nil {
		return err
	}
	return os.WriteFile(successPath, []byte{}, 0o644)
}
