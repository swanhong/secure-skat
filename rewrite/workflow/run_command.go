package workflow

import (
	"context"
	"crypto/rand"
	"flag"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"time"

	"github.com/aead/chacha20/chacha"
)

func createSharedPRGKeys(ancestries []string) (string, error) {
	directory, err := os.MkdirTemp("", "secure-rvas-keys-")
	if err != nil {
		return "", fmt.Errorf("create shared PRG key directory: %w", err)
	}

	keyNames := []string{
		"shared_key_global.bin",
		"shared_key_0_1.bin",
		"shared_key_0_2.bin",
		"shared_key_1_2.bin",
	}

	for _, ancestry := range ancestries {
		ancestryDirectory := filepath.Join(directory, ancestry)
		if err := os.Mkdir(ancestryDirectory, 0o700); err != nil {
			os.RemoveAll(directory)
			return "", fmt.Errorf(
				"create shared PRG key directory for %s: %w",
				ancestry,
				err,
			)
		}

		for _, name := range keyNames {
			key := make([]byte, chacha.KeySize)
			if _, err := rand.Read(key); err != nil {
				os.RemoveAll(directory)
				return "", fmt.Errorf(
					"generate shared PRG key for %s: %w",
					ancestry,
					err,
				)
			}
			if err := os.WriteFile(
				filepath.Join(ancestryDirectory, name),
				key,
				0o600,
			); err != nil {
				os.RemoveAll(directory)
				return "", fmt.Errorf(
					"write shared PRG key for %s: %w",
					ancestry,
					err,
				)
			}
		}
	}
	return directory, nil
}

func runPartyCommand(args []string) error {
	flags := flag.NewFlagSet(
		"secure-rvas party",
		flag.ContinueOnError,
	)
	configPath := flags.String(
		"config",
		"run.1kg.conf",
		"path to the run configuration",
	)
	partyID := flags.Int(
		"party",
		-1,
		"Internal party ID",
	)
	sharedKeysPath := flags.String(
		"shared-keys",
		"",
		"path to the shared PRG keys",
	)

	if err := flags.Parse(args); err != nil {
		return fmt.Errorf("parse party arguments: %w", err)
	}
	if flags.NArg() != 0 {
		return fmt.Errorf(
			"unexpected party arguments: %v",
			flags.Args(),
		)
	}
	if *partyID < 0 || *partyID >= partyCount {
		return fmt.Errorf(
			"party ID must be between 0 and %d",
			partyCount-1,
		)
	}

	if *sharedKeysPath == "" {
		return fmt.Errorf("shared-keys is required")
	}

	config, err := LoadConfig(*configPath)
	if err != nil {
		return err
	}

	return runParty(config, *partyID, *sharedKeysPath)
}

func runPartyProcesses(
	executable string,
	configPath string,
	sharedKeysPath string,
) ([]processMetric, error) {
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
				"--config", configPath,
				"--party", strconv.Itoa(partyID),
				"--shared-keys", sharedKeysPath,
			)
			command.Stdout = os.Stdout
			command.Stderr = os.Stderr

			startedAt := time.Now()
			runErr := command.Run()
			peakRSS, peakRSSErr := childProcessPeakRSSBytes(
				command.ProcessState,
			)
			result := partyResult{
				partyID: partyID,
				metric: processMetric{
					process:      fmt.Sprintf("party%d", partyID),
					duration:     time.Since(startedAt),
					peakRSSBytes: peakRSS,
				},
			}
			switch {
			case runErr != nil:
				result.err = fmt.Errorf("party %d failed: %w", partyID, runErr)
			case peakRSSErr != nil:
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

func Run(configPath string) error {
	startedAt := time.Now()

	config, err := LoadConfig(configPath)
	if err != nil {
		return err
	}
	if err := ValidateConfig(config); err != nil {
		return fmt.Errorf("validate config: %w", err)
	}

	successPath := filepath.Join(
		config.RunDir,
		"secure",
		"_SUCCESS",
	)
	if err := os.Remove(successPath); err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("remove _SUCCESS file: %w", err)
	}

	executable, err := os.Executable()
	if err != nil {
		return fmt.Errorf("get executable path: %w", err)
	}

	sharedKeysPath, err := createSharedPRGKeys(config.Ancestries)
	if err != nil {
		return err
	}
	defer os.RemoveAll(sharedKeysPath)

	partyMetrics, err := runPartyProcesses(
		executable,
		configPath,
		sharedKeysPath,
	)
	if err != nil {
		return err
	}

	parentPeakRSS, err := currentProcessPeakRSSBytes()
	if err != nil {
		return fmt.Errorf("read parent peak RSS: %w", err)
	}
	parentMetric := processMetric{
		process:      "parent",
		duration:     time.Since(startedAt),
		peakRSSBytes: parentPeakRSS,
	}
	processMetrics := append(
		[]processMetric{parentMetric},
		partyMetrics...,
	)
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
		return fmt.Errorf("create secure directory: %w", err)
	}
	if err := os.WriteFile(successPath, []byte{}, 0o644); err != nil {
		return fmt.Errorf("write _SUCCESS file: %w", err)
	}

	return nil
}
