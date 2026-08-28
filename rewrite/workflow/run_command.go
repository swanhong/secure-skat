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

	"github.com/aead/chacha20/chacha"
)

func createSharedPRGKeys() (string, error) {
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

	for _, name := range keyNames {
		key := make([]byte, chacha.KeySize)
		if _, err := rand.Read(key); err != nil {
			os.RemoveAll(directory)
			return "", fmt.Errorf("generate shared PRG key %s: %w", name, err)
		}
		if err := os.WriteFile(
			filepath.Join(directory, name),
			key,
			0o600,
		); err != nil {
			os.RemoveAll(directory)
			return "", fmt.Errorf("write shared PRG key %s: %w", name, err)
		}
	}
	return directory, nil
}

func runPartyCommand(args []string) error {
	flags := flag.NewFlagSet(
		"secure-rvas _party",
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
) error {
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()

	results := make(chan error, partyCount)

	for partyID := 0; partyID < partyCount; partyID++ {
		go func(partyID int) {
			command := exec.CommandContext(
				ctx,
				executable,
				"_party",
				"--config", configPath,
				"--party", strconv.Itoa(partyID),
				"--shared-keys", sharedKeysPath,
			)
			command.Stdout = os.Stdout
			command.Stderr = os.Stderr

			if err := command.Run(); err != nil {
				results <- fmt.Errorf(
					"party %d failed: %w",
					partyID,
					err,
				)
				cancel()
				return
			}
			results <- nil
		}(partyID)
	}

	var firstError error
	for party := 0; party < partyCount; party++ {
		if err := <-results; err != nil && firstError == nil {
			firstError = err
		}
	}
	return firstError
}

func Run(configPath string) error {
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

	sharedKeysPath, err := createSharedPRGKeys()
	if err != nil {
		return err
	}
	defer os.RemoveAll(sharedKeysPath)

	if err := runPartyProcesses(
		executable,
		configPath,
		sharedKeysPath,
	); err != nil {
		return err
	}

	if err := os.MkdirAll(filepath.Dir(successPath), 0o755); err != nil {
		return fmt.Errorf("create secure directory: %w", err)
	}
	if err := os.WriteFile(successPath, []byte{}, 0o644); err != nil {
		return fmt.Errorf("write _SUCCESS file: %w", err)
	}

	return nil
}
