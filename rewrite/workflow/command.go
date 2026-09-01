package workflow

import (
	"flag"
	"fmt"
)

func Main(args []string) error {
	if len(args) == 0 {
		return fmt.Errorf(
			"usage: secure-rvas <prepare|run> [--config path]",
		)
	}

	switch args[0] {
	case "prepare":
		return runPrepareCommand(args[1:])
	case "run":
		return runRunCommand(args[1:])
	case "_party":
		return runPartyCommand(args[1:])
	default:
		return fmt.Errorf(
			"unknown command %q; expected <prepare|run>",
			args[0],
		)
	}
}

func runPrepareCommand(args []string) error {
	flags := flag.NewFlagSet(
		"secure-rvas prepare",
		flag.ContinueOnError,
	)
	configPath := flags.String(
		"config",
		"run.1kg.conf",
		"path to the run configuration",
	)
	clearOutputs := flags.Bool(
		"clear",
		false,
		"remove generated outputs before preprocessing",
	)

	if err := flags.Parse(args); err != nil {
		return fmt.Errorf("parse prepare arguments: %w", err)
	}
	if flags.NArg() != 0 {
		return fmt.Errorf(
			"unexpected prepare arguments: %v",
			flags.Args(),
		)
	}

	config, err := LoadConfig(*configPath)
	if err != nil {
		return err
	}

	if err := ValidateConfig(config); err != nil {
		return fmt.Errorf("validate config: %w", err)
	}
	if *clearOutputs {
		if err := clearGeneratedOutputs(config.RunDir); err != nil {
			return err
		}
		fmt.Printf("Cleared generated outputs from %s\n", config.RunDir)
	}

	if err := writeConfig(*configPath, config.RunDir); err != nil {
		return fmt.Errorf("write config: %w", err)
	}

	fmt.Println("Running secure-rvas::prepare")
	fmt.Printf("Read configuration from %s\n", *configPath)
	fmt.Printf("Saved run configuration to %s/%s\n", config.RunDir, runConfigFilename)

	return Prepare(config)
}

func runRunCommand(args []string) error {
	flags := flag.NewFlagSet(
		"secure-rvas run",
		flag.ContinueOnError,
	)
	configPath := flags.String(
		"config",
		"run.1kg.conf",
		"path to the run configuration",
	)

	if err := flags.Parse(args); err != nil {
		return fmt.Errorf("parse run arguments: %w", err)
	}
	if flags.NArg() != 0 {
		return fmt.Errorf(
			"unexpected run arguments: %v",
			flags.Args(),
		)
	}

	fmt.Println("Running secure-rvas::run")
	fmt.Printf("Read configuration from %s\n", *configPath)

	return Run(*configPath)
}
