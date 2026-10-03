package workflow

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func Main(args []string) error {
	if len(args) == 0 {
		return fmt.Errorf(
			"usage: secure-rvas <prepare|keygen|run|party> [--config directory]",
		)
	}

	switch args[0] {
	case "prepare":
		return runPrepareCommand(args[1:])
	case "keygen":
		return runKeygenCommand(args[1:])
	case "run":
		return runRunCommand(args[1:])
	case "party":
		return runPartyCommand(args[1:])
	default:
		return fmt.Errorf(
			"unknown command %q; expected <prepare|keygen|run|party>",
			args[0],
		)
	}
}

func runPrepareCommand(args []string) error {
	flags := flag.NewFlagSet(
		"secure-rvas prepare",
		flag.ContinueOnError,
	)
	configDirectory := flags.String(
		"config",
		"config/1kg",
		"configuration directory",
	)
	clearOutputs := flags.Bool(
		"clear",
		false,
		"remove generated outputs before preprocessing",
	)

	party := flags.Int("party", 0, "prepare one site: 1=A, 2=B; 0=split A/B")
	publicVarList := flags.String("public-var-list", "", "public variant list from Party 1 (required for Party 2)")
	if err := flags.Parse(args); err != nil {
		return fmt.Errorf("parse prepare arguments: %w", err)
	}
	if flags.NArg() != 0 {
		return fmt.Errorf(
			"unexpected prepare arguments: %v",
			flags.Args(),
		)
	}

	if *party < 0 || *party > 2 || (*party == 2) != (*publicVarList != "") {
		return fmt.Errorf("prepare: --party must be 0, 1, or 2; --public-var-list is required only for Party 2")
	}

	config, err := LoadPrepareConfig(*configDirectory)
	if err != nil {
		return err
	}

	if *party == 2 {
		for _, chromosome := range config.Chromosomes {
			path := strings.ReplaceAll(*publicVarList, "{chromosome}", strconv.Itoa(chromosome))
			if _, err := os.Stat(path); err != nil {
				return fmt.Errorf("public variant list: %w", err)
			}
		}
	}
	if err := validatePrepareConfig(config); err != nil {
		return fmt.Errorf("validate config: %w", err)
	}
	if *clearOutputs {
		if err := clearGeneratedOutputs(config.RunDir); err != nil {
			return err
		}
		fmt.Printf("Cleared generated outputs from %s\n", config.RunDir)
	}

	if err := writeConfigFiles(
		*configDirectory,
		config.RunDir,
		globalConfigFilename,
		prepareConfigFilename,
	); err != nil {
		return fmt.Errorf("write config: %w", err)
	}

	fmt.Println("Running secure-rvas::prepare")
	fmt.Printf("Read configuration from %s\n", *configDirectory)
	fmt.Printf("Saved configuration to %s/config\n", config.RunDir)
	if *party == 0 {
		fmt.Printf("Shared variant rate: %g\n", config.SharedRate)
	} else {
		fmt.Printf("Preparing Party %d only\n", *party)
	}

	return prepareSite(config, *party, *publicVarList)
}

func runRunCommand(args []string) error {
	flags := flag.NewFlagSet(
		"secure-rvas run",
		flag.ContinueOnError,
	)
	configDirectory := flags.String(
		"config",
		"config/1kg",
		"configuration directory",
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
	fmt.Printf("Read configuration from %s\n", *configDirectory)

	return Run(*configDirectory)
}
