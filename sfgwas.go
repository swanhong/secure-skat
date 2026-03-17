package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"

	"github.com/BurntSushi/toml"
	"github.com/hhcho/sfgwas/gwas"
	"github.com/raulk/go-watchdog"
)

// Expects a party ID provided as an environment variable;
// e.g., run "PID=1 go run sfgwas.go"
var PID, PID_ERR = strconv.Atoi(os.Getenv("PID"))

// Default config path
var CONFIG_PATH = "config/"

func main() {
	RunProtocol()
}

func InitProtocol(configPath string) *gwas.ProtocolInfo {
	config := new(gwas.Config)

	// Import global parameters
	if _, err := toml.DecodeFile(filepath.Join(configPath, "configGlobal.toml"), config); err != nil {
		fmt.Println(err)
		return nil
	}

	// Import local parameters
	if _, err := toml.DecodeFile(filepath.Join(configPath, fmt.Sprintf("configLocal.Party%d.toml", PID)), config); err != nil {
		fmt.Println(err)
		return nil
	}

	if suffix := strings.TrimSpace(os.Getenv("SFGWAS_OUTPUT_SUFFIX")); suffix != "" {
		config.OutDir = config.OutDir + "_" + suffix
		config.CacheDir = config.CacheDir + "_" + suffix
	}

	// Create cache/output directories
	if err := os.MkdirAll(config.CacheDir, 0755); err != nil {
		panic(err)
	}
	if err := os.MkdirAll(config.OutDir, 0755); err != nil {
		panic(err)
	}

	// Set max number of threads
	runtime.GOMAXPROCS(config.LocalNumThreads)

	return gwas.InitializeGWASProtocol(config, PID, false)
}

func parseRareVariantMode() gwas.RareVariantMode {
	mode := strings.ToLower(strings.TrimSpace(os.Getenv("SFGWAS_MODE")))
	switch mode {
	case "", string(gwas.RareVariantModeGWAS):
		return gwas.RareVariantModeGWAS
	case string(gwas.RareVariantModeSKAT):
		return gwas.RareVariantModeSKAT
	case string(gwas.RareVariantModeBurden):
		return gwas.RareVariantModeBurden
	case string(gwas.RareVariantModeSKATO):
		return gwas.RareVariantModeSKATO
	default:
		panic(fmt.Sprintf("unsupported SFGWAS_MODE: %s", mode))
	}
}

func parseSKATORho() float64 {
	rhoStr := strings.TrimSpace(os.Getenv("SKATO_RHO"))
	if rhoStr == "" {
		return 0.5
	}

	rho, err := strconv.ParseFloat(rhoStr, 64)
	if err != nil {
		panic(fmt.Sprintf("invalid SKATO_RHO: %v", err))
	}
	if rho < 0.0 || rho > 1.0 {
		panic("SKATO_RHO must be in [0, 1]")
	}
	return rho
}

func RunProtocol() {
	if PID_ERR != nil {
		panic(PID_ERR)
	}

	// Initialize protocol
	prot := InitProtocol(CONFIG_PATH)

	// Invoke memory manager
	err, stopFn := watchdog.HeapDriven(prot.GetConfig().MemoryLimit, 40, watchdog.NewAdaptivePolicy(0.5))
	if err != nil {
		panic(err)
	}
	defer stopFn()

	mode := parseRareVariantMode()
	switch mode {
	case gwas.RareVariantModeGWAS:
		prot.GWAS()
	case gwas.RareVariantModeSKAT:
		prot.SKAT()
	case gwas.RareVariantModeBurden:
		prot.Burden()
	case gwas.RareVariantModeSKATO:
		prot.SKATO(parseSKATORho())
	default:
		panic(fmt.Sprintf("unsupported protocol mode: %s", mode))
	}

	prot.SyncAndTerminate(true)
}
