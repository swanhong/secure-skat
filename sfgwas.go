package main

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"

	"github.com/BurntSushi/toml"
	"github.com/hhcho/sfgwas/gwas"
	"github.com/raulk/go-watchdog"
)

// Expects a party ID provided as an environment variable;
// e.g., run "PID=1 go run sfgwas.go"
var PID, PID_ERR = strconv.Atoi(os.Getenv("PID"))

// CONFIG_PATH defaults to config/, overridable via SFGWAS_CONFIG_PATH.
var CONFIG_PATH = envOr("SFGWAS_CONFIG_PATH", "config/")

func envOr(key, def string) string {
	if v := os.Getenv(key); v != "" {
		return v
	}
	return def
}

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

	// Dispatch by SFGWAS_MODE (default gwas): single-variant GWAS or a rare-variant test.
	switch mode := envOr("SFGWAS_MODE", "gwas"); mode {
	case "gwas":
		prot.GWAS()
	case "skat":
		prot.SKAT()
	case "burden":
		prot.Burden()
	case "skato":
		rho := 0.5
		if v, err := strconv.ParseFloat(os.Getenv("SKATO_RHO"), 64); err == nil {
			rho = v
		}
		prot.SKATO(rho)
	case "skat_fed":
		prot.SKATFederatedPrivate()
	default:
		panic(fmt.Sprintf("unsupported SFGWAS_MODE: %s (gwas|skat|burden|skato|skat_fed)", mode))
	}

	prot.SyncAndTerminate(true)
}
