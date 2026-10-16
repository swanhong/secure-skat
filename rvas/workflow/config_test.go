package workflow

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestPrepareRunDir(t *testing.T) {
	directory := t.TempDir()
	runDir := filepath.Join(directory, "run")
	files := map[string]string{
		"configGlobal.toml":  "",
		"configLocal.toml":   "shared_keys_path = \"{run_dir}/secure/keys\"\n",
		"configPrepare.toml": "run_dir = \"" + runDir + "\"\n",
	}
	write := func(name, contents string) {
		t.Helper()
		if err := os.WriteFile(filepath.Join(directory, name), []byte(contents), 0o600); err != nil {
			t.Fatal(err)
		}
	}
	for name, contents := range files {
		write(name, contents)
	}
	prepared, err := LoadPrepareConfig(directory)
	if err != nil {
		t.Fatal(err)
	}
	party, err := LoadPartyConfig(directory, 1)
	if err != nil {
		t.Fatal(err)
	}
	if prepared.RunDir != runDir || party.RunDir != runDir || party.SharedKeysPath != runDir+"/secure/keys" {
		t.Fatal("prepare and party must use run_dir from Prepare")
	}
	if prepareRequestFromConfig(prepared).PreparedCacheDir != filepath.Join(runDir, "cache", "prepared") {
		t.Fatal("cache must be inside run_dir")
	}
	for _, name := range []string{"configGlobal.toml", "configLocal.toml"} {
		write(name, files[name]+"run_dir = \"old-path\"\n")
		if _, err := LoadPrepareConfig(directory); err == nil || !strings.Contains(err.Error(), "run_dir belongs in") {
			t.Fatalf("expected migration error for %s, got %v", name, err)
		}
		write(name, files[name])
	}
}
