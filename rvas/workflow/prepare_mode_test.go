package workflow

import (
	"os"
	"path/filepath"
	"testing"
)

func TestClearPrepareMode(t *testing.T) {
	joint := []string{"prepared/EUR/chr21/data", "cache/prepared/EUR/chr21/data", "secure/result"}
	preserved := []string{"public_variants.tsv", "local_cut_variants.tsv"}
	local := []string{"prepared/local/EUR/chr21/data", "cache/prepared/local/EUR/chr21/data", "local/result"}
	for _, mode := range []string{"both", "joint", "local"} {
		t.Run(mode, func(t *testing.T) {
			root := t.TempDir()
			for _, name := range append(append(append([]string{}, joint...), local...), preserved...) {
				path := filepath.Join(root, name)
				if err := os.MkdirAll(filepath.Dir(path), 0700); err != nil {
					t.Fatal(err)
				}
				if err := os.WriteFile(path, []byte("input"), 0600); err != nil {
					t.Fatal(err)
				}
			}
			if err := clearGeneratedOutputs(root, mode); err != nil {
				t.Fatal(err)
			}
			for _, name := range preserved {
				data, err := os.ReadFile(filepath.Join(root, name))
				if err != nil || string(data) != "input" {
					t.Fatalf("mode %s changed input %s: %v", mode, name, err)
				}
			}
			for kind, names := range map[string][]string{"joint": joint, "local": local} {
				for _, name := range names {
					_, err := os.Stat(filepath.Join(root, name))
					removed := os.IsNotExist(err)
					if removed != (mode == "both" || mode == kind) {
						t.Fatalf("mode %s: unexpected removal of %s: %v", mode, name, err)
					}
				}
			}
		})
	}
}
