package gwas

import (
	"encoding/csv"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"testing"

	"gonum.org/v1/gonum/mat"
)

// L2 driver: run the plaintext oracle on a CSV fixture, write Q/Burden so the same
// (G, X, y) can be compared to the R::SKAT package. Skips unless LOWRANK_L2_DIR
// points at the fixture dir.
//
//	LOWRANK_L2_DIR=<dir> go test -mod=mod ./gwas/ -run TestSKATPlainL2 -count=1
//
// reads:  <dir>/G.csv (n×m), <dir>/X.csv (n×c, incl. intercept), <dir>/y.csv (n)
// writes: <dir>/go_qburden.txt  ("Q=<...>\nBurden=<...>\n")

func readCSVMatrix(t *testing.T, path string) *mat.Dense {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	rows, err := csv.NewReader(f).ReadAll()
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	nr := len(rows)
	nc := len(rows[0])
	data := make([]float64, 0, nr*nc)
	for _, r := range rows {
		for _, v := range r {
			x, err := strconv.ParseFloat(v, 64)
			if err != nil {
				t.Fatalf("parse %q in %s: %v", v, path, err)
			}
			data = append(data, x)
		}
	}
	return mat.NewDense(nr, nc, data)
}

func TestSKATPlainL2(t *testing.T) {
	dir := os.Getenv("LOWRANK_L2_DIR")
	if dir == "" {
		t.Skip("set LOWRANK_L2_DIR to the CSV fixture dir")
	}

	G := readCSVMatrix(t, filepath.Join(dir, "G.csv"))
	X := readCSVMatrix(t, filepath.Join(dir, "X.csv"))
	yMat := readCSVMatrix(t, filepath.Join(dir, "y.csv"))
	n, _ := yMat.Dims()
	y := make([]float64, n)
	for i := 0; i < n; i++ {
		y[i] = yMat.At(i, 0)
	}

	r := SKATPlain(G, X, y)

	out := fmt.Sprintf("Q=%.12e\nBurden=%.12e\nRSS=%.12e\ndof=%d\nsigma2=%.12e\n",
		r.Q, r.Burden, r.RSS, r.Dof, r.Sigma2)
	if err := os.WriteFile(filepath.Join(dir, "go_qburden.txt"), []byte(out), 0644); err != nil {
		t.Fatalf("write output: %v", err)
	}
	t.Logf("oracle on %s:\n%s", dir, out)
}
