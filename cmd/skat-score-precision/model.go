package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"

	"github.com/BurntSushi/toml"
	"gonum.org/v1/gonum/mat"
)

const ridgeRel = 1e-6

type manifest struct {
	NGenes    int      `json:"n_genes"`
	GeneNames []string `json:"gene_symbols"`
	PublicM   []int    `json:"pub_m"`
	PrivateM  []int    `json:"priv_m"`
}

type globalConfig struct {
	CKKSParams   string `toml:"ckks_params"`
	MPCFracBits  int    `toml:"mpc_frac_bits"`
	MPCFieldSize int    `toml:"mpc_field_size"`
	NumCovs      int    `toml:"num_covs"`
	BinaryPheno  bool   `toml:"binary_pheno"`
}

type cohort struct {
	N           int
	C           int
	X           []float64 // row-major N x C, including intercept
	Y           []float64 // centered exactly as the live null path
	Residual    []float64
	SumX        []float64
	SumY        float64
	SumResidual float64
}

type contraction struct {
	Gty    []float64
	Stable []float64
	Dosage []float64
	Gtx    [][]float64 // covariate-major: C slices, each length M
}

type scoreRecord struct {
	Gene      int
	Weight    float64
	Reference float64
	Cancel64  float64
	SSLike    float64
	TermBound float64
	Gty1      float64
	Gty2      float64
	Gtx1      []float64
	Gtx2      []float64
}

func loadManifest(root string) (manifest, error) {
	var m manifest
	f, err := os.Open(filepath.Join(root, "manifest.json"))
	if err != nil {
		return m, err
	}
	defer f.Close()
	if err := json.NewDecoder(f).Decode(&m); err != nil {
		return m, err
	}
	if m.NGenes <= 0 || len(m.PublicM) != m.NGenes || len(m.PrivateM) != m.NGenes {
		return m, fmt.Errorf("manifest dimensions are inconsistent: n_genes=%d pub_m=%d priv_m=%d", m.NGenes, len(m.PublicM), len(m.PrivateM))
	}
	if len(m.GeneNames) != 0 && len(m.GeneNames) != m.NGenes {
		return m, fmt.Errorf("manifest gene_symbols has %d entries, want 0 or %d", len(m.GeneNames), m.NGenes)
	}
	return m, nil
}

func loadGlobalConfig(root string) (globalConfig, error) {
	cfg := globalConfig{CKKSParams: "PN14QP438", MPCFracBits: 30, MPCFieldSize: 256}
	path := filepath.Join(root, "config", "configGlobal.toml")
	if _, err := toml.DecodeFile(path, &cfg); err != nil {
		return cfg, fmt.Errorf("read %s: %w", path, err)
	}
	if cfg.CKKSParams == "" {
		cfg.CKKSParams = "PN14QP438"
	}
	if cfg.MPCFracBits <= 0 {
		cfg.MPCFracBits = 30
	}
	if cfg.MPCFieldSize <= 0 {
		cfg.MPCFieldSize = 256
	}
	return cfg, nil
}

func readFloatTokens(path string) ([]float64, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	s := bufio.NewScanner(f)
	s.Split(bufio.ScanWords)
	// Covariate files can be large; permit long tokens even though normal values are tiny.
	s.Buffer(make([]byte, 1024), 1024*1024)
	out := make([]float64, 0)
	for s.Scan() {
		v, err := strconv.ParseFloat(s.Text(), 64)
		if err != nil {
			return nil, fmt.Errorf("parse %s token %q: %w", path, s.Text(), err)
		}
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return nil, fmt.Errorf("%s contains a non-finite value", path)
		}
		out = append(out, v)
	}
	if err := s.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

func loadCohort(root, name string, binaryPheno bool, configuredCovs int) (cohort, error) {
	dir := filepath.Join(root, name)
	y, err := readFloatTokens(filepath.Join(dir, "pheno.txt"))
	if err != nil {
		return cohort{}, err
	}
	if len(y) == 0 {
		return cohort{}, fmt.Errorf("%s phenotype is empty", name)
	}
	cov, err := readFloatTokens(filepath.Join(dir, "cov.txt"))
	if err != nil {
		return cohort{}, err
	}
	if len(cov)%len(y) != 0 {
		return cohort{}, fmt.Errorf("%s covariate values=%d are not divisible by n=%d", name, len(cov), len(y))
	}
	n, nCov := len(y), len(cov)/len(y)
	if configuredCovs > 0 && nCov != configuredCovs {
		return cohort{}, fmt.Errorf("%s covariate columns=%d, config num_covs=%d", name, nCov, configuredCovs)
	}
	c := nCov + 1
	x := make([]float64, n*c)
	sumX := make([]float64, c)
	center := 0.0
	if binaryPheno {
		center = 0.5
	}
	for i := 0; i < n; i++ {
		x[i*c] = 1
		copy(x[i*c+1:(i+1)*c], cov[i*nCov:(i+1)*nCov])
		y[i] -= center
		sumX[0]++
		for k := 1; k < c; k++ {
			sumX[k] += x[i*c+k]
		}
	}
	sumY := 0.0
	for _, v := range y {
		sumY += v
	}
	return cohort{N: n, C: c, X: x, Y: y, SumX: sumX, SumY: sumY}, nil
}

// fitNull matches the live null model's ridge placement. Each party contributes
// ridgeRel*trace(local covariate block)/c to every non-intercept diagonal; summing
// those contributions equals the value applied below.
func fitNull(a, b *cohort, fracBits int) ([]float64, error) {
	if a.C != b.C {
		return nil, fmt.Errorf("cohort covariate widths differ: A=%d B=%d", a.C, b.C)
	}
	c := a.C
	xtx := make([]float64, c*c)
	xty := make([]float64, c)
	eps := 0.0
	for _, co := range []*cohort{a, b} {
		for i := 0; i < co.N; i++ {
			xi := co.X[i*c : (i+1)*c]
			for j := 0; j < c; j++ {
				xty[j] += xi[j] * co.Y[i]
				for k := 0; k < c; k++ {
					xtx[j*c+k] += xi[j] * xi[k]
				}
			}
		}
		trace := 0.0
		for k := 1; k < c; k++ {
			trace += xtxLocalDiag(co, k)
		}
		eps += ridgeRel * trace / float64(c)
	}
	for k := 1; k < c; k++ {
		xtx[k*c+k] += eps
	}
	var beta mat.VecDense
	if err := beta.SolveVec(mat.NewDense(c, c, xtx), mat.NewVecDense(c, xty)); err != nil {
		return nil, fmt.Errorf("solve ridge null model: %w", err)
	}
	out := make([]float64, c)
	for k := range out {
		out[k] = quantize(beta.AtVec(k), fracBits)
	}
	for _, co := range []*cohort{a, b} {
		co.Residual = make([]float64, co.N)
		co.SumResidual = 0
		for i := 0; i < co.N; i++ {
			r := co.Y[i]
			for k := 0; k < c; k++ {
				r -= co.X[i*c+k] * out[k]
			}
			co.Residual[i] = r
			co.SumResidual += r
		}
	}
	return out, nil
}

func xtxLocalDiag(co *cohort, k int) float64 {
	s := 0.0
	for i := 0; i < co.N; i++ {
		v := co.X[i*co.C+k]
		s += v * v
	}
	return s
}

func quantize(x float64, fracBits int) float64 {
	return math.Ldexp(math.Round(math.Ldexp(x, fracBits)), -fracBits)
}

// contractBlock reads one row-major int8 block and performs the same local
// minor-allele orientation as orientGenotypeLocal before computing sufficient
// statistics. Negative/missing dosages are mapped to zero like denseFromStream.
// It makes one pass over the block; flipped columns are transformed afterward
// using sum(2-g)=2*n-sum(g) and the analogous contraction identities.
func contractBlock(path string, co *cohort, m int) (contraction, error) {
	out := contraction{
		Gty:    make([]float64, m),
		Stable: make([]float64, m),
		Dosage: make([]float64, m),
		Gtx:    make([][]float64, co.C),
	}
	for k := range out.Gtx {
		out.Gtx[k] = make([]float64, m)
	}
	if m == 0 {
		return out, nil
	}
	raw, err := os.ReadFile(path)
	if err != nil {
		return out, err
	}
	want := co.N * m
	if len(raw) != want {
		return out, fmt.Errorf("%s has %d bytes, want n*m=%d", path, len(raw), want)
	}
	sumY, sumResidual, sumX := co.SumY, co.SumResidual, co.SumX
	// Unit tests and other direct callers may construct a cohort by hand.
	if len(sumX) != co.C {
		sumY, sumResidual = 0, 0
		sumX = make([]float64, co.C)
		for i := 0; i < co.N; i++ {
			sumY += co.Y[i]
			sumResidual += co.Residual[i]
			for k := 0; k < co.C; k++ {
				sumX[k] += co.X[i*co.C+k]
			}
		}
	}
	for i := 0; i < co.N; i++ {
		xi := co.X[i*co.C : (i+1)*co.C]
		for j := 0; j < m; j++ {
			iv := int8(raw[i*m+j])
			if iv <= 0 {
				continue
			}
			v := float64(iv)
			out.Dosage[j] += v
			out.Gty[j] += v * co.Y[i]
			out.Stable[j] += v * co.Residual[i]
			for k := 0; k < co.C; k++ {
				out.Gtx[k][j] += v * xi[k]
			}
		}
	}
	for j := 0; j < m; j++ {
		if out.Dosage[j] <= float64(co.N) {
			continue
		}
		out.Dosage[j] = 2*float64(co.N) - out.Dosage[j]
		out.Gty[j] = 2*sumY - out.Gty[j]
		out.Stable[j] = 2*sumResidual - out.Stable[j]
		for k := 0; k < co.C; k++ {
			out.Gtx[k][j] = 2*sumX[k] - out.Gtx[k][j]
		}
	}
	return out, nil
}

func recordsForGene(gene int, one, two contraction, beta []float64, totalN, fracBits int) ([]scoreRecord, []float64) {
	m := len(one.Gty)
	if len(two.Gty) != m {
		panic("recordsForGene: party dimensions differ")
	}
	records := make([]scoreRecord, m)
	kappas := make([]float64, m)
	for j := 0; j < m; j++ {
		r := scoreRecord{
			Gene:      gene,
			Reference: one.Stable[j] + two.Stable[j],
			Gty1:      one.Gty[j],
			Gty2:      two.Gty[j],
			Gtx1:      make([]float64, len(beta)),
			Gtx2:      make([]float64, len(beta)),
		}
		// Match scoreHE's evaluation order: each party subtracts all beta terms,
		// then the two encrypted party scores are aggregated.
		partyOne, partyTwo := r.Gty1, r.Gty2
		ssGty := quantize(r.Gty1, fracBits) + quantize(r.Gty2, fracBits)
		ssProd := 0.0
		r.TermBound = math.Abs(r.Gty1) + math.Abs(r.Gty2)
		for k, bk := range beta {
			r.Gtx1[k], r.Gtx2[k] = one.Gtx[k][j], two.Gtx[k][j]
			partyOne -= r.Gtx1[k] * bk
			partyTwo -= r.Gtx2[k] * bk
			ssGtx := quantize(r.Gtx1[k], fracBits) + quantize(r.Gtx2[k], fracBits)
			ssProd += ssGtx * bk
			r.TermBound += math.Abs(r.Gtx1[k]*bk) + math.Abs(r.Gtx2[k]*bk)
		}
		r.Cancel64 = partyOne + partyTwo
		// This mirrors scoreSS's input quantization and one post-matvec
		// truncation, but not the fixed-point Cholesky that produced betaSS.
		r.SSLike = ssGty - quantize(ssProd, fracBits)
		p := (one.Dosage[j] + two.Dosage[j]) / float64(2*totalN)
		// Each local block has already been minor-oriented, so p is in [0, 1/2]
		// up to malformed input/roundoff. Clamp only to keep the diagnostic finite.
		p = math.Max(0, math.Min(0.5, p))
		r.Weight = 25 * math.Pow(1-p, 24)
		kappas[j] = r.TermBound / math.Max(math.Abs(r.Reference), 1e-12)
		records[j] = r
	}
	return records, kappas
}

func emptyContraction(c, m int) contraction {
	out := contraction{Gty: make([]float64, m), Stable: make([]float64, m), Dosage: make([]float64, m), Gtx: make([][]float64, c)}
	for k := range out.Gtx {
		out.Gtx[k] = make([]float64, m)
	}
	return out
}

func chooseGenes(n, maxGenes int) []int {
	if maxGenes <= 0 || maxGenes >= n {
		out := make([]int, n)
		for i := range out {
			out[i] = i
		}
		return out
	}
	if maxGenes == 1 {
		return []int{0}
	}
	out := make([]int, 0, maxGenes)
	seen := make(map[int]bool, maxGenes)
	for i := 0; i < maxGenes; i++ {
		idx := int(math.Round(float64(i) * float64(n-1) / float64(maxGenes-1)))
		if !seen[idx] {
			seen[idx] = true
			out = append(out, idx)
		}
	}
	sort.Ints(out)
	return out
}

func percentile(values []float64, p float64) float64 {
	clean := make([]float64, 0, len(values))
	for _, v := range values {
		if !math.IsNaN(v) && !math.IsInf(v, 0) {
			clean = append(clean, v)
		}
	}
	if len(clean) == 0 {
		return math.NaN()
	}
	sort.Float64s(clean)
	if p <= 0 {
		return clean[0]
	}
	if p >= 1 {
		return clean[len(clean)-1]
	}
	x := p * float64(len(clean)-1)
	lo := int(math.Floor(x))
	hi := int(math.Ceil(x))
	if lo == hi {
		return clean[lo]
	}
	return clean[lo] + (x-float64(lo))*(clean[hi]-clean[lo])
}
