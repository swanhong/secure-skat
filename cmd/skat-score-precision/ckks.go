package main

import (
	"fmt"
	"runtime"

	securecrypto "github.com/hhcho/sfgwas/crypto"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

// These custom N14 parameter sets are diagnostic-only. They keep the same
// ring dimension and 8192-slot capacity as PN14QP438 while matching each
// rescaling prime to a higher default scale. Production use still requires the
// collective SS<->CKKS, weight, moment, and end-to-end p-value gates described
// in README.md.
var (
	ckksComplexParamsN14QP431S40 = ckks.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{51, 40, 40, 40, 40, 40, 40, 40},
		LogP:            []int{50, 50},
		LogDefaultScale: 40,
	}
)

type heEngine struct {
	CP              *securecrypto.CryptoParams
	Beta            securecrypto.CipherVector // one full-slot replicated ciphertext per coefficient
	C               int
	Slots           int
	Params          string
	LogN            int
	LogQP           float64
	LogQ            []int
	LogP            []int
	LogDefaultScale int
	MaxLevel        int
	DiagnosticOnly  bool
}

func diagnosticOnlyParameter(name string) bool {
	return name == "PN14QP431S40" || name == "PN14QP436S45"
}

func parameterLiteral(name string) (ckks.ParametersLiteral, error) {
	if name == "PN14QP431S40" {
		return ckksComplexParamsN14QP431S40, nil
	}
	return securecrypto.ResolveCKKSParametersLiteral(name)
}

func newHEEngine(name string, beta []float64, precision uint, threads int) (*heEngine, error) {
	lit, err := parameterLiteral(name)
	if err != nil {
		return nil, err
	}
	params, err := ckks.NewParametersFromLiteral(lit)
	if err != nil {
		return nil, err
	}
	if threads <= 0 {
		threads = runtime.GOMAXPROCS(0)
	}
	kgen := ckks.NewKeyGenerator(params)
	sk := kgen.GenSecretKeyNew()
	pk := kgen.GenPublicKeyNew(sk)
	// This diagnostic uses only plaintext-ciphertext multiplication, which keeps
	// ciphertext degree one and needs neither relinearization nor rotation keys.
	cp := securecrypto.NewCryptoParams(params, sk, sk.CopyNew(), pk, nil, precision, threads)
	slots := cp.GetSlots()
	betaCT := make(securecrypto.CipherVector, len(beta))
	repeated := make([]float64, slots)
	for k, b := range beta {
		for i := range repeated {
			repeated[i] = b
		}
		ct, _ := securecrypto.EncryptFloatVector(cp, repeated)
		betaCT[k] = ct[0]
	}
	return &heEngine{
		CP:              cp,
		Beta:            betaCT,
		C:               len(beta),
		Slots:           slots,
		Params:          name,
		LogN:            params.LogN(),
		LogQP:           params.LogQP(),
		LogQ:            params.LogQi(),
		LogP:            params.LogPi(),
		LogDefaultScale: params.LogDefaultScale(),
		MaxLevel:        params.MaxLevel(),
		DiagnosticOnly:  diagnosticOnlyParameter(name),
	}, nil
}

func cipherLevel(v securecrypto.CipherVector) int {
	if len(v) == 0 || v[0] == nil {
		return 0
	}
	level := v[0].Level()
	for _, ct := range v[1:] {
		if ct != nil && ct.Level() < level {
			level = ct.Level()
		}
	}
	return level
}

func alignLevels(cp *securecrypto.CryptoParams, a, b securecrypto.CipherVector) (securecrypto.CipherVector, securecrypto.CipherVector) {
	target := cipherLevel(a)
	if l := cipherLevel(b); l < target {
		target = l
	}
	if cipherLevel(a) != target {
		a = securecrypto.DropLevel(cp, securecrypto.CipherMatrix{a}, target)[0]
	}
	if cipherLevel(b) != target {
		b = securecrypto.DropLevel(cp, securecrypto.CipherMatrix{b}, target)[0]
	}
	return a, b
}

// partyScore mirrors gwas.scoreHE exactly: fresh Enc(G^T y), then one
// plaintext-column times the common replicated beta ciphertext per covariate.
func (e *heEngine) partyScore(gty []float64, gtx [][]float64) securecrypto.CipherVector {
	if len(gty) == 0 {
		return nil
	}
	score, _ := securecrypto.EncryptFloatVector(e.CP, gty)
	for k := 0; k < e.C; k++ {
		pt, _ := securecrypto.EncodeFloatVector(e.CP, gtx[k])
		term := securecrypto.CPMult(e.CP, securecrypto.CipherVector{e.Beta[k]}, pt)
		score, term = alignLevels(e.CP, score, term)
		score = securecrypto.CSub(e.CP, score, term)
	}
	return score
}

func (e *heEngine) evaluate(gty1 []float64, gtx1 [][]float64, gty2 []float64, gtx2 [][]float64, twoParties bool) []float64 {
	one := e.partyScore(gty1, gtx1)
	out := one
	if twoParties {
		two := e.partyScore(gty2, gtx2)
		one, two = alignLevels(e.CP, one, two)
		out = securecrypto.CAdd(e.CP, one, two)
	}
	return securecrypto.DecryptFloatVector(e.CP, out, len(gty1))
}

type scoreChunker struct {
	Engine     *heEngine
	Metrics    []categoryMetrics
	TwoParties bool
	C          int
	Slots      int

	Genes   []int
	Weights []float64
	Refs    []float64
	Gty1    []float64
	Gty2    []float64
	Gtx1    [][]float64
	Gtx2    [][]float64
}

func newScoreChunker(engine *heEngine, metrics []categoryMetrics, twoParties bool) *scoreChunker {
	c := &scoreChunker{Engine: engine, Metrics: metrics, TwoParties: twoParties, C: engine.C, Slots: engine.Slots}
	c.reset()
	return c
}

func (c *scoreChunker) reset() {
	c.Genes = make([]int, 0, c.Slots)
	c.Weights = make([]float64, 0, c.Slots)
	c.Refs = make([]float64, 0, c.Slots)
	c.Gty1 = make([]float64, 0, c.Slots)
	c.Gty2 = make([]float64, 0, c.Slots)
	c.Gtx1 = make([][]float64, c.C)
	c.Gtx2 = make([][]float64, c.C)
	for k := 0; k < c.C; k++ {
		c.Gtx1[k] = make([]float64, 0, c.Slots)
		c.Gtx2[k] = make([]float64, 0, c.Slots)
	}
}

func (c *scoreChunker) Add(r scoreRecord) error {
	c.Genes = append(c.Genes, r.Gene)
	c.Weights = append(c.Weights, r.Weight)
	c.Refs = append(c.Refs, r.Reference)
	c.Gty1 = append(c.Gty1, r.Gty1)
	c.Gty2 = append(c.Gty2, r.Gty2)
	for k := 0; k < c.C; k++ {
		c.Gtx1[k] = append(c.Gtx1[k], r.Gtx1[k])
		c.Gtx2[k] = append(c.Gtx2[k], r.Gtx2[k])
	}
	if len(c.Genes) == c.Slots {
		return c.Flush()
	}
	return nil
}

func (c *scoreChunker) Flush() error {
	if len(c.Genes) == 0 {
		return nil
	}
	got := c.Engine.evaluate(c.Gty1, c.Gtx1, c.Gty2, c.Gtx2, c.TwoParties)
	if len(got) != len(c.Genes) {
		return fmt.Errorf("CKKS decode returned %d values, want %d", len(got), len(c.Genes))
	}
	for i, v := range got {
		c.Metrics[c.Genes[i]].addHE(c.Weights[i], c.Refs[i], v)
	}
	c.reset()
	return nil
}
