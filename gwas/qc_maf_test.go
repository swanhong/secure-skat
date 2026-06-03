package gwas

import (
	"math"
	"testing"

	mpc_core "github.com/hhcho/mpc-core"
)

// TestQCMAFRareFilter exercises the SF-GWAS secure MAF filter (the division-free
// (2a-c)^2 vs c^2(2b-1)^2 test from gwas/qualcontrol.go) on a small synthetic
// "dataset" of allele counts with KNOWN frequencies, in BOTH directions:
//   - lower bound  (existing SF-GWAS): keep COMMON variants, MAF >= maf_lb
//   - upper bound  (rare flip):        keep RARE   variants, MAF <  tau
// It confirms the algebra is reusable for rare-variant selection by flipping a
// single inequality, and surfaces the a=0 (monomorphic) boundary that needs an
// explicit a>=1 guard. Runs under the 3-party harness via run_skat_tests.sh.
func TestQCMAFRareFilter(t *testing.T) {
	prot := InitProtocolForTest(t)
	if prot == nil {
		return
	}
	defer prot.SyncAndTerminate(true)

	mpcPar := prot.mpcObj
	mpcObj := prot.mpcObj[0]
	pid := mpcObj.GetPid()
	rtype := mpcObj.GetRType()
	useBoolean := mpcObj.GetBooleanShareFlag()
	const prec = 20

	// Synthetic variants: global (ALT allele count a, total observed alleles c).
	// c = 200 => 100 individuals. Mix of rare and common, plus monomorphic edge.
	variants := []struct {
		name string
		a, c int
	}{
		{"rare_singleton", 1, 200},  // AF 0.005
		{"rare_doubleton", 2, 200},  // AF 0.010
		{"rare_triple", 3, 200},     // AF 0.015
		{"common_0.20", 40, 200},    // AF 0.20
		{"common_0.30", 60, 200},    // AF 0.30
		{"common_0.50", 100, 200},   // AF 0.50
		{"monomorphic", 0, 200},     // AF 0.000 (no ALT)
	}
	m := len(variants)

	// Each party contributes its local counts as an additive secret share (the
	// SF-GWAS trick). For this test party1 holds the full global counts and the
	// others hold zero, so the secret-shared global sum equals the designed counts.
	locA := make([]int, m)
	locC := make([]int, m)
	if pid == 1 {
		for i, v := range variants {
			locA[i] = v.a
			locC[i] = v.c
		}
	}
	xSumRV := mpc_core.IntToRVec(rtype, locA)   // global ALT count a (shared)
	xCountRV := mpc_core.IntToRVec(rtype, locC) // global total alleles c (shared)

	// (2a - c), then (2a-c)^2 and c^2  (mirrors qualcontrol.go:491-495)
	twoAmC := xSumRV.Copy()
	twoAmC.MulScalar(rtype.FromInt(2))
	twoAmC.Sub(xCountRV)
	sSq := mpcPar.SSMultElemVec(twoAmC, twoAmC)     // (2a-c)^2
	cSq := mpcPar.SSMultElemVec(xCountRV, xCountRV) // c^2

	scale := rtype.FromUint64(uint64(1) << prec)

	// ---- Lower bound (existing): keep COMMON  <=>  (2a-c)^2 <= c^2 (2b-1)^2 ----
	const mafLB = 0.1
	boundLow := rtype.FromFloat64(math.Pow(2*mafLB-1.0, 2), prec)
	lo := cSq.Copy()
	lo.MulScalar(boundLow) // c^2 (2b-1)^2
	loS := sSq.Copy()
	loS.MulScalar(scale) // (2a-c)^2 * 2^prec
	lo.Sub(loS)          // c^2(2b-1)^2 - (2a-c)^2   (>0 => inside band => common)
	commonFilt := mpcPar.IsPositive(lo, useBoolean)

	// ---- Upper bound (rare flip): keep RARE  <=>  (2a-c)^2 > c^2 (1-2tau)^2 ----
	const tau = 0.05
	boundUp := rtype.FromFloat64(math.Pow(1.0-2*tau, 2), prec)
	hi := sSq.Copy()
	hi.MulScalar(scale) // (2a-c)^2 * 2^prec
	hiC := cSq.Copy()
	hiC.MulScalar(boundUp) // c^2 (1-2tau)^2
	hi.Sub(hiC)            // (2a-c)^2 - c^2(1-2tau)^2  (>0 => outside band => rare/extreme)
	rareFilt := mpcPar.IsPositive(hi, useBoolean)

	common := mpcPar.RevealSymVec(commonFilt)
	rare := mpcPar.RevealSymVec(rareFilt)

	if pid != 1 {
		return
	}

	t.Logf("variant            AF      MAF    common(>=%.2f)  rare(<%.2f)", mafLB, tau)
	failures := 0
	for i, v := range variants {
		af := float64(v.a) / float64(v.c)
		maf := math.Min(af, 1-af)
		gotCommon := common[i].Uint64() != 0
		gotRare := rare[i].Uint64() != 0
		t.Logf("%-16s  %.4f  %.4f      %v          %v", v.name, af, maf, gotCommon, gotRare)

		wantCommon := maf >= mafLB
		if gotCommon != wantCommon {
			t.Errorf("%s: common filter got %v want %v (MAF %.4f)", v.name, gotCommon, wantCommon, maf)
			failures++
		}
		// Bare MAF<tau (matches the flipped algebra). NOTE: monomorphic a=0 has MAF 0
		// and therefore passes the bare rare test — a real pipeline must AND with a>=1.
		wantRareBare := maf < tau
		if gotRare != wantRareBare {
			t.Errorf("%s: rare filter got %v want %v (MAF %.4f)", v.name, gotRare, wantRareBare, maf)
			failures++
		}
	}
	if failures == 0 {
		t.Logf("OK: secure MAF reusable for rare selection by flipping one inequality; "+
			"lower keeps common, upper keeps rare. Monomorphic (a=0) passes bare MAF<tau => needs a>=1 guard.")
	}
}
