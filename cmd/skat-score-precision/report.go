package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
	"strconv"
)

type categoryMetrics struct {
	Variants int

	RefQ    float64
	CancelQ float64
	SSLikeQ float64
	HEQ     float64
	RefL    float64
	CancelL float64
	SSLikeL float64
	HEL     float64

	WeightedTermSq       float64
	AbsWeightedTerm      float64
	AbsWeightedReference float64
	ReferenceSq          float64
	FloatErrorSq         float64
	SSLikeErrorSq        float64
	HEErrorSq            float64
	FloatMaxAbs          float64
	SSLikeMaxAbs         float64
	HEMaxAbs             float64

	KappaP99  float64
	KappaMax  float64
	HasHE     bool
	HEInvalid int
}

func (m *categoryMetrics) addPlain(r scoreRecord) {
	m.Variants++
	ws := r.Weight * r.Reference
	wc := r.Weight * r.Cancel64
	wt := r.Weight * r.SSLike
	m.RefQ += ws * ws
	m.CancelQ += wc * wc
	m.SSLikeQ += wt * wt
	m.RefL += ws
	m.CancelL += wc
	m.SSLikeL += wt
	m.WeightedTermSq += math.Pow(r.Weight*r.TermBound, 2)
	m.AbsWeightedTerm += math.Abs(r.Weight) * r.TermBound
	m.AbsWeightedReference += math.Abs(ws)
	m.ReferenceSq += r.Reference * r.Reference
	d := r.Cancel64 - r.Reference
	m.FloatErrorSq += d * d
	m.FloatMaxAbs = math.Max(m.FloatMaxAbs, math.Abs(d))
	dSS := r.SSLike - r.Reference
	m.SSLikeErrorSq += dSS * dSS
	m.SSLikeMaxAbs = math.Max(m.SSLikeMaxAbs, math.Abs(dSS))
}

func (m *categoryMetrics) setKappas(k []float64) {
	if len(k) == 0 {
		return
	}
	m.KappaP99 = percentile(k, 0.99)
	m.KappaMax = percentile(k, 1)
}

func (m *categoryMetrics) addHE(weight, reference, got float64) {
	m.HasHE = true
	if !isFinite(weight) || !isFinite(reference) || !isFinite(got) {
		m.HEInvalid++
		return
	}
	w := weight * got
	q, l := m.HEQ+w*w, m.HEL+w
	d := got - reference
	errSq := m.HEErrorSq + d*d
	if !isFinite(q) || !isFinite(l) || !isFinite(errSq) {
		m.HEInvalid++
		return
	}
	m.HEQ, m.HEL, m.HEErrorSq = q, l, errSq
	m.HEMaxAbs = math.Max(m.HEMaxAbs, math.Abs(d))
}

func combineMetrics(a, b categoryMetrics) categoryMetrics {
	return categoryMetrics{
		Variants:             a.Variants + b.Variants,
		RefQ:                 a.RefQ + b.RefQ,
		CancelQ:              a.CancelQ + b.CancelQ,
		SSLikeQ:              a.SSLikeQ + b.SSLikeQ,
		HEQ:                  a.HEQ + b.HEQ,
		RefL:                 a.RefL + b.RefL,
		CancelL:              a.CancelL + b.CancelL,
		SSLikeL:              a.SSLikeL + b.SSLikeL,
		HEL:                  a.HEL + b.HEL,
		WeightedTermSq:       a.WeightedTermSq + b.WeightedTermSq,
		AbsWeightedTerm:      a.AbsWeightedTerm + b.AbsWeightedTerm,
		AbsWeightedReference: a.AbsWeightedReference + b.AbsWeightedReference,
		ReferenceSq:          a.ReferenceSq + b.ReferenceSq,
		FloatErrorSq:         a.FloatErrorSq + b.FloatErrorSq,
		SSLikeErrorSq:        a.SSLikeErrorSq + b.SSLikeErrorSq,
		HEErrorSq:            a.HEErrorSq + b.HEErrorSq,
		FloatMaxAbs:          math.Max(a.FloatMaxAbs, b.FloatMaxAbs),
		SSLikeMaxAbs:         math.Max(a.SSLikeMaxAbs, b.SSLikeMaxAbs),
		HEMaxAbs:             math.Max(a.HEMaxAbs, b.HEMaxAbs),
		KappaP99:             math.Max(a.KappaP99, b.KappaP99),
		KappaMax:             math.Max(a.KappaMax, b.KappaMax),
		HasHE:                a.HasHE || b.HasHE,
		HEInvalid:            a.HEInvalid + b.HEInvalid,
	}
}

func relErr(got, want float64) float64 {
	if !isFinite(got) || !isFinite(want) {
		return math.Inf(1)
	}
	if math.Abs(want) <= 1e-30 {
		if math.Abs(got) <= 1e-30 {
			return 0
		}
		return math.Inf(1)
	}
	return math.Abs(got-want) / math.Abs(want)
}

func (m categoryMetrics) plainQRel() float64  { return relErr(m.CancelQ, m.RefQ) }
func (m categoryMetrics) ssLikeQRel() float64 { return relErr(m.SSLikeQ, m.RefQ) }
func (m categoryMetrics) heQRel() float64 {
	if m.HEInvalid > 0 || (m.Variants > 0 && !m.HasHE) {
		return math.Inf(1)
	}
	if !m.HasHE && m.Variants == 0 {
		return 0
	}
	if !m.HasHE {
		return math.NaN()
	}
	return relErr(m.HEQ, m.RefQ)
}

func (m categoryMetrics) lDenominator() float64 {
	return math.Max(math.Max(m.AbsWeightedReference, math.Sqrt(math.Max(0, m.RefQ))), 1e-30)
}
func (m categoryMetrics) plainLAbs() float64   { return math.Abs(m.CancelL - m.RefL) }
func (m categoryMetrics) plainLNorm() float64  { return m.plainLAbs() / m.lDenominator() }
func (m categoryMetrics) ssLikeLAbs() float64  { return math.Abs(m.SSLikeL - m.RefL) }
func (m categoryMetrics) ssLikeLNorm() float64 { return m.ssLikeLAbs() / m.lDenominator() }
func (m categoryMetrics) heLAbs() float64 {
	if m.HEInvalid > 0 || (m.Variants > 0 && !m.HasHE) {
		return math.Inf(1)
	}
	return math.Abs(m.HEL - m.RefL)
}
func (m categoryMetrics) heLNorm() float64 { return m.heLAbs() / m.lDenominator() }
func (m categoryMetrics) kappaQ() float64 {
	if m.RefQ <= 0 {
		return math.NaN()
	}
	return math.Sqrt(m.WeightedTermSq / m.RefQ)
}
func (m categoryMetrics) kappaL() float64 {
	if math.Abs(m.RefL) <= 1e-30 {
		return math.NaN()
	}
	return m.AbsWeightedTerm / math.Abs(m.RefL)
}
func (m categoryMetrics) requiredEpsilon1PctQ() float64 {
	k := m.kappaQ()
	if math.IsNaN(k) || k == 0 {
		return math.NaN()
	}
	return (math.Sqrt(1.01) - 1) / k
}
func (m categoryMetrics) heScoreNRMSE() float64 {
	if m.HEInvalid > 0 || (m.Variants > 0 && !m.HasHE) {
		return math.Inf(1)
	}
	if m.ReferenceSq <= 0 {
		if m.HEErrorSq == 0 {
			return 0
		}
		return math.Inf(1)
	}
	return math.Sqrt(m.HEErrorSq / m.ReferenceSq)
}

func isFinite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

type geneResult struct {
	Index   int
	Name    string
	Public  categoryMetrics
	Private categoryMetrics
}

func (g geneResult) Total() categoryMetrics { return combineMetrics(g.Public, g.Private) }

func formatFloat(v float64) string {
	if math.IsNaN(v) {
		return ""
	}
	return strconv.FormatFloat(v, 'e', 8, 64)
}

func writeCSV(path string, rows []geneResult, includeHE bool) error {
	f, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0o600)
	if err != nil {
		return err
	}
	defer f.Close()
	w := csv.NewWriter(f)
	header := []string{
		"gene_index", "gene_symbol", "public_variants",
		"plain_cancel_q_rel", "plain_cancel_l_abs", "plain_cancel_l_norm",
		"sslike_q_rel", "sslike_l_abs", "sslike_l_norm",
		"public_kappa_q", "private_kappa_q", "total_kappa_q",
		"required_term_rel_error_for_1pct_q", "public_cancel_kappa_p99", "public_cancel_kappa_max",
	}
	if includeHE {
		header = append(header,
			"public_ckks_q_rel", "private_ckks_q_rel", "total_ckks_q_rel",
			"public_ckks_l_norm", "private_ckks_l_norm", "total_ckks_l_abs", "total_ckks_l_norm",
			"total_ckks_score_nrmse", "total_ckks_score_max_abs", "ckks_invalid_slots")
	}
	if err := w.Write(header); err != nil {
		return err
	}
	for _, g := range rows {
		t := g.Total()
		rec := []string{
			strconv.Itoa(g.Index), g.Name, strconv.Itoa(g.Public.Variants),
			formatFloat(t.plainQRel()), formatFloat(t.plainLAbs()), formatFloat(t.plainLNorm()),
			formatFloat(t.ssLikeQRel()), formatFloat(t.ssLikeLAbs()), formatFloat(t.ssLikeLNorm()),
			formatFloat(g.Public.kappaQ()), formatFloat(g.Private.kappaQ()), formatFloat(t.kappaQ()),
			formatFloat(t.requiredEpsilon1PctQ()), formatFloat(g.Public.KappaP99), formatFloat(g.Public.KappaMax),
		}
		if includeHE {
			rec = append(rec,
				formatFloat(g.Public.heQRel()), formatFloat(g.Private.heQRel()), formatFloat(t.heQRel()),
				formatFloat(g.Public.heLNorm()), formatFloat(g.Private.heLNorm()), formatFloat(t.heLAbs()), formatFloat(t.heLNorm()),
				formatFloat(t.heScoreNRMSE()), formatFloat(t.HEMaxAbs), strconv.Itoa(t.HEInvalid))
		}
		if err := w.Write(rec); err != nil {
			return err
		}
	}
	w.Flush()
	return w.Error()
}

func metricValues(rows []geneResult, fn func(geneResult) float64) []float64 {
	out := make([]float64, 0, len(rows))
	for _, g := range rows {
		v := fn(g)
		if !math.IsNaN(v) && !math.IsInf(v, 0) {
			out = append(out, v)
		}
	}
	return out
}

func printQuantiles(label string, v []float64) {
	fmt.Printf("  %-35s p50=%9.2e  p90=%9.2e  p95=%9.2e  p99=%9.2e  max=%9.2e\n",
		label, percentile(v, .5), percentile(v, .9), percentile(v, .95), percentile(v, .99), percentile(v, 1))
}

const qZeroFloor = 1e-30

type qAudit struct {
	Total         int
	Relative      []float64
	Degenerate    int
	DegenerateAbs []float64
	Invalid       int
}

func auditQ(rows []geneResult, selectMetrics func(geneResult) categoryMetrics, got func(categoryMetrics) float64, requireHE bool) qAudit {
	a := qAudit{Total: len(rows), Relative: make([]float64, 0, len(rows))}
	for _, g := range rows {
		m := selectMetrics(g)
		if requireHE && (m.HEInvalid > 0 || (m.Variants > 0 && !m.HasHE)) {
			a.Invalid++
			continue
		}
		v := got(m)
		if !isFinite(m.RefQ) || !isFinite(v) {
			a.Invalid++
			continue
		}
		if math.Abs(m.RefQ) <= qZeroFloor {
			a.Degenerate++
			a.DegenerateAbs = append(a.DegenerateAbs, math.Abs(v-m.RefQ))
			continue
		}
		a.Relative = append(a.Relative, math.Abs(v-m.RefQ)/math.Abs(m.RefQ))
	}
	return a
}

func printQAudit(label string, a qAudit) {
	printQuantiles(label, a.Relative)
	fmt.Printf("  %-35s relative=%d/%d  zero-Q=%d  invalid=%d\n", "Q coverage (fail-closed)", len(a.Relative), a.Total, a.Degenerate, a.Invalid)
	if a.Degenerate > 0 {
		fmt.Printf("  %-35s max=%9.2e\n", "zero-Q absolute error", percentile(a.DegenerateAbs, 1))
	}
	if a.Invalid > 0 {
		fmt.Printf("  %-35s %d gene(s); do not treat this run as passing\n", "INVALID/non-finite genes", a.Invalid)
	}
}

type valueAudit struct {
	Total   int
	Values  []float64
	Invalid int
}

func auditValues(rows []geneResult, fn func(geneResult) float64) valueAudit {
	a := valueAudit{Total: len(rows), Values: make([]float64, 0, len(rows))}
	for _, g := range rows {
		v := fn(g)
		if isFinite(v) {
			a.Values = append(a.Values, v)
		} else {
			a.Invalid++
		}
	}
	return a
}

func printValueAudit(label string, a valueAudit) {
	printQuantiles(label, a.Values)
	if a.Invalid > 0 {
		fmt.Printf("  %-35s finite=%d/%d  invalid=%d\n", "metric coverage (fail-closed)", len(a.Values), a.Total, a.Invalid)
	}
}

func printThresholdCounts(label string, a qAudit) {
	within1e3, within1e2 := 0, 0
	for _, v := range a.Relative {
		if v <= 1e-3 {
			within1e3++
		}
		if v <= 1e-2 {
			within1e2++
		}
	}
	fmt.Printf("  %s Q rel <=1e-3: %d/%d; <=1e-2: %d/%d; zero-Q=%d invalid=%d\n",
		label, within1e3, a.Total, within1e2, a.Total, a.Degenerate, a.Invalid)
}

func printReport(rows []geneResult, includeHE bool, top int) {
	hasPrivate := false
	for _, g := range rows {
		if g.Private.Variants > 0 {
			hasPrivate = true
			break
		}
	}
	primaryName := "public"
	primary := func(g geneResult) categoryMetrics { return g.Public }
	if hasPrivate {
		primaryName = "total"
		primary = func(g geneResult) categoryMetrics { return g.Total() }
	}

	fmt.Println("\n=== Plain cancellation screen ===")
	printQAudit("float64 party-order Q rel. error", auditQ(rows, primary, func(m categoryMetrics) float64 { return m.CancelQ }, false))
	printValueAudit("float64 burden L normalized error", auditValues(rows, func(g geneResult) float64 { return primary(g).plainLNorm() }))
	printValueAudit("float64 burden L absolute error", auditValues(rows, func(g geneResult) float64 { return primary(g).plainLAbs() }))
	printQAudit("SS-like fixed-point Q rel. error", auditQ(rows, primary, func(m categoryMetrics) float64 { return m.SSLikeQ }, false))
	printValueAudit("SS-like burden L normalized error", auditValues(rows, func(g geneResult) float64 { return primary(g).ssLikeLNorm() }))
	printQuantiles("Q cancellation amplification kappa", metricValues(rows, func(g geneResult) float64 { return primary(g).kappaQ() }))
	required := metricValues(rows, func(g geneResult) float64 { return primary(g).requiredEpsilon1PctQ() })
	fmt.Printf("  %-35s min=%9.2e  p01=%9.2e  p05=%9.2e  p50=%9.2e\n",
		"term rel. error allowed for 1% Q", percentile(required, 0), percentile(required, .01), percentile(required, .05), percentile(required, .5))
	fmt.Printf("  scope used above: %s score\n", primaryName)

	if includeHE {
		fmt.Println("\n=== Single-process packed CKKS vs stable plaintext ===")
		publicQ := auditQ(rows, func(g geneResult) categoryMetrics { return g.Public }, func(m categoryMetrics) float64 { return m.HEQ }, true)
		printQAudit("public Q relative error", publicQ)
		printValueAudit("public burden L normalized error", auditValues(rows, func(g geneResult) float64 { return g.Public.heLNorm() }))
		printThresholdCounts("public", publicQ)
		if hasPrivate {
			privateQ := auditQ(rows, func(g geneResult) categoryMetrics { return g.Private }, func(m categoryMetrics) float64 { return m.HEQ }, true)
			totalQ := auditQ(rows, func(g geneResult) categoryMetrics { return g.Total() }, func(m categoryMetrics) float64 { return m.HEQ }, true)
			printQAudit("private Q relative error", privateQ)
			printQAudit("total Q relative error", totalQ)
			printValueAudit("total burden L normalized error", auditValues(rows, func(g geneResult) float64 { return g.Total().heLNorm() }))
			printThresholdCounts("private", privateQ)
			printThresholdCounts("total", totalQ)
		}
		printValueAudit(primaryName+" score normalized RMSE", auditValues(rows, func(g geneResult) float64 { return primary(g).heScoreNRMSE() }))
		printValueAudit(primaryName+" max score absolute error", auditValues(rows, func(g geneResult) float64 { return primary(g).HEMaxAbs }))
	}

	if top <= 0 {
		return
	}
	ordered := append([]geneResult(nil), rows...)
	if includeHE {
		sort.Slice(ordered, func(i, j int) bool {
			return nanLow(primary(ordered[i]).heQRel()) > nanLow(primary(ordered[j]).heQRel())
		})
		fmt.Printf("\n=== Worst %d genes by %s CKKS Q relative error ===\n", minInt(top, len(ordered)), primaryName)
		fmt.Printf("  %6s  %-20s  %11s  %11s  %11s  %7s\n", "index", "gene", "Q_rel", "L_norm", "kappa_Q", "invalid")
		for _, g := range ordered[:minInt(top, len(ordered))] {
			m := primary(g)
			fmt.Printf("  %6d  %-20s  %11.3e  %11.3e  %11.3e  %7d\n",
				g.Index, g.Name, m.heQRel(), m.heLNorm(), m.kappaQ(), m.HEInvalid)
		}
		orderedL := append([]geneResult(nil), rows...)
		sort.Slice(orderedL, func(i, j int) bool {
			return nanLow(primary(orderedL[i]).heLNorm()) > nanLow(primary(orderedL[j]).heLNorm())
		})
		fmt.Printf("\n=== Worst %d genes by %s CKKS burden L normalized error ===\n", minInt(top, len(orderedL)), primaryName)
		fmt.Printf("  %6s  %-20s  %11s  %11s  %11s\n", "index", "gene", "L_norm", "L_abs", "Q_rel")
		for _, g := range orderedL[:minInt(top, len(orderedL))] {
			m := primary(g)
			fmt.Printf("  %6d  %-20s  %11.3e  %11.3e  %11.3e\n", g.Index, g.Name, m.heLNorm(), m.heLAbs(), m.heQRel())
		}
	} else {
		sort.Slice(ordered, func(i, j int) bool {
			return nanLow(primary(ordered[i]).kappaQ()) > nanLow(primary(ordered[j]).kappaQ())
		})
		fmt.Printf("\n=== Worst %d genes by cancellation amplification ===\n", minInt(top, len(ordered)))
		fmt.Printf("  %6s  %-20s  %11s  %14s  %11s  %11s\n", "index", "gene", "kappa_Q", "required_eps", "L_norm", "SSlike_Q")
		for _, g := range ordered[:minInt(top, len(ordered))] {
			m := primary(g)
			fmt.Printf("  %6d  %-20s  %11.3e  %14.3e  %11.3e  %11.3e\n",
				g.Index, g.Name, m.kappaQ(), m.requiredEpsilon1PctQ(), m.plainLNorm(), m.ssLikeQRel())
		}
	}
}

func nanLow(v float64) float64 {
	if math.IsNaN(v) {
		return -1
	}
	return v
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}
