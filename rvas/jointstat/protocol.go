package jointstat

import (
	"fmt"

	core "github.com/hhcho/mpc-core"
	"github.com/hhcho/sfgwas/mpc"
	"github.com/hhcho/sfgwas/rvas/protocol"
)

// Counts and flags are exact unscaled field elements. Compare 2*x with 2*t-1
// to avoid equality ambiguity without fixed-point scaling or CKKS rounding.
type engine struct{ obj *mpc.MPC }

func (e engine) zeros(n int) core.RVec        { return core.InitRVec(e.obj.GetRType().Zero(), n) }
func (e engine) share(values []int) core.RVec { return protocol.ShareSumInt(e.obj, values) }
func (e engine) mul(a, b core.RVec) core.RVec {
	if len(a) == 0 {
		return e.zeros(0)
	}
	return e.obj.SSMultElemVec(a, b)
}
func (e engine) ge(a core.RVec, threshold int) core.RVec {
	if len(a) == 0 {
		return e.zeros(0)
	}
	scaled := a.Copy()
	scaled.MulScalar(e.obj.GetRType().FromInt(2))
	return e.obj.NotLessThanPublic(scaled, e.obj.GetRType().FromInt(2*threshold-1), e.obj.GetBooleanShareFlag())
}
func (e engine) sum(a core.RVec) core.RElem { return protocol.SumShares(e.obj.GetRType(), a) }

type layout struct {
	Genes   []string
	Keys    []string
	Mapping [][]int
}

func publicLayout(p *prepared) layout {
	l := layout{Mapping: p.PublicGenes}
	for _, g := range p.Genes {
		l.Genes = append(l.Genes, g.ID)
	}
	for _, v := range p.Public {
		l.Keys = append(l.Keys, v.Key)
	}
	return l
}

type dimensions struct {
	NA, NB, CohortN int
	Layout          layout
}
type row struct {
	Ancestry                    string
	MAF                         float64
	Annotation                  string
	Chromosome                  int
	Panel, Series, Lower, Upper string
	Count                       int
}
type calculation struct {
	Labels []row
	Shares core.RVec
}

type statistics struct {
	kAlt, kRef, minorA, minorB       core.RVec
	both, onlyA, onlyB               core.RVec
	privateMinor                     []int
	privateTotal                     core.RElem
	carriersA, carriersB, macA, macB core.RVec
}

// compute reveals nothing. Counts and genotype indicators are shared once per
// chromosome; all MAFs use the same shares. Private identities stay at B.
func compute(obj *mpc.MPC, p *prepared, site string, d dimensions, mafs []float64, thresholds, bins []int, batchCells int) ([]calculation, error) {
	e := engine{obj}
	pid, n, ng := obj.GetPid(), len(d.Layout.Keys), len(d.Layout.Genes)
	states := make([]statistics, len(mafs))
	local := make([]int, 5*n+len(mafs)*(ng+1))
	if pid > 0 {
		for j, v := range p.Public {
			local[j] = v.CohortAC
			offset := n
			if pid == 2 {
				offset = 3 * n
			}
			local[offset+j], local[offset+n+j] = v.AC, v.Called
		}
	}
	for m, maf := range mafs {
		st := &states[m]
		st.carriersA, st.carriersB, st.macA, st.macB = e.zeros(ng), e.zeros(ng), e.zeros(ng), e.zeros(ng)
		if pid != 2 {
			continue
		}
		st.privateMinor = make([]int, len(p.Private))
		offset := 5*n + m*(ng+1)
		for j, v := range p.Private {
			if !protocol.GlobalMAFKeepPlain(float64(v.CohortAC), d.CohortN, maf) {
				continue
			}
			ac := v.AC
			if v.AC > d.NA+d.NB {
				ac = 2*v.Called - ac
			}
			st.privateMinor[j] = ac
			if ac > 0 {
				local[offset]++
			}
		}
		for g, indices := range p.PrivateGenes {
			for _, j := range indices {
				local[offset+1+g] += st.privateMinor[j]
			}
		}
	}
	shared := e.share(local)
	cohort, acA, calledA := shared[:n], shared[n:2*n], shared[2*n:3*n]
	acB, calledB := shared[3*n:4*n], shared[4*n:5*n]
	jointAC := acA.Copy()
	jointAC.Add(acB)
	// Main's current MPC comparison returns zero at an exact tie.
	flip := e.ge(jointAC, d.NA+d.NB+1)
	for m, maf := range mafs {
		st := &states[m]
		keep := protocol.GlobalMAFKeep(obj, cohort, d.CohortN, maf, 0)
		st.kRef = e.mul(keep, flip)
		st.kAlt = keep.Copy()
		st.kAlt.Sub(st.kRef)
		minor := func(ac, called core.RVec) core.RVec {
			ref := called.Copy()
			ref.MulScalar(obj.GetRType().FromInt(2))
			ref.Sub(ac)
			result := e.mul(st.kAlt, ac)
			result.Add(e.mul(st.kRef, ref))
			return result
		}
		st.minorA, st.minorB = minor(acA, calledA), minor(acB, calledB)
		observedA, observedB := e.ge(st.minorA, 1), e.ge(st.minorB, 1)
		st.both = e.mul(observedA, observedB)
		st.onlyA = observedA.Copy()
		st.onlyA.Sub(st.both)
		st.onlyB = observedB.Copy()
		st.onlyB.Sub(st.both)
		offset := 5*n + m*(ng+1)
		st.privateTotal = shared[offset]
		for g, indices := range d.Layout.Mapping {
			for _, j := range indices {
				st.macA[g] = st.macA[g].Add(st.minorA[j])
				st.macB[g] = st.macB[g].Add(st.minorB[j])
			}
			st.macB[g] = st.macB[g].Add(shared[offset+1+g])
		}
	}
	if err := carriers(e, p, site, d, states, batchCells); err != nil {
		return nil, err
	}
	results := make([]calculation, len(mafs))
	for i, st := range states {
		results[i] = aggregate(e, st, thresholds, bins)
	}
	return results, nil
}

func carriers(e engine, p *prepared, site string, d dimensions, states []statistics, batchCells int) error {
	ng, nm := len(d.Layout.Genes), len(states)
	// One matrix per gene: [ALT-present, REF-present] x [MAF weights].
	// SSMultMats batches their openings without a huge block-diagonal matrix.
	var geneMatrices []int
	var weights []core.RMat
	cellsPerRow := ng * nm
	for g, indices := range d.Layout.Mapping {
		cellsPerRow += 2 * len(indices)
		if len(indices) == 0 {
			continue
		}
		geneMatrices = append(geneMatrices, g)
		w := core.InitRMat(e.obj.GetRType().Zero(), 2*len(indices), nm)
		for col, j := range indices {
			for m, st := range states {
				w[2*col][m], w[2*col+1][m] = st.kAlt[j], st.kRef[j]
			}
		}
		weights = append(weights, w)
	}
	batchRows := max(1, min(1024, batchCells/max(1, cellsPerRow)))
	pid := e.obj.GetPid()
	for start := 0; start < d.NA+d.NB; start += batchRows {
		end := min(start+batchRows, d.NA+d.NB)
		width := end - start
		// Fixed public layout includes private-presence slots even at empty genes.
		local := make([]int, width*cellsPerRow)
		privateSize := ng * width * nm
		offset := privateSize
		for g, indices := range d.Layout.Mapping {
			if pid > 0 {
				base, limit := 0, d.NA
				if pid == 2 {
					base, limit = d.NA, d.NA+d.NB
				}
				lo, hi := max(start, base), min(end, limit)
				if lo < hi {
					pub, err := readBlockRows(site, g, p.Samples, len(indices), false, lo-base, hi-base)
					if err != nil {
						return err
					}
					for i := lo; i < hi; i++ {
						for col := range indices {
							at := (i-lo)*pub.Columns + col
							if !pub.called(at) {
								continue
							}
							x := offset + (i-start)*2*len(indices) + 2*col
							if pub.Dosage[at] > 0 {
								local[x] = 1
							}
							if pub.Dosage[at] < 2 {
								local[x+1] = 1
							}
						}
					}
					if pid == 2 {
						priv, err := readBlockRows(site, g, p.Samples, len(p.PrivateGenes[g]), true, lo-base, hi-base)
						if err != nil {
							return err
						}
						for i := lo; i < hi; i++ {
							for col, j := range p.PrivateGenes[g] {
								at := (i-lo)*priv.Columns + col
								flip := p.Private[j].AC > d.NA+d.NB
								if !priv.called(at) || (!flip && priv.Dosage[at] == 0) || (flip && priv.Dosage[at] == 2) {
									continue
								}
								for m, st := range states {
									if st.privateMinor[j] > 0 {
										local[(g*width+i-start)*nm+m] = 1
									}
								}
							}
						}
					}
				}
			}
			offset += width * 2 * len(indices)
		}
		shared := e.share(local)
		burden := shared[:privateSize]
		var indicators []core.RMat
		offset = privateSize
		for _, g := range geneMatrices {
			cols := 2 * len(d.Layout.Mapping[g])
			matrix := make(core.RMat, width)
			for i := range matrix {
				matrix[i] = shared[offset+i*cols : offset+(i+1)*cols]
			}
			indicators = append(indicators, matrix)
			offset += width * cols
		}
		if len(indicators) > 0 {
			products := e.obj.SSMultMats(indicators, weights)
			for j, g := range geneMatrices {
				for i := 0; i < width; i++ {
					for m := range states {
						at := (g*width+i)*nm + m
						burden[at] = burden[at].Add(products[j][i][m])
					}
				}
			}
		}
		present := e.ge(burden, 1)
		for g := 0; g < ng; g++ {
			for i := 0; i < width; i++ {
				for m, st := range states {
					counts := st.carriersA
					if start+i >= d.NA {
						counts = st.carriersB
					}
					counts[g] = counts[g].Add(present[(g*width+i)*nm+m])
				}
			}
		}
	}
	return nil
}

func aggregate(e engine, st statistics, thresholds, bins []int) calculation {
	carriersA, carriersB, macA, macB := st.carriersA, st.carriersB, st.macA, st.macB
	jointCarrier := carriersA.Copy()
	jointCarrier.Add(carriersB)
	jointMAC := macA.Copy()
	jointMAC.Add(macB)
	result := calculation{}
	appendCount := func(panel, series, lo, hi string, count core.RElem) {
		result.Labels = append(result.Labels, row{Panel: panel, Series: series, Lower: lo, Upper: hi})
		result.Shares = append(result.Shares, count)
	}
	for _, t := range thresholds {
		a, b, j := e.ge(carriersA, t), e.ge(carriersB, t), e.ge(jointCarrier, t)
		any := a.Copy()
		any.Add(b)
		any.Sub(e.mul(a, b))
		for i, values := range []core.RVec{a, b, any, j} {
			appendCount("A", []string{"MVP", "AoU", "local_any", "joint"}[i], fmt.Sprint(t), "NA", e.sum(values))
		}
	}
	eligible := e.mul(e.ge(carriersA, 20), e.ge(carriersB, 20))
	// Every gene appears in exactly one bin, including a final unbounded bin.
	for s, values := range []core.RVec{macA, macB, jointMAC} {
		for i, lo := range bins {
			selected := e.ge(values, lo)
			upper := "Inf"
			if i+1 < len(bins) {
				hi := e.ge(values, bins[i+1])
				selected.Sub(hi)
				upper = fmt.Sprint(bins[i+1])
			}
			appendCount("B", []string{"MVP", "AoU", "joint"}[s], fmt.Sprint(lo), upper, e.sum(e.mul(eligible, selected)))
		}
	}
	appendCount("C", "shared", "NA", "NA", e.sum(st.both))
	appendCount("C", "MVP_only", "NA", "NA", e.sum(st.onlyA))
	appendCount("C", "AoU_only", "NA", "NA", e.sum(st.onlyB).Add(st.privateTotal))
	return result
}
