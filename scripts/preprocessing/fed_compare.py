#!/usr/bin/env python3
"""Compare the secure skat_fed per-gene SKAT (Q) and Burden p-value against the plaintext
federated_skat_burden_from_blocks on the SAME fed_prep blocks. Run on the workbench after the secure run.

    python3 fed_compare.py        # reads $FED_OUT (default ~/fed_prep_out)

Loads the int8 genotype blocks + cov(5 PCs)+pheno written by fed_prep, recomputes per-gene SKAT and
Burden p-value in plaintext (PART A = A + B-aligned, PART B = B private), and diffs against the secure
outputs out/party2/skat_fed_out.txt (SKAT) and skat_fed_burden_p_out.txt (Burden p). The secure run
never releases the Burden statistic or zᵀPz, so there is nothing to diff for those. A CKKS match is
~1e-3; the per-gene 'within' column shows O/X against THRES (default 0.01).
"""
import json
import math
import os
import time

import numpy as np

from skat_plain_local import federated_skat_burden_from_blocks, federated_skat_p_from_blocks

OUT = os.path.expanduser(os.environ.get("FED_OUT", "~/fed_prep_out"))
RIDGE_REL = 1e-6  # match secure computeBetaHatEnc Tikhonov ridge (gwas/skat.go ridgeRel)
THRES = 0.01  # per-gene relative-gap threshold for the 'within' (O/X) column


def emit_timing(scope, milliseconds, kind="phase", status="done", count=1):
    """Emit one machine-readable timing record for timing_summary.py."""
    print(f"[timing] scope={scope} parent=compare party=driver kind={kind} "
          f"status={status} milliseconds={milliseconds:.3f} count={count}", flush=True)


def load_blocks(sub, kind, n, ng):
    blocks = []
    for b in range(ng):
        g = np.fromfile(f"{OUT}/{sub}/{kind}.{b}.bin", dtype=np.int8)
        blocks.append(g.reshape(n, len(g) // n))  # m = size/n; empty -> (n,0)
    return blocks


def load_xy(sub, n):
    y = np.loadtxt(f"{OUT}/{sub}/pheno.txt")
    cov = np.loadtxt(f"{OUT}/{sub}/cov.txt")
    if cov.ndim == 1:
        cov = cov.reshape(n, -1)
    return np.column_stack([np.ones(n), cov]), y  # X = [1 | PCs]


def compare(name, secure, plain, ng, thres=THRES, atol=1.0):
    """Per-gene secure-vs-plaintext table for one statistic. The 'within' column is O when the
    relative gap is <= thres (X otherwise). atol guards near-zero genes for large-magnitude stats
    (SKAT/Burden); pass atol=0 for O(1) stats like the p-value. Returns the count within thres."""
    secure, plain = np.atleast_1d(secure), np.asarray(plain, float)
    print(f"\n=== {name}  (threshold rel = {thres}) ===")
    print(f"{'gene':>4} {'secure':>16} {'plain':>16} {'rel':>10}  within")
    n_within = 0
    for b in range(ng):
        d = abs(secure[b] - plain[b])
        rel = d / max(abs(plain[b]), 1e-9)
        # atol guards near-zero genes (rel meaningless); rtol=thres covers CKKS precision at scale.
        within = d <= atol + thres * abs(plain[b])
        n_within += within
        print(f"{b:>4} {secure[b]:>16.4f} {plain[b]:>16.4f} {rel:>10.2e}  {'O' if within else 'X'}")
    ss_tot = float(np.sum((plain[:ng] - plain[:ng].mean()) ** 2))
    r2 = 1 - float(np.sum((secure[:ng] - plain[:ng]) ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    maxrel = max(abs(secure[b] - plain[b]) / max(abs(plain[b]), 1e-9) for b in range(ng))
    print(f"R^2 (secure vs plaintext) = {r2:.6f}   |   max rel = {maxrel:.2e}   |   within: {n_within}/{ng}")
    return n_within


compare_started = time.perf_counter()
phase_started = time.perf_counter()
ng = json.load(open(f"{OUT}/manifest.json"))["n_genes"]
nA = len(np.loadtxt(f"{OUT}/A/pheno.txt"))
nB = len(np.loadtxt(f"{OUT}/B/pheno.txt"))
XA, yA = load_xy("A", nA)
XB, yB = load_xy("B", nB)
emit_timing("compare.load_inputs", 1000.0 * (time.perf_counter() - phase_started))
t0 = time.perf_counter()
Splain, _, Pplain = federated_skat_burden_from_blocks(
    load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
    load_blocks("B", "priv", nB, ng), XA, yA, XB, yB, ridge_rel=RIDGE_REL)
plain_skat_burden_ms = 1000.0 * (time.perf_counter() - t0)
print(f"  plaintext federated SKAT+Burden-p: {plain_skat_burden_ms / 1000.0:.2f}s")
emit_timing("compare.plain_skat_burden", plain_skat_burden_ms)
# Reveal set depends on the mode: with skat_pvalue_probes=0 the run reveals the SKAT statistic Q
# (skat_fed_out.txt) + Burden p; with skat_pvalue_probes>0 Q is withheld (only the WH pivot z leaves,
# so skat_fed_out.txt is absent) and the SKAT p-value section below runs instead. Burden always
# reveals only √(T/2) (→ p); the Burden statistic and zᵀPz stay secret-shared.
skat_q_file = f"{OUT}/out/party2/skat_fed_out.txt"
phase_started = time.perf_counter()
Ssec = np.atleast_1d(np.loadtxt(skat_q_file)) if os.path.exists(skat_q_file) else None
Psec = np.atleast_1d(np.loadtxt(f"{OUT}/out/party2/skat_fed_burden_p_out.txt"))
emit_timing("compare.load_secure_outputs", 1000.0 * (time.perf_counter() - phase_started))

def _erfcinv(y):
    """Inverse complementary error function via bisection (erfc is monotone-decreasing on [0,40],
    covering p down to ~1e-300). Pure numpy; recovers √(T/2) from the revealed Burden p."""
    y = float(min(max(y, 1e-300), 1.0))
    lo, hi = 0.0, 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if math.erfc(mid) > y:  # erfc(mid) too large ⇒ mid too small
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print(f"  nA={nA} nB={nB} genes={ng}")
if Ssec is not None:  # SKAT statistic Q — only revealed when skat_pvalue_probes=0 (else Q is hidden)
    skat_within = compare("SKAT (Q)", Ssec, Splain, ng)
else:
    skat_within = None
    print("\n  SKAT statistic Q not released (skat_pvalue_probes>0 → only the WH pivot z leaves)")
burdenp_within = compare("Burden p-value", Psec, Pplain, ng, atol=0.0)
# Burden T = 2·(√(T/2))² = 2·erfcinv(p)²: the χ²₁ statistic behind the Burden p. Derived from the
# revealed p (T ↔ p bijective), so no new information leaves; shown for interpretability.
Tsec = np.array([2 * _erfcinv(p) ** 2 for p in np.atleast_1d(Psec)])
Tplain = np.array([2 * _erfcinv(p) ** 2 for p in np.atleast_1d(Pplain)])
# atol=0.02: T=2·erfcinv(p)² amplifies the relative error at high p (T→0), so a small-T gene can miss
# a pure-rel threshold while its authoritative Burden p is within — the absolute floor absorbs that.
burdenT_within = compare("Burden T (=2·erfcinv(p)²)", Tsec, Tplain, ng, atol=0.02)
summary = f"\nwithin threshold ({THRES}):  "
summary += f"SKAT Q {skat_within}/{ng},  " if skat_within is not None else ""
summary += f"Burden p {burdenp_within}/{ng},  Burden T {burdenT_within}/{ng}"

# SKAT p-value (only if the secure run was configured with skat_pvalue_probes > 0). The secure output
# is the Wilson-Hilferty screening p; the plaintext oracle also computes the Liu and Davies references
# so the WH approximation can be judged against them (WH vs Liu, WH vs Davies).
SkatPsec = SkatPplain = None
skat_p_file = f"{OUT}/out/party2/skat_fed_skat_p_out.txt"
if os.path.exists(skat_p_file):
    skat_p_started = time.perf_counter()
    aB = load_blocks("A", "geno", nA, ng)
    bB = load_blocks("B", "geno", nB, ng)
    pB = load_blocks("B", "priv", nB, ng)
    WHplain, LiuPlain, DaviesPlain = (np.asarray(x, float)
                                      for x in federated_skat_p_from_blocks(aB, bB, pB, XA, yA, XB, yB, ridge_rel=RIDGE_REL))
    emit_timing("compare.plain_skat_p_refs", 1000.0 * (time.perf_counter() - skat_p_started))
    SkatPsec = np.atleast_1d(np.loadtxt(skat_p_file))  # secure Wilson-Hilferty p
    SkatPplain = WHplain  # plain WH — the CSV/scatter reference for the secure output

    print("\n=== SKAT p-value: secure WH vs plaintext references (Liu, Davies) ===")
    print(f"{'gene':>4} {'secure_WH':>12} {'plain_WH':>12} {'Liu':>12} {'Davies':>12}")
    for b in range(ng):
        print(f"{b:>4} {SkatPsec[b]:>12.4e} {WHplain[b]:>12.4e} {LiuPlain[b]:>12.4e} {DaviesPlain[b]:>12.4e}")

    def _r2(a, ref):
        ss = float(np.sum((ref - ref.mean()) ** 2))
        return 1 - float(np.sum((a - ref) ** 2)) / ss if ss > 0 else float("nan")

    def _maxrel(a, ref):
        return max(abs(a[i] - ref[i]) / max(abs(ref[i]), 1e-12) for i in range(len(ref)))

    def _mae(a, ref):
        return float(np.mean(np.abs(a - ref)))

    print(f"\n  {'comparison':<24} {'R^2':>10} {'max rel':>10} {'mean abs err':>13}")
    for name, ref in [("secure-WH vs plain-WH", WHplain), ("WH vs Liu", LiuPlain), ("WH vs Davies", DaviesPlain)]:
        print(f"  {name:<24} {_r2(SkatPsec, ref):>10.6f} {_maxrel(SkatPsec, ref):>10.2e} {_mae(SkatPsec, ref):>13.2e}")
    summary += f",  SKAT p (WH vs Liu R^2={_r2(SkatPsec, LiuPlain):.4f}, vs Davies R^2={_r2(SkatPsec, DaviesPlain):.4f})"
else:
    emit_timing("compare.plain_skat_p_refs", 0.0, status="skipped", count=0)
print(summary)


def gene_positions(ng):
    """Per-gene (chrom, median pos) from block_sizes.txt (m per gene) + pos.txt (chrom<TAB>pos per
    public-list variant, gene-block order). Public genomic position only — no sample-level data."""
    sizes = [int(x) for x in open(f"{OUT}/block_sizes.txt")]
    rows = [ln.split() for ln in open(f"{OUT}/pos.txt")]
    chrom, pos, i = [], [], 0
    for g in range(ng):
        chunk = rows[i:i + sizes[g]]
        i += sizes[g]
        chrom.append(int(chunk[0][0]))
        pos.append(int(np.median([int(r[1]) for r in chunk])))
    return chrom, pos


# FED_CSV=1: dump per-gene aggregate results (positions + p-values only) for fed_plot.py. The SKAT
# p-value columns appear only when the secure run computed it (skat_pvalue_probes > 0).
if os.environ.get("FED_CSV"):
    csv_started = time.perf_counter()
    chrom, pos = gene_positions(ng)
    has_skatp = SkatPsec is not None
    csv_path = f"{OUT}/fed_results.csv"
    with open(csv_path, "w") as f:
        header = "gene,chrom,pos,burden_p_secure,burden_p_plain"
        if has_skatp:
            header += ",skat_p_secure,skat_p_plain"
        f.write(header + "\n")
        for b in range(ng):
            row = f"{b},{chrom[b]},{pos[b]},{Psec[b]:.6e},{Pplain[b]:.6e}"
            if has_skatp:
                row += f",{SkatPsec[b]:.6e},{SkatPplain[b]:.6e}"
            f.write(row + "\n")
    print(f"  FED_CSV: wrote {csv_path}  ({'burden+skat p' if has_skatp else 'burden p'})")
    emit_timing("compare.result_csv_write", 1000.0 * (time.perf_counter() - csv_started))
else:
    emit_timing("compare.result_csv_write", 0.0, status="skipped", count=0)

emit_timing("compare.total", 1000.0 * (time.perf_counter() - compare_started), kind="total")
