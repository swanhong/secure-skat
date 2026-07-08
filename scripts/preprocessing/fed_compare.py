#!/usr/bin/env python3
"""Compare the secure skat_fed per-gene SKAT (Q) and Burden against the plaintext
federated_skat_burden_from_blocks on the SAME fed_prep blocks. Run on the workbench after the secure run.

    python3 fed_compare.py        # reads $FED_OUT (default ~/fed_prep_out)

Loads the int8 genotype blocks + cov(5 PCs)+pheno written by fed_prep, recomputes per-gene SKAT and
Burden in plaintext (PART A = A + B-aligned, PART B = B private), and diffs against the secure
outputs out/party2/skat_fed_out.txt (SKAT) and skat_fed_burden_out.txt (Burden). A CKKS match is ~1e-3;
the per-gene 'within' column shows O/X against THRES (default 0.01).
"""
import json
import os
import time

import numpy as np

from skat_plain_local import federated_skat_burden_from_blocks

OUT = os.path.expanduser(os.environ.get("FED_OUT", "~/fed_prep_out"))
RIDGE_REL = 1e-6  # match secure computeBetaHatEnc Tikhonov ridge (gwas/skat.go ridgeRel)
THRES = 0.01  # per-gene relative-gap threshold for the 'within' (O/X) column


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


def compare(name, secure, plain, ng, thres=THRES):
    """Per-gene secure-vs-plaintext table for one statistic. The 'within' column is O when the
    relative gap is <= thres (X otherwise). Returns the count of genes within thres."""
    secure, plain = np.atleast_1d(secure), np.asarray(plain, float)
    print(f"\n=== {name}  (threshold rel = {thres}) ===")
    print(f"{'gene':>4} {'secure':>16} {'plain':>16} {'rel':>10}  within")
    n_within = 0
    for b in range(ng):
        d = abs(secure[b] - plain[b])
        rel = d / max(abs(plain[b]), 1e-9)
        # atol 1 covers near-zero genes (rel meaningless); rtol=thres covers CKKS precision at scale.
        within = d <= 1.0 + thres * abs(plain[b])
        n_within += within
        print(f"{b:>4} {secure[b]:>16.4f} {plain[b]:>16.4f} {rel:>10.2e}  {'O' if within else 'X'}")
    ss_tot = float(np.sum((plain[:ng] - plain[:ng].mean()) ** 2))
    r2 = 1 - float(np.sum((secure[:ng] - plain[:ng]) ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    maxrel = max(abs(secure[b] - plain[b]) / max(abs(plain[b]), 1e-9) for b in range(ng))
    print(f"R^2 (secure vs plaintext) = {r2:.6f}   |   max rel = {maxrel:.2e}   |   within: {n_within}/{ng}")
    return n_within


ng = json.load(open(f"{OUT}/manifest.json"))["n_genes"]
nA = len(np.loadtxt(f"{OUT}/A/pheno.txt"))
nB = len(np.loadtxt(f"{OUT}/B/pheno.txt"))
XA, yA = load_xy("A", nA)
XB, yB = load_xy("B", nB)
t0 = time.perf_counter()
Splain, Bplain = federated_skat_burden_from_blocks(
    load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
    load_blocks("B", "priv", nB, ng), XA, yA, XB, yB, ridge_rel=RIDGE_REL)
print(f"  plaintext federated SKAT+Burden: {time.perf_counter() - t0:.2f}s")
Ssec = np.loadtxt(f"{OUT}/out/party2/skat_fed_out.txt")
Bsec = np.loadtxt(f"{OUT}/out/party2/skat_fed_burden_out.txt")

print(f"  nA={nA} nB={nB} genes={ng}")
skat_within = compare("SKAT", Ssec, Splain, ng)
burden_within = compare("Burden", Bsec, Bplain, ng)
print(f"\nwithin threshold ({THRES}):  SKAT {skat_within}/{ng},  Burden {burden_within}/{ng}")
