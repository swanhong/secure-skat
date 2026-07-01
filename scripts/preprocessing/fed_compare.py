#!/usr/bin/env python3
"""Compare the secure skat_fed per-gene Q against the plaintext federated_Q_from_blocks on the
SAME fed_prep blocks. Run on the workbench after the secure run.

    python3 fed_compare.py        # reads $FED_OUT (default ~/fed_prep_out)

Loads the int8 genotype blocks + cov(5 PCs)+pheno written by fed_prep, recomputes per-gene Q in
plaintext (PART A = A + B-aligned, PART B = B private), and diffs against out/party2/skat_fed_out.txt.
A CKKS match is ~1e-3; we flag rel > 1e-2.
"""
import json
import os

import numpy as np

from skat_plain_local import federated_Q_from_blocks

OUT = os.path.expanduser(os.environ.get("FED_OUT", "~/fed_prep_out"))


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


ng = json.load(open(f"{OUT}/manifest.json"))["n_genes"]
nA = len(np.loadtxt(f"{OUT}/A/pheno.txt"))
nB = len(np.loadtxt(f"{OUT}/B/pheno.txt"))
XA, yA = load_xy("A", nA)
XB, yB = load_xy("B", nB)
Qplain = federated_Q_from_blocks(
    load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
    load_blocks("B", "priv", nB, ng), XA, yA, XB, yB)
Qsec = np.atleast_1d(np.loadtxt(f"{OUT}/out/party2/skat_fed_out.txt"))

print(f"  nA={nA} nB={nB} genes={ng}")
print(f"{'gene':>4} {'Q_secure':>16} {'Q_plain':>16} {'rel':>10}  status")
ok = True
for b in range(ng):
    d = abs(Qsec[b] - Qplain[b])
    rel = d / max(abs(Qplain[b]), 1e-9)
    # np.isclose semantics: atol handles near-zero Q (e.g. monomorphic 1-variant gene, Q~0,
    # where rel is meaningless); rtol handles the CKKS precision at scale.
    match = d <= 1.0 + 1e-2 * abs(Qplain[b])
    ok &= match
    print(f"{b:>4} {Qsec[b]:>16.4f} {Qplain[b]:>16.4f} {rel:>10.2e}  {'ok' if match else 'MISMATCH'}")
print("MATCH (secure == plaintext)" if ok else "MISMATCH -- investigate")
