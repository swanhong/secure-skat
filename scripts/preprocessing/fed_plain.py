#!/usr/bin/env python3
"""Plaintext federated SKAT/Burden on fed_prep output — NO secure run needed. Runs all genes and
prints/saves per-gene SKAT Q, SKAT p (Wilson-Hilferty / Liu / Davies) and Burden p, using the same
plaintext oracle (skat_plain_local.py) that fed_compare.py uses for its "plain" columns.

    FED_OUT=~/runs/out260713034114 python3 fed_plain.py
    python3 fed_plain.py ~/runs/out260713034114        # dir as arg also works
"""
import json
import os
import sys
import time

import numpy as np

# import the plaintext oracle (this file lives next to skat_plain_local.py; fallback = repo copy)
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.expanduser("~/secure-skat/scripts/preprocessing")):
    if os.path.exists(os.path.join(p, "skat_plain_local.py")):
        sys.path.insert(0, p)
        break
from skat_plain_local import federated_skat_burden_from_blocks, federated_skat_p_from_blocks

OUT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FED_OUT", "~/fed_prep_out"))
RIDGE_REL = 1e-6  # match secure computeBetaHatEnc Tikhonov ridge

ng = json.load(open(f"{OUT}/manifest.json"))["n_genes"]
nA = len(np.loadtxt(f"{OUT}/A/pheno.txt"))
nB = len(np.loadtxt(f"{OUT}/B/pheno.txt"))


def load_xy(sub, n):
    cov = np.loadtxt(f"{OUT}/{sub}/cov.txt")
    if cov.ndim == 1:
        cov = cov.reshape(n, -1)
    return np.column_stack([np.ones(n), cov]), np.loadtxt(f"{OUT}/{sub}/pheno.txt")  # X=[1|PCs], y


def load_blocks(sub, kind, n):
    out = []
    for b in range(ng):
        g = np.fromfile(f"{OUT}/{sub}/{kind}.{b}.bin", dtype=np.int8)
        out.append(g.reshape(n, len(g) // n))  # m = size/n; empty priv -> (n, 0)
    return out


XA, yA = load_xy("A", nA)
XB, yB = load_xy("B", nB)
A  = load_blocks("A", "geno", nA)
Bp = load_blocks("B", "geno", nB)
Bv = load_blocks("B", "priv", nB)

t0 = time.perf_counter()
skatQ, _, burden_p = federated_skat_burden_from_blocks(A, Bp, Bv, XA, yA, XB, yB, ridge_rel=RIDGE_REL)
wh, liu, davies    = federated_skat_p_from_blocks(A, Bp, Bv, XA, yA, XB, yB, ridge_rel=RIDGE_REL)
dt = time.perf_counter() - t0

print(f"nA={nA} nB={nB} genes={ng}   ({dt:.1f}s)")
print(f"{'gene':>4} {'SKAT_Q':>14} {'SKAT_p(WH)':>12} {'Liu':>12} {'Davies':>12} {'Burden_p':>12}")
for b in range(ng):
    print(f"{b:>4} {skatQ[b]:>14.4f} {wh[b]:>12.4e} {liu[b]:>12.4e} {davies[b]:>12.4e} {burden_p[b]:>12.4e}")

csv = f"{OUT}/fed_plain.csv"
with open(csv, "w") as f:
    f.write("gene,skat_Q,skat_p_wh,skat_p_liu,skat_p_davies,burden_p\n")
    for b in range(ng):
        f.write(f"{b},{skatQ[b]:.8e},{wh[b]:.8e},{liu[b]:.8e},{davies[b]:.8e},{burden_p[b]:.8e}\n")
print(f"\nwrote {csv}")
