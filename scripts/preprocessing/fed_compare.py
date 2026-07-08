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
import time

import numpy as np

from skat_plain_local import federated_Q_from_blocks, federated_Q_split_from_blocks, _fed_null

OUT = os.path.expanduser(os.environ.get("FED_OUT", "~/fed_prep_out"))
RIDGE_REL = 1e-6  # match secure computeBetaHatEnc Tikhonov ridge (gwas/skat.go ridgeRel)


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
t0 = time.perf_counter()
Qplain = federated_Q_from_blocks(
    load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
    load_blocks("B", "priv", nB, ng), XA, yA, XB, yB, ridge_rel=RIDGE_REL)
print(f"  plaintext federated_Q: {time.perf_counter() - t0:.2f}s")
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

# R^2 of secure vs plaintext (plaintext = truth): 1 - SS_res/SS_tot.
qs, qp = np.asarray(Qsec[:ng], float), np.asarray(Qplain[:ng], float)
ss_tot = float(np.sum((qp - qp.mean()) ** 2))
r2 = 1 - float(np.sum((qs - qp) ** 2)) / ss_tot if ss_tot > 0 else float("nan")
maxrel = max(abs(qs[b] - qp[b]) / max(abs(qp[b]), 1e-9) for b in range(ng))
print(f"\nR^2 (secure vs plaintext) = {r2:.6f}   |   max rel = {maxrel:.2e}")
print("MATCH (secure == plaintext)" if ok else "MISMATCH -- investigate")

# --- diagnostic: PART A vs PART B error localization (needs skat_fed_A/B.txt) ---
fA, fB = f"{OUT}/out/party2/skat_fed_A.txt", f"{OUT}/out/party2/skat_fed_B.txt"
if os.path.exists(fA) and os.path.exists(fB):
    QAp, QBp = federated_Q_split_from_blocks(
        load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
        load_blocks("B", "priv", nB, ng), XA, yA, XB, yB, ridge_rel=RIDGE_REL)
    QAs = np.atleast_1d(np.loadtxt(fA))
    QBs = np.atleast_1d(np.loadtxt(fB))
    print(f"\n--- PART A vs PART B split (which half errs) ---")
    print(f"{'gene':>4} {'A_sec':>13} {'A_plain':>13} {'A_rel':>9} | {'B_sec':>12} {'B_plain':>12} {'B_rel':>9}")
    for b in range(ng):
        arel = abs(QAs[b] - QAp[b]) / max(abs(QAp[b]), 1e-9)
        brel = abs(QBs[b] - QBp[b]) / max(abs(QBp[b]), 1e-9)
        print(f"{b:>4} {QAs[b]:>13.1f} {QAp[b]:>13.1f} {arel:>9.2e} | "
              f"{QBs[b]:>12.1f} {QBp[b]:>12.1f} {brel:>9.2e}")

# --- diagnostic: secure β̂ vs plaintext β̂ (is the null-model solve the error source?) ---
fbeta = f"{OUT}/out/party2/skat_fed_beta.txt"
if os.path.exists(fbeta):
    bsec = np.atleast_1d(np.loadtxt(fbeta))
    b_true = _fed_null(XA, yA, XB, yB, 0.0)[0]       # exact (no ridge) = SKAT truth
    b_ridge = _fed_null(XA, yA, XB, yB, RIDGE_REL)[0]  # secure design (with ridge)
    print(f"\n--- β̂ secure vs plaintext (β0 = intercept; its error × colsum(G)=O(n) amplifies) ---")
    for k in range(len(bsec)):
        rel = abs(bsec[k] - b_ridge[k]) / max(abs(b_ridge[k]), 1e-12)
        print(f"  β{k}: secure={bsec[k]:+.8f}  plain_ridge={b_ridge[k]:+.8f}  plain_true={b_true[k]:+.8f}"
              f"   |sec-ridge|/ridge={rel:.2e}")

    # Decisive test: plaintext Q built on the SECURE β̂ vs the secure Q. If they agree, the whole
    # secure Q error is explained by the wrong β̂ (everything downstream is faithful).
    Qsb = federated_Q_from_blocks(
        load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
        load_blocks("B", "priv", nB, ng), XA, yA, XB, yB, beta_override=bsec)
    print(f"\n--- plaintext-Q(secure β̂) vs secure Q  (match ⇒ β̂ is the SOLE cause) ---")
    print(f"{'gene':>4} {'Q_sec':>15} {'Q_plain(secβ̂)':>15} {'rel':>9}")
    for b in range(ng):
        rel = abs(Qsec[b] - Qsb[b]) / max(abs(Qsb[b]), 1e-9)
        print(f"{b:>4} {Qsec[b]:>15.1f} {Qsb[b]:>15.1f} {rel:>9.2e}")

    # β0-only replay: secure β̂ with just the intercept swapped for the true intercept. If Q snaps
    # back near truth, the whole error is the β0-error × colsum(G) amplification → genotype centering
    # (which zeroes colsum(G̃)) is the fix.
    b_b0fix = bsec.copy(); b_b0fix[0] = b_true[0]
    Qb0 = federated_Q_from_blocks(
        load_blocks("A", "geno", nA, ng), load_blocks("B", "geno", nB, ng),
        load_blocks("B", "priv", nB, ng), XA, yA, XB, yB, beta_override=b_b0fix)
    print(f"\n--- Q(secure β̂ but true β0) vs Q_plain (near 0 ⇒ β0 amplification is the whole error) ---")
    print(f"{'gene':>4} {'Q(secβ̂,trueβ0)':>16} {'Q_plain':>15} {'rel':>9}")
    for b in range(ng):
        rel = abs(Qb0[b] - Qplain[b]) / max(abs(Qplain[b]), 1e-9)
        print(f"{b:>4} {Qb0[b]:>16.1f} {Qplain[b]:>15.1f} {rel:>9.2e}")
