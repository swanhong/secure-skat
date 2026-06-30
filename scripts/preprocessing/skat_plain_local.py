#!/usr/bin/env python3
"""Plaintext local SKAT computation for the federated-private scenario. Pure numpy.

Given full cohort genotype matrices + gene/role structure, assemble per-gene block tensors and
compute the per-gene SKAT Q two ways, checking they agree:
  - federated_Q_from_blocks: PART A (public list = A + B aligned) + PART B (B private), per gene.
  - pooled_Q: ground truth on the union genotype (0-filled for party-unique), per gene.
build_party_blocks does the alignment (B -> public list, public_only cols = 0) + the union.

These mirror the secure path exactly; the secure SKAT is tested against them. No plink2/AoU here --
fed_prep.py does the real genotype extraction and imports build_party_blocks from this module.
A real federated==pooled check needs real cov/pheno; call verify_blocks with them.

    python3 skat_plain_local.py    # self-check (federated == pooled on a tiny fixed example)
"""
import numpy as np


def beta_density_weight(count, two_n):
    p = count / two_n
    maf = np.minimum(p, 1 - p)
    return 25.0 * (1 - maf) ** 24


def _variant_q(score, count, two_n, sigma2):
    """one variant's SKAT contribution: w²s²/(2σ̂²), w = beta-density weight."""
    w = beta_density_weight(count, two_n)
    return w * w * score * score / (2 * sigma2)


def build_party_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol):
    """Assemble per-gene blocks from full cohort genotype matrices.
      gene_keys[g]: ordered public-list keys (shared+public_only); priv_keys[g]: private keys (B only)
      roles: dict key -> 'shared'|'public_only'|'private'
      A_geno: nA x ?  ; B_geno: nB x ?  ; keycol[key] -> column index into A_geno/B_geno
    Returns A_blocks, B_aligned, B_priv, union_blocks."""
    nA, nB = A_geno.shape[0], B_geno.shape[0]
    A_blocks, B_aligned, B_priv, union_blocks = [], [], [], []
    for g in range(len(gene_keys)):
        pub = gene_keys[g]
        priv = priv_keys[g]
        # PART A blocks, columns in public-list order
        A_blk = np.zeros((nA, len(pub)), dtype=np.float64)
        B_alg = np.zeros((nB, len(pub)), dtype=np.float64)
        for k, key in enumerate(pub):
            A_blk[:, k] = A_geno[:, keycol[key]]               # A has the whole public list
            if roles[key] == "shared":
                B_alg[:, k] = B_geno[:, keycol[key]]           # shared -> B data; public_only -> stays 0
        # PART B block (B private)
        B_prv = np.zeros((nB, len(priv)), dtype=np.float64)
        for k, key in enumerate(priv):
            B_prv[:, k] = B_geno[:, keycol[key]]
        # pooled union (shared both; public_only A only; private B only)
        union = pub + priv
        U = np.zeros((nA + nB, len(union)), dtype=np.float64)
        for k, key in enumerate(union):
            if roles[key] in ("shared", "public_only"):
                U[:nA, k] = A_geno[:, keycol[key]]
            if roles[key] in ("shared", "private"):
                U[nA:, k] = B_geno[:, keycol[key]]
        A_blocks.append(A_blk); B_aligned.append(B_alg); B_priv.append(B_prv)
        union_blocks.append(U)
    return A_blocks, B_aligned, B_priv, union_blocks


def federated_Q_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB):
    """Per-gene Q from the block tensors, mirroring the secure path exactly.
    A_blocks[g]: nA x m_pub ; B_aligned[g]: nB x m_pub (public_only cols already 0) ;
    B_priv[g]: nB x m_privg . Returns list of per-gene Q."""
    nA, nB = len(yA), len(yB)
    two_n = 2.0 * (nA + nB)
    XtX = XA.T @ XA + XB.T @ XB
    Xty = XA.T @ yA + XB.T @ yB
    beta = np.linalg.solve(XtX, Xty)
    sigma2 = (yA @ yA + yB @ yB - Xty @ beta) / (nA + nB - XA.shape[1])
    rA, rB = yA - XA @ beta, yB - XB @ beta

    Q = []
    for g in range(len(A_blocks)):
        q = 0.0
        # PART A: public list (A + B aligned, summed column-wise)
        for k in range(A_blocks[g].shape[1]):
            score = A_blocks[g][:, k] @ rA + B_aligned[g][:, k] @ rB
            count = A_blocks[g][:, k].sum() + B_aligned[g][:, k].sum()
            q += _variant_q(score, count, two_n, sigma2)
        # PART B: private (B only)
        for k in range(B_priv[g].shape[1]):
            score = B_priv[g][:, k] @ rB
            count = B_priv[g][:, k].sum()
            q += _variant_q(score, count, two_n, sigma2)
        Q.append(float(q))
    return Q


def pooled_Q(union_blocks, XA, yA, XB, yB):
    """Ground truth: per-gene SKAT on the pooled union genotype (0-filled for party-unique).
    union_blocks[g]: (nA+nB) x m_union , rows 0..nA-1 = A, nA.. = B."""
    nA, nB = len(yA), len(yB)
    two_n = 2.0 * (nA + nB)
    X = np.vstack([XA, XB])
    y = np.concatenate([yA, yB])
    XtX = X.T @ X
    Xty = X.T @ y
    beta = np.linalg.solve(XtX, Xty)
    sigma2 = (y @ y - Xty @ beta) / (len(y) - X.shape[1])
    r = y - X @ beta
    Q = []
    for g in range(len(union_blocks)):
        G = union_blocks[g]
        q = 0.0
        for k in range(G.shape[1]):
            score = G[:, k] @ r
            count = G[:, k].sum()
            q += _variant_q(score, count, two_n, sigma2)
        Q.append(float(q))
    return Q


def verify_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, XA, yA, XB, yB):
    """federated(from blocks) == pooled per gene; raises on mismatch. Returns (qfed, qpool)."""
    A_blocks, B_aligned, B_priv, union_blocks = build_party_blocks(
        gene_keys, priv_keys, roles, A_geno, B_geno, keycol)
    qfed = federated_Q_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB)
    qpool = pooled_Q(union_blocks, XA, yA, XB, yB)
    for a, b in zip(qfed, qpool):
        assert abs(a - b) / max(abs(b), 1e-12) < 1e-9, "federated != pooled -- alignment bug"
    return qfed, qpool


def _selfcheck():
    """1 gene, 3 variants (shared / public_only / private); checks federated == pooled."""
    roles = {"v_sh": "shared", "v_pub": "public_only", "v_prv": "private"}
    gene_keys = [["v_sh", "v_pub"]]
    priv_keys = [["v_prv"]]
    keycol = {"v_sh": 0, "v_pub": 1, "v_prv": 2}
    A_geno = np.array([[1, 0, 0], [2, 1, 0], [0, 0, 0], [1, 2, 0]], float)  # private col unused (A lacks it)
    B_geno = np.array([[0, 0, 1], [1, 0, 2], [2, 0, 0]], float)             # public_only col unused (B lacks it)
    XA = np.array([[1, 0.2], [1, -0.5], [1, 1.0], [1, 0.3]]); yA = np.array([0.5, -0.2, 1.1, 0.0])
    XB = np.array([[1, -0.1], [1, 0.7], [1, 0.4]]); yB = np.array([0.3, -0.6, 0.8])
    qfed, _ = verify_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, XA, yA, XB, yB)
    print(f"selfcheck OK: federated == pooled (gene Q = {qfed[0]:.4f})")


if __name__ == "__main__":
    _selfcheck()
