#!/usr/bin/env python3
"""Plaintext local SKAT computation for the federated-private scenario. Pure numpy.

Given full cohort genotype matrices + gene/role structure, assemble per-gene block tensors and
compute the per-gene SKAT Q two ways, checking they agree:
  - federated_skat_burden_from_blocks: PART A (public list = A + B aligned) + PART B (B private), per gene.
  - pooled_skat_burden: ground truth on the union genotype (0-filled for party-unique), per gene.
build_party_blocks does the alignment (B -> public list, public_only cols = 0) + the union.

These mirror the secure path exactly; the secure SKAT is tested against them. No plink2/AoU here --
fed_prep.py does the real genotype extraction and imports build_party_blocks from this module.
A real federated==pooled check needs real cov/pheno; call verify_blocks with them.

    python3 skat_plain_local.py    # self-check (federated == pooled on a tiny fixed example)
"""
import math

import numpy as np


def beta_density_weight(count, two_n):
    p = count / two_n
    maf = np.minimum(p, 1 - p)
    return 25.0 * (1 - maf) ** 24


def _variant_q(score, count, two_n, sigma2):
    """one variant's SKAT contribution: w²s²/(2σ̂²), w = beta-density weight."""
    w = beta_density_weight(count, two_n)
    return w * w * score * score / (2 * sigma2)


def _signed_weight(count, two_n):
    """ŵ = minor-allele-oriented weight (−w iff p̄>½), for the burden collapse z = Σ ŵ_j G_j."""
    w = beta_density_weight(count, two_n)
    return -w if count / two_n > 0.5 else w


def _variant_burden_term(score, count, two_n):
    """one variant's Burden linear term ŵ·s, ŵ = signed weight (−w iff p̄>½), oriented to the minor
    allele (R::SKAT convention). Summed per gene, then squared/scaled for the Burden statistic."""
    return _signed_weight(count, two_n) * score


def _burden_pvalue(blin, zpz, sigma2):
    """Burden p-value = erfc(√(T/2)), T/2 = Burden/zᵀPz (Burden = blin²/(2σ̂²)); χ²₁ (R::SKAT r.corr=1)."""
    if zpz <= 0:
        return 1.0
    burden = blin * blin / (2 * sigma2)
    return math.erfc(math.sqrt(burden / zpz))


def build_party_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, with_union=True):
    """Assemble per-gene blocks from full cohort genotype matrices (vectorized column gather).
      gene_keys[g]: ordered public-list keys (shared+public_only); priv_keys[g]: private keys (B only)
      roles: dict key -> 'shared'|'public_only'|'private'
      A_geno: nA x ?  ; B_geno: nB x ?  ; keycol[key] -> column index into A_geno/B_geno
    Returns A_blocks, B_aligned, B_priv, union_blocks (union_blocks empty if with_union=False)."""
    nA, nB = A_geno.shape[0], B_geno.shape[0]
    A_blocks, B_aligned, B_priv, union_blocks = [], [], [], []
    for g in range(len(gene_keys)):
        pub, priv = gene_keys[g], priv_keys[g]
        pub_cols = np.fromiter((keycol[k] for k in pub), dtype=int, count=len(pub))
        priv_cols = np.fromiter((keycol[k] for k in priv), dtype=int, count=len(priv))
        shared = np.fromiter((roles[k] == "shared" for k in pub), dtype=bool, count=len(pub))

        A_blk = A_geno[:, pub_cols].copy()                     # A has the whole public list
        B_alg = np.zeros((nB, len(pub)), dtype=A_geno.dtype)
        B_alg[:, shared] = B_geno[:, pub_cols[shared]]         # shared -> B data; public_only -> stays 0
        B_prv = B_geno[:, priv_cols].copy()                    # B private variants
        A_blocks.append(A_blk); B_aligned.append(B_alg); B_priv.append(B_prv)

        if with_union:                                         # pooled ground truth (fed_prep skips this)
            U = np.zeros((nA + nB, len(pub) + len(priv)), dtype=A_geno.dtype)
            U[:nA, :len(pub)] = A_geno[:, pub_cols]            # A: all public-list cols (shared + public_only)
            U[nA:, :len(pub)][:, shared] = B_geno[:, pub_cols[shared]]  # B: shared public-list cols
            if len(priv):
                U[nA:, len(pub):] = B_geno[:, priv_cols]       # B: private cols
            union_blocks.append(U)
    return A_blocks, B_aligned, B_priv, union_blocks


def _fed_null(XA, yA, XB, yB, ridge_rel=0.0):
    """Null model (β̂, σ², rA, rB, two_n) from the pooled normal equations. ridge_rel>0 mirrors the
    secure Tikhonov ridge (skat.go computeBetaHatEnc), intercept excluded, hub = party 1 = cohort A."""
    nA, nB = len(yA), len(yB)
    XtX = XA.T @ XA + XB.T @ XB
    Xty = XA.T @ yA + XB.T @ yB
    if ridge_rel:
        c = XA.shape[1]
        XtX_A = XA.T @ XA
        eps = ridge_rel * sum(XtX_A[k, k] for k in range(1, c)) / c
        for k in range(1, c):
            XtX[k, k] += eps
    beta = np.linalg.solve(XtX, Xty)
    sigma2 = (yA @ yA + yB @ yB - Xty @ beta) / (nA + nB - XA.shape[1])
    return beta, sigma2, yA - XA @ beta, yB - XB @ beta, 2.0 * (nA + nB), XtX


def federated_skat_burden_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB, ridge_rel=0.0):
    """Per-gene SKAT (Q), Burden, and Burden p-value from the block tensors, mirroring the secure path.
    A_blocks[g]: nA x m_pub ; B_aligned[g]: nB x m_pub (public_only cols already 0) ;
    B_priv[g]: nB x m_privg . Returns (skat, burden, burden_p), each a list of per-gene values.
    Burden = (Σ ŵ·s)²/(2σ̂²); Burden p from zᵀPz = z_fullᵀP z_full, z_full = Σ ŵ_j G_j (cohort-additive)."""
    beta, sigma2, rA, rB, two_n, XtX = _fed_null(XA, yA, XB, yB, ridge_rel)

    skat, burden, burden_p = [], [], []
    for g in range(len(A_blocks)):
        q, blin = 0.0, 0.0
        zA, zB = np.zeros(len(yA)), np.zeros(len(yB))
        # PART A: public list (A + B aligned, summed column-wise)
        for k in range(A_blocks[g].shape[1]):
            acol, bcol = A_blocks[g][:, k], B_aligned[g][:, k]
            score = acol @ rA + bcol @ rB
            count = acol.sum() + bcol.sum()
            q += _variant_q(score, count, two_n, sigma2)
            blin += _variant_burden_term(score, count, two_n)
            sw = _signed_weight(count, two_n)
            zA += sw * acol
            zB += sw * bcol
        # PART B: private (B only)
        for k in range(B_priv[g].shape[1]):
            pcol = B_priv[g][:, k]
            score = pcol @ rB
            count = pcol.sum()
            q += _variant_q(score, count, two_n, sigma2)
            blin += _variant_burden_term(score, count, two_n)
            zB += _signed_weight(count, two_n) * pcol
        zz = zA @ zA + zB @ zB
        xtz = XA.T @ zA + XB.T @ zB
        zpz = zz - xtz @ np.linalg.solve(XtX, xtz)
        skat.append(float(q))
        burden.append(float(blin * blin / (2 * sigma2)))
        burden_p.append(float(_burden_pvalue(blin, zpz, sigma2)))
    return skat, burden, burden_p


def pooled_skat_burden(union_blocks, XA, yA, XB, yB):
    """Ground truth: per-gene SKAT, Burden, Burden p on the pooled union genotype (0-filled for
    party-unique). union_blocks[g]: (nA+nB) x m_union , rows 0..nA-1 = A, nA.. = B.
    Returns (skat, burden, burden_p)."""
    nA, nB = len(yA), len(yB)
    two_n = 2.0 * (nA + nB)
    X = np.vstack([XA, XB])
    y = np.concatenate([yA, yB])
    XtX = X.T @ X
    Xty = X.T @ y
    beta = np.linalg.solve(XtX, Xty)
    sigma2 = (y @ y - Xty @ beta) / (len(y) - X.shape[1])
    r = y - X @ beta
    skat, burden, burden_p = [], [], []
    for g in range(len(union_blocks)):
        G = union_blocks[g]
        q, blin = 0.0, 0.0
        z = np.zeros(nA + nB)
        for k in range(G.shape[1]):
            score = G[:, k] @ r
            count = G[:, k].sum()
            q += _variant_q(score, count, two_n, sigma2)
            blin += _variant_burden_term(score, count, two_n)
            z += _signed_weight(count, two_n) * G[:, k]
        zpz = z @ z - (X.T @ z) @ np.linalg.solve(XtX, X.T @ z)
        skat.append(float(q))
        burden.append(float(blin * blin / (2 * sigma2)))
        burden_p.append(float(_burden_pvalue(blin, zpz, sigma2)))
    return skat, burden, burden_p


def verify_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, XA, yA, XB, yB):
    """federated(from blocks) == pooled per gene for SKAT and Burden; raises on mismatch.
    Returns (skat_fed, skat_pool)."""
    A_blocks, B_aligned, B_priv, union_blocks = build_party_blocks(
        gene_keys, priv_keys, roles, A_geno, B_geno, keycol)
    sfed, bfed, pfed = federated_skat_burden_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB)
    spool, bpool, ppool = pooled_skat_burden(union_blocks, XA, yA, XB, yB)
    for name, fed, pool in (("SKAT", sfed, spool), ("Burden", bfed, bpool), ("BurdenP", pfed, ppool)):
        for a, b in zip(fed, pool):
            assert abs(a - b) / max(abs(b), 1e-12) < 1e-9, f"federated != pooled {name} -- alignment bug"
    return sfed, spool


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
    sfed, _ = verify_blocks(gene_keys, priv_keys, roles, A_geno, B_geno, keycol, XA, yA, XB, yB)
    print(f"selfcheck OK: federated == pooled SKAT+Burden+p (gene SKAT = {sfed[0]:.4f})")


if __name__ == "__main__":
    _selfcheck()
