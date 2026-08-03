#!/usr/bin/env python3
"""Plaintext local SKAT computation for the federated-private scenario.

Given full cohort genotype matrices + gene/role structure, assemble per-gene block tensors and
compute the per-gene SKAT Q two ways, checking they agree:
  - federated_skat_burden_from_blocks: PART A (public list = A + B aligned) + PART B (B private), per gene.
  - pooled_skat_burden: ground truth on the union genotype (0-filled for party-unique), per gene.
build_party_blocks does the alignment (B -> public list, public_only cols = 0) + the union.

These mirror the secure path exactly; the secure SKAT is tested against them. No plink2/AoU here --
fed_prep.py does the real genotype extraction and imports build_party_blocks from this module.
A real federated==pooled check needs real cov/pheno; call verify_blocks with them.
The core algebra uses NumPy; the exact general-rank Davies reference is batched through R::SKAT.

    python3 skat_plain_local.py    # self-check (federated == pooled on a tiny fixed example)
"""
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import warnings

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


# --- SKAT p-value methods for Q ~ Σλχ²₁ from the kernel eigenvalues λ. WH is the secure path's
# approximation; Liu is moment matching; the exact-mixture reference uses closed forms/a bounded
# Ruben lower-tail series plus the actual batched R::SKAT Davies/qfc routine. ---
_LANCZOS = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
            12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]


def _gammln(x):
    x -= 1.0
    a, t = 0.99999999999980993, x + 7.5
    for i, gi in enumerate(_LANCZOS):
        a += gi / (x + i + 1)
    return 0.5 * math.log(2 * math.pi) + (x + 0.5) * math.log(t) - t + math.log(a)


def _gammq(a, x):  # regularized upper incomplete gamma Q(a,x) (Numerical Recipes)
    if x <= 0:
        return 1.0
    if x < a + 1.0:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(2000):
            ap += 1.0
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-15:
                break
        else:
            raise ArithmeticError("incomplete-gamma series did not converge")
        return 1.0 - s * math.exp(-x + a * math.log(x) - _gammln(a))
    b, c, d = x + 1.0 - a, 1e300, 0.0
    d, h = 1.0 / b, 1.0 / b
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        delt = d * c
        h *= delt
        if abs(delt - 1.0) < 1e-15:
            break
    else:
        raise ArithmeticError("incomplete-gamma continued fraction did not converge")
    return math.exp(-x + a * math.log(x) - _gammln(a)) * h


def _gammp(a, x):  # regularized lower incomplete gamma P(a,x)
    if x <= 0:
        return 0.0
    if x >= a + 1.0:
        return 1.0 - _gammq(a, x)
    ap, s, d = a, 1.0 / a, 1.0 / a
    for _ in range(2000):
        ap += 1.0
        d *= x / ap
        s += d
        if abs(d) < abs(s) * 1e-15:
            break
    else:
        raise ArithmeticError("incomplete-gamma series did not converge")
    return s * math.exp(-x + a * math.log(x) - _gammln(a))


def _chi2_sf(x, df):
    return _gammq(df / 2.0, x / 2.0)


def _chi2_cdf(x, df):
    return _gammp(df / 2.0, x / 2.0)


def _ncx2_sf(x, df, nc):  # non-central chi² survival = Poisson-weighted central survivals
    if nc <= 0:
        return _chi2_sf(x, df)
    lam, total, logw = nc / 2.0, 0.0, -nc / 2.0
    for j in range(2000):
        total += math.exp(logw) * _chi2_sf(x, df + 2 * j)
        if j > lam and math.exp(logw) < 1e-17 * (total + 1e-300):
            break
        logw += math.log(lam) - math.log(j + 1)
    return min(max(total, 0.0), 1.0)


def skat_wh_p(lam, Q):
    """Wilson-Hilferty screening p — same formula as the secure path (from S1,S2,S3=tr(Kᵏ))."""
    lam = np.asarray(lam, float)
    S1, S2, S3 = lam.sum(), float((lam**2).sum()), float((lam**3).sum())
    if S2 <= 0 or S3 <= 0:
        return 1.0
    u = (Q - S1) * S3 / S2**2
    h = 2 * S3**2 / (9 * S2**3)
    z = (np.cbrt(1 + u) - 1 + h) / math.sqrt(h)
    return 0.5 * math.erfc(z / math.sqrt(2))


def skat_liu_p(lam, Q):
    """Liu et al. (2009) moment-matching to a (non-central) chi-square (R::SKAT method='liu')."""
    lam = np.asarray(lam, float)
    c1, c2, c3, c4 = (float((lam**k).sum()) for k in (1, 2, 3, 4))
    if c2 <= 0:
        return 1.0
    s1, s2 = c3 / c2**1.5, c4 / c2**2
    if s1**2 > s2:
        a = 1.0 / (s1 - math.sqrt(s1**2 - s2))
        delta = s1 * a**3 - a**2
        l = a**2 - 2 * delta
    else:
        a, delta = (1.0 / s1 if s1 > 0 else 0.0), 0.0
        l = 1.0 / s2 if s2 > 0 else 0.0
    if l <= 0:
        return 1.0
    q_norm = (Q - c1) / math.sqrt(2 * c2) * (math.sqrt(2.0) * a) + (l + delta)
    return _ncx2_sf(q_norm, l, delta)


def _positive_mixture_weights(lam):
    """Remove numerical zero eigenvalues while rejecting a materially non-PSD SKAT kernel."""
    lam = np.asarray(lam, dtype=float).reshape(-1)
    if not np.all(np.isfinite(lam)):
        raise ValueError("SKAT mixture eigenvalues must be finite")
    if lam.size == 0:
        return lam
    scale = float(np.max(np.abs(lam)))
    if scale == 0.0:
        return lam[:0]
    numerical_zero = 1e-12 * scale
    if np.any(lam < -numerical_zero):
        raise ValueError("SKAT kernel has a materially negative eigenvalue")
    # Keep every positive eigenvalue. Dropping many individually-small modes can move the aggregate
    # tail by much more than the requested p-value tolerance.
    return lam[lam > 0.0]


def _normalized_mixture_case(lam, Q):
    if not math.isfinite(Q):
        raise ValueError("SKAT statistic Q must be finite")
    lam = _positive_mixture_weights(lam)
    if lam.size == 0 or Q <= 0.0:
        return lam, 0.0

    scale = float(np.max(lam))
    weights = lam / scale
    q = float(Q / scale)
    if not math.isfinite(q):
        raise OverflowError("normalized SKAT statistic Q/lambda_max is not finite")
    return weights, q


def _ruben_lower_tail_p(weights, q, tol, max_terms=4096):
    """Bounded lower-tail series for a positive central chi-square mixture.

    With beta=min(lambda), Ruben's expansion has non-negative coefficients. After term k, the
    omitted CDF is bounded by its remaining coefficient mass times the next chi-square CDF. This is
    used only below the mixture mean, where it normally converges in a handful of terms; difficult
    spectra return None and are delegated to Davies instead of returning an unchecked value.
    """
    beta = float(np.min(weights))
    z = float(q / beta)
    if beta <= 0.0 or not math.isfinite(z):
        return None

    gamma = 1.0 - beta / weights
    log_d0 = 0.5 * float(np.log(beta / weights).sum())
    if log_d0 < math.log(np.finfo(float).tiny):
        return None

    d0 = math.exp(log_d0)
    coeffs = [d0]
    mass = d0
    cdf = d0 * _chi2_cdf(z, float(weights.size))
    roundoff = 64.0 * np.finfo(float).eps

    def finish():
        pvalue = 1.0 - cdf
        if not math.isfinite(pvalue) or pvalue < -tol or pvalue > 1.0 + tol:
            return None
        return min(max(pvalue, 0.0), 1.0)

    next_cdf = _chi2_cdf(z, float(weights.size + 2))
    if max(0.0, 1.0 - mass) * next_cdf + roundoff <= tol:
        return finish()

    gamma_power = np.ones_like(gamma)
    b_terms = []
    for k in range(1, max_terms + 1):
        gamma_power *= gamma
        b_terms.append(0.5 * float(gamma_power.sum()))
        dk = float(np.dot(b_terms, coeffs[::-1])) / k
        if not math.isfinite(dk) or dk < 0.0:
            return None
        coeffs.append(dk)
        mass += dk
        cdf += dk * _chi2_cdf(z, float(weights.size + 2 * k))
        if mass > 1.0 + 1e-8 or cdf > 1.0 + 1e-8:
            return None
        next_cdf = _chi2_cdf(z, float(weights.size + 2 * (k + 1)))
        error_bound = max(0.0, 1.0 - mass) * next_cdf + roundoff * (k + 1)
        if error_bound <= tol:
            return finish()
    return None


def _run_r_davies_batch(cases, tol, on_unresolved):
    """Evaluate normalized general-rank cases in one R::SKAT qfc/Davies process."""
    if not cases:
        return []
    rscript = os.environ.get("SKAT_RSCRIPT") or shutil.which("Rscript")
    helper = Path(__file__).with_name("skat_davies_batch.R")
    if not rscript:
        if on_unresolved == "nan":
            return [float("nan")] * len(cases)
        raise RuntimeError(
            "Rscript with the R package SKAT is required for general-rank Davies p-values"
        )
    if not helper.is_file():
        raise RuntimeError(f"missing Davies helper: {helper}")

    with tempfile.TemporaryDirectory(prefix="secure-skat-davies-") as tmp:
        input_path = Path(tmp) / "cases.csv"
        output_path = Path(tmp) / "results.csv"
        with input_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("case_id", "eigen_index", "q", "lambda"))
            for case_id, (weights, q) in enumerate(cases):
                for eigen_index, value in enumerate(weights):
                    writer.writerow((case_id, eigen_index, format(q, ".17g"),
                                     format(float(value), ".17g")))

        try:
            run = subprocess.run(
                [rscript, str(helper), str(input_path), str(output_path), format(tol, ".17g")],
                text=True, capture_output=True, check=False, timeout=300)
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise RuntimeError(f"could not run the R::SKAT Davies batch: {exc}") from exc
        if run.returncode != 0:
            detail = (run.stderr or run.stdout).strip().splitlines()
            suffix = f": {' | '.join(detail[-4:])}" if detail else ""
            raise RuntimeError(f"R::SKAT Davies batch failed{suffix}")

        with output_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))

    if len(rows) != len(cases):
        raise RuntimeError(
            f"R::SKAT Davies returned {len(rows)} rows for {len(cases)} cases")
    result = []
    for expected_id, row in enumerate(rows):
        case_id = int(row["case_id"])
        pvalue = float(row["pvalue"])
        ifault = int(row["ifault"])
        acc_used = float(row["acc_used"])
        if case_id != expected_id:
            raise RuntimeError("R::SKAT Davies returned cases out of order")
        if ifault != 0 or not math.isfinite(pvalue) or pvalue <= 0.0 or pvalue > 1.0:
            detail = f"case {case_id} (p={pvalue}, ifault={ifault}, acc={acc_used:g})"
            if on_unresolved == "raise":
                raise RuntimeError(f"R::SKAT Davies did not resolve {detail}")
            result.append(float("nan"))
        else:
            result.append(pvalue)
    return result


def skat_davies_p_batch(cases, tol=1e-10, on_unresolved="raise"):
    """Exact upper tails for multiple Q~sum(lambda_i*chi-square_1) cases.

    Rank-one/equal spectra use their closed form, convergent lower tails use a bounded Ruben series,
    and all remaining cases are evaluated together by the actual R::SKAT Davies/qfc routine. Invalid
    qfc status and unresolved extreme tails are never clipped into 0, 0.5, or 1. By default they
    raise; on_unresolved="nan" preserves the other batch results and emits a warning.
    """
    if not (1e-10 <= tol < 1.0):
        raise ValueError("Davies tolerance must lie in [1e-10,1)")
    if on_unresolved not in ("raise", "nan"):
        raise ValueError("on_unresolved must be 'raise' or 'nan'")
    cases = list(cases)
    result = [None] * len(cases)
    pending, pending_indices = [], []
    for index, (lam, Q) in enumerate(cases):
        weights, q = _normalized_mixture_case(lam, Q)
        if weights.size == 0 or Q <= 0.0:
            result[index] = 1.0
            continue

        # These are exact and avoid asking a general-purpose routine to resolve a special case. If
        # their survival underflows (or the gamma iteration fails), delegate to qfc so the configured
        # unresolved policy is applied instead of silently storing p=0.
        if weights.size == 1:
            try:
                special = _chi2_sf(q / float(weights[0]), 1.0)
            except ArithmeticError:
                special = 0.0
            if math.isfinite(special) and special > 0.0:
                result[index] = special
                continue
        elif float(np.max(np.abs(weights - weights[0]))) <= 16.0 * np.finfo(float).eps:
            try:
                special = _chi2_sf(q / float(weights[0]), float(weights.size))
            except ArithmeticError:
                special = 0.0
            if math.isfinite(special) and special > 0.0:
                result[index] = special
                continue

        try:
            ruben = _ruben_lower_tail_p(weights, q, tol) if q <= float(weights.sum()) else None
        except ArithmeticError:
            ruben = None
        if ruben is not None:
            result[index] = ruben
        else:
            pending_indices.append(index)
            pending.append((weights, q))

    unresolved_indices = []
    for index, pvalue in zip(
            pending_indices, _run_r_davies_batch(pending, tol, on_unresolved)):
        result[index] = pvalue
        if math.isnan(pvalue):
            unresolved_indices.append(index)
    if unresolved_indices:
        preview = ", ".join(str(index) for index in unresolved_indices[:10])
        suffix = ", ..." if len(unresolved_indices) > 10 else ""
        warnings.warn(
            f"R::SKAT Davies left {len(unresolved_indices)}/{len(cases)} cases unresolved "
            f"(case indices: {preview}{suffix}); their reference p-values are NaN",
            RuntimeWarning,
            stacklevel=2,
        )
    return [float(value) for value in result]


def skat_davies_p(lam, Q, tol=1e-10):
    """Single-case wrapper around skat_davies_p_batch."""
    return skat_davies_p_batch([(lam, Q)], tol=tol)[0]


def federated_skat_p_from_blocks(A_blocks, B_aligned, B_priv, XA, yA, XB, yB, ridge_rel=0.0):
    """Per-gene SKAT p-value from the block tensors, three ways: Wilson-Hilferty (the secure path's
    approximation), Liu moment-matching, and exact Davies. Kernel K=½D(GᵀPG)D, D=diag(w).
    Returns (wh, liu, davies)."""
    beta, sigma2, rA, rB, two_n, XtX = _fed_null(XA, yA, XB, yB, ridge_rel)
    nA, nB = len(yA), len(yB)
    X = np.vstack([XA, XB])
    r = np.concatenate([rA, rB])
    wh, liu, davies_cases = [], [], []
    for g in range(len(A_blocks)):
        mpub, mpriv = A_blocks[g].shape[1], B_priv[g].shape[1]
        G = np.zeros((nA + nB, mpub + mpriv))          # union genotype: A private cols stay 0
        G[:nA, :mpub] = A_blocks[g]
        G[nA:, :mpub] = B_aligned[g]
        if mpriv:
            G[nA:, mpub:] = B_priv[g]
        cnt = G.sum(0)
        w = 25.0 * (1 - np.minimum(cnt / two_n, 1 - cnt / two_n)) ** 24
        s = G.T @ r
        Q = float(np.sum(w**2 * s**2) / (2 * sigma2))   # SKAT statistic
        M = G.T @ G - (G.T @ X) @ np.linalg.solve(X.T @ X, (X.T @ G))
        lam = np.linalg.eigvalsh(0.5 * (w[:, None] * M) * w[None, :])   # eig of K = ½ D M D
        wh.append(skat_wh_p(lam, Q))
        liu.append(skat_liu_p(lam, Q))
        davies_cases.append((lam, Q))
    # A single unresolved extreme tail must not discard every other gene. Preserve it as NaN; summary
    # metrics and plots exclude non-finite references and the batch emits an explicit warning.
    return wh, liu, skat_davies_p_batch(davies_cases, on_unresolved="nan")


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
