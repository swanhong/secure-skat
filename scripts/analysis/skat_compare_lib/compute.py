"""Core SKAT/Burden math with secure-compatible semantics."""

from __future__ import annotations

import numpy as np


def fit_null_model(arg_X: np.ndarray, arg_y: np.ndarray) -> dict:
    X = np.asarray(arg_X, dtype=float)
    y = np.asarray(arg_y, dtype=float).reshape(-1)
    n_total = int(y.size)

    # Build Z = [1, X] in R^(n x k), where k includes the intercept.
    Z = np.column_stack([np.ones(n_total, dtype=float), X])
    k_total = int(Z.shape[1])

    # Solve beta = argmin ||y - Z beta||_2^2.
    beta_hat, *_ = np.linalg.lstsq(Z, y, rcond=None)

    # Compute y_resid = y - Z beta, dim(y_resid) = n.
    y_resid = y - (Z @ beta_hat)
    if np.isnan(y_resid).any():
        raise RuntimeError("Null-model residual calculation produced NaN values")

    # Compute RSS = ||y_resid||_2^2 and rare_variant_scale = (n - k) / (2 * RSS).
    null_rss = float(np.sum(y_resid * y_resid))
    dof = n_total - k_total
    if dof <= 0:
        raise RuntimeError("Null-model degrees of freedom must be positive")
    if null_rss <= 0.0:
        raise RuntimeError("Null-model residual sum of squares must be positive")

    return {
        "X": X,
        "y": y,
        "y_resid": y_resid,
        "n_total": n_total,
        "k_total": k_total,
        "dof": dof,
        "null_rss": null_rss,
        "rare_variant_scale": float(dof / (2.0 * null_rss)),
    }


def compute_beta_weight(arg_beta_base_vec: np.ndarray) -> np.ndarray:
    return 25.0 * np.power(np.asarray(arg_beta_base_vec, dtype=float), 24)


def compute_manual_block(arg_G_parts: list[np.ndarray], arg_y_resid: np.ndarray, arg_n_total: int) -> dict:
    G_parts = [np.asarray(G_part, dtype=float) for G_part in arg_G_parts]
    y_resid = np.asarray(arg_y_resid, dtype=float).reshape(-1)
    n_total = int(arg_n_total)

    # Stack party matrices into G in R^(n x m).
    G = np.vstack(G_parts)
    if G.shape[0] != n_total:
        raise RuntimeError("Sample count mismatch while stacking block genotype matrix")

    # Compute alternative-allele counts d in R^m and secure-oriented counts d_bar in R^m.
    alt_dosage_sum_vec = G.sum(axis=0)
    secure_dosage_sum_vec = (2.0 * n_total) - alt_dosage_sum_vec

    # Compute p_bar = d_bar / (2n), p = 1 - p_bar, and w in R^m.
    p_bar_vec = secure_dosage_sum_vec / (2.0 * n_total)
    p_vec = 1.0 - p_bar_vec
    beta_base_vec = np.maximum(p_vec, p_bar_vec)
    weight_vec = compute_beta_weight(beta_base_vec)

    # Compute s = -G^T y_resid, dim(s) = m.
    score_vec = -(G.T @ y_resid)

    # Aggregate q_skat_block = sum((w^2) * (s^2)) and q_burden_block = sum(w * s).
    q_skat_block_raw = float(np.sum((weight_vec * weight_vec) * (score_vec * score_vec)))
    q_burden_block_raw = float(np.sum(weight_vec * score_vec))

    return {
        "weight_vec": weight_vec,
        "q_skat_block_raw": q_skat_block_raw,
        "q_burden_block_raw": q_burden_block_raw,
    }


def compute_manual_results(arg_ctx: dict, arg_block_inputs: list[dict]) -> dict:
    y_resid = arg_ctx["model"]["y_resid"]
    n_total = arg_ctx["model"]["n_total"]
    rare_variant_scale = arg_ctx["model"]["rare_variant_scale"]

    blocks = []
    all_q_skat_raw_total = 0.0
    all_q_burden_raw_total = 0.0

    for block_input in arg_block_inputs:
        block_math = compute_manual_block(
            block_input["local_alt_genotypes_by_party"],
            y_resid,
            n_total,
        )

        block_result = {
            "block_index": block_input["block_index"],
            "n_variants": len(block_input["variant_ids"]),
            "variant_ids": block_input["variant_ids"],
            "positions": block_input["positions"],
            "raw_paths_by_party": block_input["raw_paths_by_party"],
            "weight_vec": block_math["weight_vec"],
            "q_skat_block_raw": block_math["q_skat_block_raw"],
            "q_burden_block_raw": block_math["q_burden_block_raw"],
        }
        blocks.append(block_result)
        all_q_skat_raw_total += block_result["q_skat_block_raw"]
        all_q_burden_raw_total += block_result["q_burden_block_raw"]

    selected_blocks = [blocks[block_index - 1] for block_index in arg_ctx["analysis_blocks"]]
    analysis_q_skat_raw_total = float(sum(block["q_skat_block_raw"] for block in selected_blocks))
    analysis_q_burden_raw_total = float(sum(block["q_burden_block_raw"] for block in selected_blocks))

    return {
        "blocks": blocks,
        "analysis_q_skat_raw_total": analysis_q_skat_raw_total,
        "analysis_q_burden_raw_total": analysis_q_burden_raw_total,
        "analysis_skat_q": analysis_q_skat_raw_total * rare_variant_scale,
        "analysis_burden_q": (analysis_q_burden_raw_total**2) * rare_variant_scale,
        "all_q_skat_raw_total": all_q_skat_raw_total,
        "all_q_burden_raw_total": all_q_burden_raw_total,
        "all_skat_q": all_q_skat_raw_total * rare_variant_scale,
        "all_burden_q": (all_q_burden_raw_total**2) * rare_variant_scale,
    }
