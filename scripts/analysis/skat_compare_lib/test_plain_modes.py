"""Experimental plain-mode variants kept separate from the main SKAT math path."""

from __future__ import annotations

import numpy as np

from .compute import compute_beta_weight

PLAIN_MODE_LOCAL_WEIGHT_BURDEN = "local-weight-burden"
LOCAL_WEIGHT_MODE_DIRECT_TOTAL = "direct-total"
LOCAL_WEIGHT_MODE_PRODUCT_APPROX = "product-approx"


def compute_local_weight_burden_direct_total_block(
    arg_G_parts: list[np.ndarray],
    arg_y_resid: np.ndarray,
    arg_n_total: int,
) -> dict:
    G_parts = [np.asarray(G_part, dtype=float) for G_part in arg_G_parts]
    y_resid = np.asarray(arg_y_resid, dtype=float).reshape(-1)
    n_total = int(arg_n_total)
    if n_total <= 0:
        raise RuntimeError("Total sample count must be positive for local-weight burden mode")

    local_burden_terms = []
    local_details = []
    offset = 0

    for party_index, G_part in enumerate(G_parts, start=1):
        n_local = int(G_part.shape[0])
        y_local = y_resid[offset : offset + n_local]
        if y_local.size != n_local:
            raise RuntimeError("Residual vector length mismatch while slicing party-local block rows")
        offset += n_local

        secure_dosage_sum_local = G_part.sum(axis=0)
        p_bar_local = secure_dosage_sum_local / (2.0 * n_total)
        p_local = 1.0 - p_bar_local
        weight_local = compute_beta_weight(np.maximum(p_local, p_bar_local))
        score_local = G_part.T @ y_local
        burden_local = float(np.sum(weight_local * score_local))

        local_burden_terms.append(burden_local)
        local_details.append(
            {
                "party_index": party_index,
                "n_local": n_local,
                "n_total": n_total,
                "p_local": p_local,
                "p_bar_local": p_bar_local,
                "weight_local": weight_local,
                "score_local": score_local,
                "burden_local_raw": burden_local,
            }
        )

    if offset != y_resid.size:
        raise RuntimeError("Residual vector contains unexpected extra rows after local burden slicing")

    return {
        "q_burden_block_raw": float(sum(local_burden_terms)),
        "local_weight_details": local_details,
    }


def compute_local_weight_burden_product_approx_block(
    arg_G_parts: list[np.ndarray],
    arg_y_resid: np.ndarray,
    arg_n_total: int,
) -> dict:
    G_parts = [np.asarray(G_part, dtype=float) for G_part in arg_G_parts]
    y_resid = np.asarray(arg_y_resid, dtype=float).reshape(-1)
    n_total = int(arg_n_total)
    if n_total <= 0:
        raise RuntimeError("Total sample count must be positive for local-weight burden mode")

    local_factor_vec = None
    local_details = []
    offset = 0
    score_sum = None

    for party_index, G_part in enumerate(G_parts, start=1):
        n_local = int(G_part.shape[0])
        y_local = y_resid[offset : offset + n_local]
        if y_local.size != n_local:
            raise RuntimeError("Residual vector length mismatch while slicing party-local block rows")
        offset += n_local

        secure_dosage_sum_local = G_part.sum(axis=0)
        p_bar_local = secure_dosage_sum_local / (2.0 * n_total)
        factor_local = np.power(p_bar_local, 24)
        score_local = G_part.T @ y_local

        if local_factor_vec is None:
            local_factor_vec = factor_local.copy()
            score_sum = score_local.copy()
        else:
            local_factor_vec *= factor_local
            score_sum += score_local

        local_details.append(
            {
                "party_index": party_index,
                "n_local": n_local,
                "n_total": n_total,
                "p_bar_local": p_bar_local,
                "weight_factor_local": factor_local,
                "score_local": score_local,
            }
        )

    if offset != y_resid.size:
        raise RuntimeError("Residual vector contains unexpected extra rows after local burden slicing")
    if local_factor_vec is None or score_sum is None:
        raise RuntimeError("No party-local genotype blocks were provided for product-approx burden mode")

    weight_approx = 25.0 * local_factor_vec
    return {
        "weight_vec": weight_approx,
        "score_vec": score_sum,
        "q_burden_block_raw": float(np.sum(weight_approx * score_sum)),
        "local_weight_details": local_details,
    }


def apply_plain_test_mode(
    arg_block_result: dict,
    arg_G_parts: list[np.ndarray],
    arg_y_resid: np.ndarray,
    arg_n_total: int,
    arg_plain_mode: str,
    arg_local_weight_mode: str,
) -> dict:
    if arg_plain_mode != PLAIN_MODE_LOCAL_WEIGHT_BURDEN:
        raise RuntimeError(f"Unsupported plain mode in test helper: {arg_plain_mode}")

    if arg_local_weight_mode == LOCAL_WEIGHT_MODE_DIRECT_TOTAL:
        local_weight_result = compute_local_weight_burden_direct_total_block(arg_G_parts, arg_y_resid, arg_n_total)
    elif arg_local_weight_mode == LOCAL_WEIGHT_MODE_PRODUCT_APPROX:
        local_weight_result = compute_local_weight_burden_product_approx_block(arg_G_parts, arg_y_resid, arg_n_total)
    else:
        raise RuntimeError(f"Unsupported local weight mode: {arg_local_weight_mode}")

    out = dict(arg_block_result)
    out["q_burden_block_raw"] = local_weight_result["q_burden_block_raw"]
    out["burden_mode"] = f"{PLAIN_MODE_LOCAL_WEIGHT_BURDEN}/{arg_local_weight_mode}"
    out["local_weight_details"] = local_weight_result["local_weight_details"]
    if "weight_vec" in local_weight_result:
        out["weight_vec"] = local_weight_result["weight_vec"]
    if "score_vec" in local_weight_result:
        out["score_vec"] = local_weight_result["score_vec"]
    return out
