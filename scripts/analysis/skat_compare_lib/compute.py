"""Secure-compatible SKAT/Burden math and block/window aggregation."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .common import safe_rel_diff
from .models import BlockComparisonRow, CompareContext, DatasetInputs, ManualBlockResult, ManualResults, WindowComparisonRow


def compute_beta_weight(beta_base_vec: np.ndarray) -> np.ndarray:
    return 25.0 * np.power(beta_base_vec, 24)


def compute_block_statistics(weight_vec: np.ndarray, score_vec: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    score_sq_vec = score_vec**2
    weight_sq_vec = weight_vec**2
    w2s2_vec = weight_sq_vec * score_sq_vec
    wS_vec = weight_vec * score_vec
    q_skat_block_raw = float(np.sum(w2s2_vec))
    q_burden_block_raw = float(np.sum(wS_vec))
    return score_sq_vec, weight_sq_vec, w2s2_vec, wS_vec, q_skat_block_raw, q_burden_block_raw


def compute_manual_results(ctx: CompareContext, dataset_inputs: DatasetInputs) -> ManualResults:
    blocks: list[ManualBlockResult] = []
    variant_total = 0
    all_q_skat_raw_total = 0.0
    all_q_burden_raw_total = 0.0

    for block_input in dataset_inputs.block_inputs:
        local_alt_genotypes_by_party = [np.asarray(geno, dtype=float) for geno in block_input.local_alt_genotypes_by_party]
        stacked_alt_geno = np.vstack(local_alt_genotypes_by_party)
        if stacked_alt_geno.shape[0] != ctx.model.n_total:
            raise RuntimeError(f"Sample count mismatch at block {block_input.block_index}")

        local_alt_dosage_sum_vecs = [geno.sum(axis=0) for geno in local_alt_genotypes_by_party]
        local_secure_dosage_sum_vecs = [
            (2.0 * geno.shape[0]) - local_alt_dosage_sum_vec
            for geno, local_alt_dosage_sum_vec in zip(local_alt_genotypes_by_party, local_alt_dosage_sum_vecs)
        ]

        alt_dosage_sum_vec = stacked_alt_geno.sum(axis=0)
        secure_dosage_sum_vec = (2.0 * ctx.model.n_total) - alt_dosage_sum_vec
        p_bar_vec = secure_dosage_sum_vec / (2.0 * ctx.model.n_total)
        p_vec = 1.0 - p_bar_vec
        beta_base_vec = np.maximum(p_vec, p_bar_vec)
        weight_vec = compute_beta_weight(beta_base_vec)

        score_alt_vec = stacked_alt_geno.T @ ctx.model.y_resid
        score_vec = -score_alt_vec
        score_sq_vec, weight_sq_vec, w2s2_vec, wS_vec, q_skat_block_raw, q_burden_block_raw = compute_block_statistics(
            weight_vec,
            score_vec,
        )

        n_variants = int(stacked_alt_geno.shape[1])
        variant_total += n_variants
        all_q_skat_raw_total += q_skat_block_raw
        all_q_burden_raw_total += q_burden_block_raw

        blocks.append(
            ManualBlockResult(
                block_index=block_input.block_index,
                raw_paths_by_party=block_input.raw_paths_by_party,
                n_variants=n_variants,
                variant_ids=block_input.variant_ids,
                positions=block_input.positions,
                local_alt_dosage_sum_vecs=local_alt_dosage_sum_vecs,
                local_secure_dosage_sum_vecs=local_secure_dosage_sum_vecs,
                alt_dosage_sum_vec=alt_dosage_sum_vec,
                secure_dosage_sum_vec=secure_dosage_sum_vec,
                p_vec=p_vec,
                p_bar_vec=p_bar_vec,
                beta_base_vec=beta_base_vec,
                weight_vec=weight_vec,
                score_vec=score_vec,
                score_sq_vec=score_sq_vec,
                weight_sq_vec=weight_sq_vec,
                w2s2_vec=w2s2_vec,
                wS_vec=wS_vec,
                q_skat_block_raw=q_skat_block_raw,
                q_burden_block_raw=q_burden_block_raw,
            )
        )

    selected_blocks = [blocks[block - 1] for block in ctx.analysis_blocks]
    analysis_variant_total = sum(block.n_variants for block in selected_blocks)
    analysis_q_skat_raw_total = float(sum(block.q_skat_block_raw for block in selected_blocks))
    analysis_q_burden_raw_total = float(sum(block.q_burden_block_raw for block in selected_blocks))

    return ManualResults(
        blocks=blocks,
        variant_total=variant_total,
        analysis_variant_total=analysis_variant_total,
        analysis_q_skat_raw_total=analysis_q_skat_raw_total,
        analysis_skat_q=analysis_q_skat_raw_total * ctx.model.rare_variant_scale,
        analysis_q_burden_raw_total=analysis_q_burden_raw_total,
        analysis_burden_q=(analysis_q_burden_raw_total**2) * ctx.model.rare_variant_scale,
        all_q_skat_raw_total=all_q_skat_raw_total,
        all_skat_q=all_q_skat_raw_total * ctx.model.rare_variant_scale,
        all_q_burden_raw_total=all_q_burden_raw_total,
        all_burden_q=(all_q_burden_raw_total**2) * ctx.model.rare_variant_scale,
    )


def make_bp_windows(pos_vec: np.ndarray, block_index: int, window_bp: int, step_bp: int, min_variants: int) -> pd.DataFrame:
    if pos_vec.size == 0:
        return pd.DataFrame()
    starts = np.arange(int(pos_vec[0]), int(pos_vec[-1]) + 1, step_bp, dtype=int)
    left_idx = np.searchsorted(pos_vec, starts, side="left")
    right_idx = np.searchsorted(pos_vec, starts + window_bp - 1, side="right") - 1
    n_variants = np.maximum(0, right_idx - left_idx + 1)
    keep = (left_idx < pos_vec.size) & (right_idx >= left_idx) & (n_variants >= min_variants)
    if not np.any(keep):
        return pd.DataFrame()
    return pd.DataFrame(
        {
            "block": np.full(int(np.sum(keep)), block_index, dtype=int),
            "block_window_index": np.nonzero(keep)[0] + 1,
            "window_start_bp": starts[keep],
            "window_end_bp": starts[keep] + window_bp - 1,
            "start_index": left_idx[keep],
            "end_index": right_idx[keep],
            "n_variants": n_variants[keep],
        }
    )


def build_block_comparison_rows(ctx: CompareContext, manual: ManualResults) -> list[BlockComparisonRow]:
    rows: list[BlockComparisonRow] = []

    for block_index in ctx.analysis_blocks:
        block = manual.blocks[block_index - 1]
        secure = block.secure_artifacts

        plain_skat_q = block.q_skat_block_raw * ctx.model.rare_variant_scale
        secure_skat_q = float("nan")
        if secure is not None and secure.q_skat_block_raw is not None:
            secure_skat_q = secure.q_skat_block_raw * ctx.model.rare_variant_scale

        plain_burden_q = (block.q_burden_block_raw**2) * ctx.model.rare_variant_scale
        secure_burden_q = float("nan")
        secure_burden_sum = float("nan")
        if secure is not None and secure.q_burden_block_raw is not None:
            secure_burden_sum = secure.q_burden_block_raw
            secure_burden_q = (secure.q_burden_block_raw**2) * ctx.model.rare_variant_scale

        rows.append(
            BlockComparisonRow(
                block=block.block_index,
                n_variants=block.n_variants,
                block_start_bp=int(block.positions[0]) if block.n_variants else math.nan,
                block_end_bp=int(block.positions[-1]) if block.n_variants else math.nan,
                start_variant_id=block.variant_ids[0] if block.n_variants else "",
                end_variant_id=block.variant_ids[-1] if block.n_variants else "",
                plain_skat_q=plain_skat_q,
                secure_skat_q=secure_skat_q,
                skat_abs_diff=abs(plain_skat_q - secure_skat_q) if np.isfinite(secure_skat_q) else math.nan,
                skat_rel_diff=safe_rel_diff(plain_skat_q, secure_skat_q),
                plain_burden_q=plain_burden_q,
                secure_burden_q=secure_burden_q,
                burden_abs_diff=abs(plain_burden_q - secure_burden_q) if np.isfinite(secure_burden_q) else math.nan,
                burden_rel_diff=safe_rel_diff(plain_burden_q, secure_burden_q),
                plain_burden_sum=block.q_burden_block_raw,
                secure_burden_sum=secure_burden_sum,
            )
        )

    return rows


def build_window_comparison_rows(ctx: CompareContext, manual: ManualResults) -> list[WindowComparisonRow]:
    if ctx.window_bp is None:
        return []

    window_rows: list[WindowComparisonRow] = []
    row_count = 0
    for block_index in ctx.analysis_blocks:
        block = manual.blocks[block_index - 1]
        secure = block.secure_artifacts
        block_windows = make_bp_windows(
            block.positions,
            block.block_index,
            ctx.window_bp,
            ctx.step_bp or ctx.window_bp,
            ctx.min_window_variants,
        )
        if block_windows.empty:
            continue

        if ctx.window_limit is not None:
            remaining = ctx.window_limit - row_count
            if remaining <= 0:
                break
            block_windows = block_windows.iloc[:remaining].copy()

        for _, window in block_windows.iterrows():
            start_idx = int(window["start_index"])
            end_idx = int(window["end_index"])
            plain_skat_q = float(np.sum(block.w2s2_vec[start_idx : end_idx + 1]) * ctx.model.rare_variant_scale)
            plain_burden_sum = float(np.sum(block.wS_vec[start_idx : end_idx + 1]))
            plain_burden_q = float((plain_burden_sum**2) * ctx.model.rare_variant_scale)

            secure_skat_q = float("nan")
            secure_burden_sum = float("nan")
            secure_burden_q = float("nan")
            if secure is not None and secure.w2s2_vec is not None and secure.w2s2_vec.size > end_idx:
                secure_skat_q = float(np.sum(secure.w2s2_vec[start_idx : end_idx + 1]) * ctx.model.rare_variant_scale)
            if secure is not None and secure.wS_vec is not None and secure.wS_vec.size > end_idx:
                secure_burden_sum = float(np.sum(secure.wS_vec[start_idx : end_idx + 1]))
                secure_burden_q = float((secure_burden_sum**2) * ctx.model.rare_variant_scale)

            window_rows.append(
                WindowComparisonRow(
                    window_index=row_count + 1,
                    block=int(window["block"]),
                    block_window_index=int(window["block_window_index"]),
                    window_start_bp=int(window["window_start_bp"]),
                    window_end_bp=int(window["window_end_bp"]),
                    n_variants=int(window["n_variants"]),
                    start_variant_id=block.variant_ids[start_idx],
                    end_variant_id=block.variant_ids[end_idx],
                    plain_skat_q=plain_skat_q,
                    secure_skat_q=secure_skat_q,
                    skat_abs_diff=abs(plain_skat_q - secure_skat_q) if np.isfinite(secure_skat_q) else math.nan,
                    skat_rel_diff=safe_rel_diff(plain_skat_q, secure_skat_q),
                    plain_burden_q=plain_burden_q,
                    secure_burden_q=secure_burden_q,
                    burden_abs_diff=abs(plain_burden_q - secure_burden_q) if np.isfinite(secure_burden_q) else math.nan,
                    burden_rel_diff=safe_rel_diff(plain_burden_q, secure_burden_q),
                    plain_burden_sum=plain_burden_sum,
                    secure_burden_sum=secure_burden_sum,
                )
            )
            row_count += 1

    return window_rows
