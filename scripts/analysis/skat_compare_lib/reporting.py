"""Reporting, CSV writing, and diagnostic console output."""

from __future__ import annotations

import math
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .common import format_float, safe_corr, safe_rel_diff, sanitize_path_tag, trim_or_none
from .models import BlockComparisonRow, CompareContext, ManualBlockResult, ManualResults, ReferenceResult, SecureSummary, WindowComparisonRow
from .secure_io import read_comma_numeric_vector, read_secure_matrix, secure_party_dir


def vector_diff_stats(secure_vec: np.ndarray | None, plain_vec: np.ndarray) -> tuple[float, float, float] | None:
    secure_vec = trim_or_none(secure_vec, plain_vec.size)
    if secure_vec is None:
        return None
    diff = np.abs(secure_vec - plain_vec)
    return float(np.max(diff)), float(np.mean(diff)), safe_corr(secure_vec, plain_vec)


def print_preflight(ctx: CompareContext) -> None:
    print("\n--- Preflight ---")
    print(f"Command: {ctx.command}")
    print(f"Resolved dataset: {ctx.dataset.root}")
    print(f"Run root: {ctx.run_root if ctx.run_root else 'none'}")
    print(f"Blocks: {ctx.dataset.n_blocks}")
    print(f"Variants: {ctx.dataset.total_variants}")
    print(f"Analysis blocks: {ctx.analysis_blocks[0]}..{ctx.analysis_blocks[-1]} ({len(ctx.analysis_blocks)} selected)")
    print(f"Detail blocks: {', '.join(str(block) for block in ctx.detail_blocks)}")
    if ctx.command in {"compare", "reference"}:
        print(f"Reference step: {'skipped' if ctx.skip_reference else 'enabled'}")
    else:
        print("Reference step: n/a for this subcommand")
    print(f"Output directory: {ctx.cache_dir}")
    sys.stdout.flush()


def print_intermediate_diagnostics(ctx: CompareContext, manual: ManualResults) -> None:
    if ctx.run_root is None:
        return

    print("\n--- Intermediate Diagnostics ---")
    q_matrix = ctx.model.q_matrix
    plain_qty_scaled = np.sqrt(ctx.model.n_total) * (q_matrix.T @ ctx.model.y)

    qty_vectors = []
    for party_idx in (1, 2):
        qty_vec = read_comma_numeric_vector(secure_party_dir(ctx.run_root, party_idx) / "qty.txt")
        if qty_vec is not None:
            qty_vectors.append(qty_vec)
    if qty_vectors:
        min_len = min(vec.size for vec in qty_vectors + [plain_qty_scaled])
        secure_qty = np.sum(np.vstack([vec[:min_len] for vec in qty_vectors]), axis=0)
        print(f"Q^T y (scaled) max abs diff: {np.max(np.abs(secure_qty - plain_qty_scaled[:min_len])):.10e}")

    yproj_vectors = []
    for party_idx in (1, 2):
        vec = read_comma_numeric_vector(secure_party_dir(ctx.run_root, party_idx) / "y_proj_rescaled.txt")
        if vec is not None:
            yproj_vectors.append(vec)
    if len(yproj_vectors) == 2:
        plain_yproj = q_matrix @ (q_matrix.T @ ctx.model.y)
        split = ctx.model.party_sample_counts[0]
        plain_p1 = plain_yproj[:split]
        plain_p2 = plain_yproj[split:]
        min_len_p1 = min(plain_p1.size, yproj_vectors[0].size)
        min_len_p2 = min(plain_p2.size, yproj_vectors[1].size)
        print(f"y_proj max abs diff (p1): {np.max(np.abs(yproj_vectors[0][:min_len_p1] - plain_p1[:min_len_p1])):.10e}")
        print(f"y_proj max abs diff (p2): {np.max(np.abs(yproj_vectors[1][:min_len_p2] - plain_p2[:min_len_p2])):.10e}")

    print(f"Samples: {ctx.model.n_total}")
    print(f"Variants after QC alignment: {manual.variant_total}")
    print(f"Null-model rss: {ctx.model.null_rss:.10e}")
    print(f"Null-model dof: {ctx.model.dof}")
    print(f"Null-model s2: {ctx.model.null_model_s2:.10e}")
    print(f"Rare-variant scale dof/(2*rss): {ctx.model.rare_variant_scale:.10e}")
    print(f"Manual SKAT Q (all blocks): {manual.all_skat_q:.10e}")
    print(f"Manual Burden Q (all blocks): {manual.all_burden_q:.10e}")

    secure_ynew = []
    for party_idx in (1, 2):
        vec = read_comma_numeric_vector(secure_party_dir(ctx.run_root, party_idx) / "ynew.txt")
        if vec is None:
            secure_ynew = []
            break
        secure_ynew.append(vec)
    if secure_ynew:
        joined = np.concatenate(secure_ynew)
        min_len = min(joined.size, ctx.model.y_resid.size)
        print(f"Residual max abs diff: {np.max(np.abs(joined[:min_len] - ctx.model.y_resid[:min_len])):.10e}")
        print(f"Residual mean abs diff: {np.mean(np.abs(joined[:min_len] - ctx.model.y_resid[:min_len])):.10e}")
        if joined.size != ctx.model.y_resid.size:
            print(f"Residual length mismatch: plain={ctx.model.y_resid.size}, secure={joined.size}")
    else:
        print("Residuals (ynew.txt) not found in all secure out directories.")

    secure_qcomb = []
    for party_idx in (1, 2):
        mat = read_secure_matrix(secure_party_dir(ctx.run_root, party_idx) / "Qcomb.txt")
        if mat is None:
            secure_qcomb = []
            break
        secure_qcomb.append(mat)
    if secure_qcomb:
        secure_q_rows = np.hstack(secure_qcomb)
        q_plain_scaled = (ctx.model.q_matrix * np.sqrt(ctx.model.n_total)).T
        if secure_q_rows.shape == q_plain_scaled.shape:
            gram_secure = (secure_q_rows @ secure_q_rows.T) / ctx.model.n_total
            gram_plain = (q_plain_scaled @ q_plain_scaled.T) / ctx.model.n_total
            print("\n--- Qcomb Debug ---")
            print("Secure diag(QQ^T / N): " + ", ".join(f"{value:.6f}" for value in np.diag(gram_secure)))
            print("Plain diag(QQ^T / N): " + ", ".join(f"{value:.6f}" for value in np.diag(gram_plain)))
            offdiag = gram_secure - np.diag(np.diag(gram_secure))
            print(f"Secure max |offdiag(QQ^T / N)|: {np.max(np.abs(offdiag)):.10e}")
        else:
            print(
                f"Qcomb dimension mismatch! Secure: {secure_q_rows.shape[0]}x{secure_q_rows.shape[1]}, "
                f"Plain: {q_plain_scaled.shape[0]}x{q_plain_scaled.shape[1]}"
            )


def print_detail_block_summary(block: ManualBlockResult, ctx: CompareContext) -> None:
    secure = block.secure_artifacts
    if secure is None:
        return

    print(f"\n--- Block {block.block_index:02d} ---")
    if secure.q_skat_block_raw is not None:
        print(
            f"qBlock abs diff: {abs(secure.q_skat_block_raw - block.q_skat_block_raw):.10e}, "
            f"rel diff: {safe_rel_diff(secure.q_skat_block_raw, block.q_skat_block_raw):.10e}"
        )
    else:
        print("qBlock missing in secure output")

    if secure.q_burden_block_raw is not None:
        print(
            f"qBurdenBlock abs diff: {abs(secure.q_burden_block_raw - block.q_burden_block_raw):.10e}, "
            f"rel diff: {safe_rel_diff(secure.q_burden_block_raw, block.q_burden_block_raw):.10e}"
        )
    else:
        print("qBurdenBlock missing in secure output")

    if not ctx.debug_mode:
        return

    vector_checks = [
        ("p", secure.p_vec, block.p_vec),
        ("p_bar", secure.p_bar_vec, block.p_bar_vec),
        ("score", secure.score_vec, block.score_vec),
        ("weight", secure.weight_vec, block.weight_vec),
        ("score_sq", secure.score_sq_vec, block.score_sq_vec),
        ("weight_sq", secure.weight_sq_vec, block.weight_sq_vec),
        ("w2S2", secure.w2s2_vec, block.w2s2_vec),
        ("wS", secure.wS_vec, block.wS_vec),
    ]
    for label, secure_vec, plain_vec in vector_checks:
        stats = vector_diff_stats(secure_vec, plain_vec)
        if stats is None:
            print(f"{label}: missing secure output")
            continue
        max_abs, mean_abs, corr = stats
        print(f"{label}: max={max_abs:.3e}, mean={mean_abs:.3e}, corr={corr:.6f}")


def build_variant_debug_frame(block: ManualBlockResult) -> pd.DataFrame:
    secure = block.secure_artifacts
    n = block.n_variants

    def full(vec: np.ndarray | None) -> np.ndarray:
        if vec is None:
            return np.full(n, np.nan)
        return np.asarray(vec[:n], dtype=float)

    return pd.DataFrame(
        {
            "block": np.full(n, block.block_index, dtype=int),
            "variant_index": np.arange(1, n + 1, dtype=int),
            "position": block.positions,
            "variant_id": block.variant_ids,
            "plain_dosage_sum": block.alt_dosage_sum_vec,
            "plain_dosage_sum_bar": block.secure_dosage_sum_vec,
            "secure_dosage_sum": full(None if secure is None else secure.secure_global_dosage_sum),
            "plain_p": block.p_vec,
            "secure_p": full(None if secure is None else secure.p_vec),
            "secure_p_imag": full(None if secure is None else secure.p_vec_imag),
            "plain_p_bar": block.p_bar_vec,
            "secure_p_bar": full(None if secure is None else secure.p_bar_vec),
            "secure_p_bar_imag": full(None if secure is None else secure.p_bar_vec_imag),
            "plain_score": block.score_vec,
            "secure_score": full(None if secure is None else secure.score_vec),
            "plain_score_sq": block.score_sq_vec,
            "secure_score_sq": full(None if secure is None else secure.score_sq_vec),
            "plain_weight": block.weight_vec,
            "secure_weight": full(None if secure is None else secure.weight_vec),
            "secure_weight_imag": full(None if secure is None else secure.weight_vec_imag),
            "plain_weight_sq": block.weight_sq_vec,
            "secure_weight_sq": full(None if secure is None else secure.weight_sq_vec),
            "secure_weight_sq_imag": full(None if secure is None else secure.weight_sq_vec_imag),
            "plain_w2s2": block.w2s2_vec,
            "secure_w2s2": full(None if secure is None else secure.w2s2_vec),
            "secure_w2s2_imag": full(None if secure is None else secure.w2s2_vec_imag),
            "plain_ws": block.wS_vec,
            "secure_ws": full(None if secure is None else secure.wS_vec),
        }
    )


def block_rows_to_frame(rows: list[BlockComparisonRow]) -> pd.DataFrame:
    return pd.DataFrame([asdict(row) for row in rows])


def window_rows_to_frame(rows: list[WindowComparisonRow]) -> pd.DataFrame:
    return pd.DataFrame([asdict(row) for row in rows])


def write_block_compare_csv(ctx: CompareContext, block_df: pd.DataFrame) -> Path:
    path = ctx.cache_dir / "block_compare.csv"
    block_df.to_csv(path, index=False)
    return path


def write_window_compare_csv(ctx: CompareContext, window_df: pd.DataFrame) -> tuple[Path, str]:
    default_tag = f"window_bp{ctx.window_bp}_step{ctx.step_bp}_minv{ctx.min_window_variants}"
    window_tag = sanitize_path_tag(ctx.window_output_tag) if ctx.window_output_tag else default_tag
    path = ctx.cache_dir / f"window_compare_{window_tag}.csv"
    window_df.to_csv(path, index=False)
    return path, window_tag


def write_variant_debug_csvs(ctx: CompareContext, manual: ManualResults) -> None:
    if not ctx.debug_mode:
        return
    variant_frames: list[pd.DataFrame] = []
    for block_index in ctx.analysis_blocks:
        block = manual.blocks[block_index - 1]
        frame = build_variant_debug_frame(block)
        variant_frames.append(frame)
        frame.to_csv(ctx.csv_dir / f"variant_debug_block{block.block_index:02d}.csv", index=False)
    if variant_frames:
        pd.concat(variant_frames, ignore_index=True).to_csv(ctx.csv_dir / "variant_debug_all.csv", index=False)


def print_block_comparison_summary(
    block_df: pd.DataFrame,
    block_csv_path: Path,
    skat_plot_path: Path,
    has_skat_plot: bool,
    burden_plot_path: Path,
    has_burden_plot: bool,
) -> None:
    print("\n--- Block Comparison ---")
    print(f"Blocks retained: {len(block_df)}")
    print(f"Block summary CSV: {block_csv_path}")
    if has_skat_plot:
        print(f"Block-level SKAT scatter plot: {skat_plot_path}")
    if has_burden_plot:
        print(f"Block-level Burden scatter plot: {burden_plot_path}")
    print(f"Block SKAT corr: {safe_corr(block_df['plain_skat_q'], block_df['secure_skat_q']):.10f}")
    print(f"Block Burden corr: {safe_corr(block_df['plain_burden_q'], block_df['secure_burden_q']):.10f}")


def print_window_comparison_summary(
    ctx: CompareContext,
    window_df: pd.DataFrame,
    window_csv_path: Path,
    skat_plot_path: Path,
    has_skat_plot: bool,
    burden_plot_path: Path,
    has_burden_plot: bool,
) -> None:
    print("\n--- Window Comparison ---")
    print(f"Window definition: {ctx.window_bp} bp, step {ctx.step_bp} bp, min variants {ctx.min_window_variants}")
    print(f"Windows retained: {len(window_df)}")
    print(f"Window summary CSV: {window_csv_path}")
    if has_skat_plot:
        print(f"SKAT scatter plot: {skat_plot_path}")
    if has_burden_plot:
        print(f"Burden scatter plot: {burden_plot_path}")
    print(f"SKAT window corr: {safe_corr(window_df['plain_skat_q'], window_df['secure_skat_q']):.10f}")
    print(f"Burden window corr: {safe_corr(window_df['plain_burden_q'], window_df['secure_burden_q']):.10f}")


def write_summary_csv(
    ctx: CompareContext,
    manual: ManualResults,
    secure: SecureSummary | None,
    reference: ReferenceResult,
) -> Path:
    summary_path = ctx.cache_dir / "summary.csv"
    rows = [
        {
            "metric": "skat",
            "reference": reference.skat_q,
            "manual": manual.analysis_skat_q,
            "secure": float("nan") if secure is None else secure.secure_skat_q_for_summary,
        },
        {
            "metric": "burden",
            "reference": reference.burden_q,
            "manual": manual.analysis_burden_q,
            "secure": float("nan") if secure is None else secure.secure_burden_q_for_summary,
        },
    ]
    summary_df = pd.DataFrame(rows)
    summary_df["abs_diff_manual_vs_reference"] = np.abs(summary_df["manual"] - summary_df["reference"])
    summary_df["rel_diff_manual_vs_reference"] = [
        safe_rel_diff(m, r) for m, r in zip(summary_df["manual"], summary_df["reference"])
    ]
    summary_df["abs_diff_secure_vs_reference"] = np.abs(summary_df["secure"] - summary_df["reference"])
    summary_df["rel_diff_secure_vs_reference"] = [
        safe_rel_diff(s, r) for s, r in zip(summary_df["secure"], summary_df["reference"])
    ]
    summary_df["abs_diff_manual_vs_secure"] = np.abs(summary_df["manual"] - summary_df["secure"])
    summary_df["rel_diff_manual_vs_secure"] = [
        safe_rel_diff(m, s) for m, s in zip(summary_df["manual"], summary_df["secure"])
    ]
    summary_df.to_csv(summary_path, index=False)
    return summary_path


def print_compare_summary(
    manual: ManualResults,
    secure: SecureSummary | None,
    reference: ReferenceResult,
    summary_path: Path,
    png_paths: list[Path],
) -> None:
    print("\n--- Final Summary ---")
    if reference.available:
        print(f"Reference SKAT Q: {reference.skat_q:.10e}")
        print(f"Reference Burden Q: {reference.burden_q:.10e}")
        print(f"Reference markers tested: {reference.n_markers}")
    else:
        print(f"Reference result: {reference.skipped_reason or 'unavailable'}")

    print(f"Manual SKAT Q: {manual.analysis_skat_q:.10e}")
    print(f"Manual Burden Q: {manual.analysis_burden_q:.10e}")

    if secure is not None:
        print(f"Secure SKAT Q: {format_float(secure.secure_skat_q_for_summary)}")
        print(f"Secure Burden Q: {format_float(secure.secure_burden_q_for_summary)}")
        print(f"Secure source: {'run-level scalars' if secure.used_run_scalars else 'selected block sums'}")
        if math.isfinite(secure.secure_run_skat_q) or math.isfinite(secure.secure_run_burden_q):
            print(f"Secure run-level SKAT Q: {format_float(secure.secure_run_skat_q)}")
            print(f"Secure run-level Burden Q: {format_float(secure.secure_run_burden_q)}")
    else:
        print("Secure result: unavailable (no run provided)")

    print(f"Summary CSV: {summary_path}")
    for png_path in png_paths:
        print(f"PNG: {png_path}")


def write_manual_summary(ctx: CompareContext, manual: ManualResults) -> Path:
    path = ctx.cache_dir / "manual_summary.tsv"
    pd.DataFrame(
        [
            {"metric": "skat", "value": manual.analysis_skat_q, "raw_value": manual.analysis_q_skat_raw_total},
            {"metric": "burden", "value": manual.analysis_burden_q, "raw_value": manual.analysis_q_burden_raw_total},
        ]
    ).to_csv(path, sep="\t", index=False)
    return path


def write_secure_summary(ctx: CompareContext, secure: SecureSummary) -> Path:
    path = ctx.cache_dir / "secure_summary.tsv"
    pd.DataFrame(
        [
            {
                "metric": "skat",
                "run_level_q": secure.secure_run_skat_q,
                "selected_q": secure.selected_skat_q,
                "summary_q": secure.secure_skat_q_for_summary,
            },
            {
                "metric": "burden",
                "run_level_q": secure.secure_run_burden_q,
                "selected_q": secure.selected_burden_q,
                "summary_q": secure.secure_burden_q_for_summary,
            },
        ]
    ).to_csv(path, sep="\t", index=False)
    return path
