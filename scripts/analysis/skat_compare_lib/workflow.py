"""Top-level orchestration for compare/manual/secure/reference workflows."""

from __future__ import annotations

from . import compute, dataset_io, plotting, reference, reporting, secure_io
from .common import format_float, safe_corr
from .models import CompareContext, ManualResults


def hydrate_secure_block_artifacts(ctx: CompareContext, manual: ManualResults) -> None:
    for block_index in ctx.analysis_blocks:
        block = manual.blocks[block_index - 1]
        block.secure_artifacts = secure_io.load_secure_block_artifacts(ctx, block.block_index, block.n_variants)


def emit_detail_block_diagnostics(ctx: CompareContext, manual: ManualResults) -> None:
    for block_index in ctx.analysis_blocks:
        block = manual.blocks[block_index - 1]
        if block.block_index in ctx.detail_blocks:
            reporting.print_detail_block_summary(block, ctx)


def emit_block_outputs(ctx: CompareContext, manual: ManualResults) -> list:
    block_rows = compute.build_block_comparison_rows(ctx, manual)
    block_df = reporting.block_rows_to_frame(block_rows)
    block_csv_path = reporting.write_block_compare_csv(ctx, block_df)

    skat_plot_path = ctx.cache_dir / "block_compare_skat_scatter.png"
    burden_plot_path = ctx.cache_dir / "block_compare_burden_scatter.png"
    has_skat_plot = plotting.draw_scatter_png(
        block_df["plain_skat_q"].to_numpy(dtype=float),
        block_df["secure_skat_q"].to_numpy(dtype=float),
        skat_plot_path,
        title="Block-Level SKAT Comparison",
        subtitle=f"n = {len(block_df)}, corr = {safe_corr(block_df['plain_skat_q'], block_df['secure_skat_q']):.6f}",
    )
    has_burden_plot = plotting.draw_scatter_png(
        block_df["plain_burden_q"].to_numpy(dtype=float),
        block_df["secure_burden_q"].to_numpy(dtype=float),
        burden_plot_path,
        title="Block-Level Burden Comparison",
        subtitle=f"n = {len(block_df)}, corr = {safe_corr(block_df['plain_burden_q'], block_df['secure_burden_q']):.6f}",
    )
    reporting.print_block_comparison_summary(
        block_df,
        block_csv_path,
        skat_plot_path,
        has_skat_plot,
        burden_plot_path,
        has_burden_plot,
    )
    png_paths = []
    if has_skat_plot:
        png_paths.append(skat_plot_path)
    if has_burden_plot:
        png_paths.append(burden_plot_path)
    return png_paths


def emit_window_outputs(ctx: CompareContext, manual: ManualResults) -> list:
    if ctx.window_bp is None:
        return []

    window_rows = compute.build_window_comparison_rows(ctx, manual)
    if not window_rows:
        print(
            "\n--- Window Comparison ---\n"
            f"No windows met the requested definition ({ctx.window_bp} bp, "
            f"step {ctx.step_bp} bp, min variants {ctx.min_window_variants})."
        )
        return []

    window_df = reporting.window_rows_to_frame(window_rows)
    window_csv_path, window_tag = reporting.write_window_compare_csv(ctx, window_df)
    skat_plot_path = ctx.cache_dir / f"window_compare_{window_tag}_skat_scatter.png"
    burden_plot_path = ctx.cache_dir / f"window_compare_{window_tag}_burden_scatter.png"
    has_skat_plot = plotting.draw_scatter_png(
        window_df["plain_skat_q"].to_numpy(dtype=float),
        window_df["secure_skat_q"].to_numpy(dtype=float),
        skat_plot_path,
        title=f"Windowed SKAT Comparison ({ctx.window_bp} bp, step {ctx.step_bp} bp)",
        subtitle=f"n = {len(window_df)}, corr = {safe_corr(window_df['plain_skat_q'], window_df['secure_skat_q']):.6f}",
    )
    has_burden_plot = plotting.draw_scatter_png(
        window_df["plain_burden_q"].to_numpy(dtype=float),
        window_df["secure_burden_q"].to_numpy(dtype=float),
        burden_plot_path,
        title=f"Windowed Burden Comparison ({ctx.window_bp} bp, step {ctx.step_bp} bp)",
        subtitle=f"n = {len(window_df)}, corr = {safe_corr(window_df['plain_burden_q'], window_df['secure_burden_q']):.6f}",
    )
    reporting.print_window_comparison_summary(
        ctx,
        window_df,
        window_csv_path,
        skat_plot_path,
        has_skat_plot,
        burden_plot_path,
        has_burden_plot,
    )
    png_paths = []
    if has_skat_plot:
        png_paths.append(skat_plot_path)
    if has_burden_plot:
        png_paths.append(burden_plot_path)
    return png_paths


def run_compare(ctx: CompareContext) -> int:
    reporting.print_preflight(ctx)
    dataset_inputs = dataset_io.load_dataset_inputs(ctx)
    manual = compute.compute_manual_results(ctx, dataset_inputs)
    secure = secure_io.load_secure_summary(ctx)
    reporting.print_intermediate_diagnostics(ctx, manual)
    hydrate_secure_block_artifacts(ctx, manual)
    emit_detail_block_diagnostics(ctx, manual)
    reporting.write_variant_debug_csvs(ctx, manual)
    png_paths = emit_block_outputs(ctx, manual)
    png_paths.extend(emit_window_outputs(ctx, manual))
    reference_result = reference.run_reference(ctx, manual)
    summary_path = reporting.write_summary_csv(ctx, manual, secure, reference_result)
    reporting.print_compare_summary(manual, secure, reference_result, summary_path, png_paths)
    return 0


def run_manual(ctx: CompareContext) -> int:
    reporting.print_preflight(ctx)
    dataset_inputs = dataset_io.load_dataset_inputs(ctx)
    manual = compute.compute_manual_results(ctx, dataset_inputs)
    summary_path = reporting.write_manual_summary(ctx, manual)
    print("\n--- Manual Summary ---")
    print(f"Manual SKAT Q: {manual.analysis_skat_q:.10e}")
    print(f"Manual Burden Q: {manual.analysis_burden_q:.10e}")
    print(f"Manual summary TSV: {summary_path}")
    return 0


def run_reference_only(ctx: CompareContext) -> int:
    reporting.print_preflight(ctx)
    dataset_inputs = dataset_io.load_dataset_inputs(ctx)
    manual = compute.compute_manual_results(ctx, dataset_inputs)
    reference_result = reference.run_reference(ctx, manual)
    if reference_result.available:
        print("\n--- Reference Summary ---")
        print(f"Reference SKAT Q: {reference_result.skat_q:.10e}")
        print(f"Reference Burden Q: {reference_result.burden_q:.10e}")
        print(f"Reference markers tested: {reference_result.n_markers}")
        print(f"Reference summary TSV: {reference_result.summary_path}")
    else:
        print(f"Reference skipped: {reference_result.skipped_reason}")
    return 0


def run_secure_only(ctx: CompareContext) -> int:
    reporting.print_preflight(ctx)
    secure = secure_io.load_secure_summary(ctx)
    if secure is None:
        raise RuntimeError("Secure mode requires --run-id")
    summary_path = reporting.write_secure_summary(ctx, secure)
    print("\n--- Secure Summary ---")
    print(f"Secure SKAT Q: {format_float(secure.secure_skat_q_for_summary)}")
    print(f"Secure Burden Q: {format_float(secure.secure_burden_q_for_summary)}")
    print(f"Secure summary TSV: {summary_path}")
    return 0
