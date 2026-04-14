# Block-level summaries, useful when each block already represents a selected
# biological window or locus.

# Export one row per block and draw block-level plain-vs-secure scatter plots.
run_block_comparison <- function(context, plain_blocks, secure_results) {
  secure_block_scale <- resolve_selected_secure_scale(
    secure_results$secure_scale_global,
    context$skat_package_q_scale
  )

  block_rows <- lapply(seq_along(plain_blocks), function(block_idx) {
    plain_block <- plain_blocks[[block_idx]]
    secure_block_idx <- block_idx - 1L

    secure_q_raw <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("qBlock_block%d.txt", secure_block_idx)
      )),
      1L
    )
    secure_burden_raw <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("qBurdenBlock_block%d.txt", secure_block_idx)
      )),
      1L
    )

    data.frame(
      block = block_idx,
      n_variants = plain_block$n_variants,
      block_start_bp = if (plain_block$n_variants > 0L) plain_block$positions[[1]] else NA_integer_,
      block_end_bp = if (plain_block$n_variants > 0L) plain_block$positions[[plain_block$n_variants]] else NA_integer_,
      start_variant_id = if (plain_block$n_variants > 0L) plain_block$variant_ids[[1]] else NA_character_,
      end_variant_id = if (plain_block$n_variants > 0L) plain_block$variant_ids[[plain_block$n_variants]] else NA_character_,
      plain_skat_q = plain_block$q_block * context$skat_package_q_scale,
      secure_skat_q = if (is.null(secure_q_raw)) NA_real_ else secure_q_raw[[1]] * secure_block_scale,
      skat_abs_diff = if (is.null(secure_q_raw)) NA_real_ else abs((plain_block$q_block * context$skat_package_q_scale) - (secure_q_raw[[1]] * secure_block_scale)),
      skat_rel_diff = if (is.null(secure_q_raw)) NA_real_ else safe_rel_diff(
        plain_block$q_block * context$skat_package_q_scale,
        secure_q_raw[[1]] * secure_block_scale
      ),
      plain_burden_q = compute_burden_q(plain_block$burden_block, context$skat_package_q_scale),
      secure_burden_q = if (is.null(secure_burden_raw)) NA_real_ else compute_burden_q(secure_burden_raw[[1]], secure_block_scale),
      burden_abs_diff = if (is.null(secure_burden_raw)) NA_real_ else abs(
        compute_burden_q(plain_block$burden_block, context$skat_package_q_scale) -
          compute_burden_q(secure_burden_raw[[1]], secure_block_scale)
      ),
      burden_rel_diff = if (is.null(secure_burden_raw)) NA_real_ else safe_rel_diff(
        compute_burden_q(plain_block$burden_block, context$skat_package_q_scale),
        compute_burden_q(secure_burden_raw[[1]], secure_block_scale)
      ),
      plain_burden_sum = plain_block$burden_block,
      secure_burden_sum = if (is.null(secure_burden_raw)) NA_real_ else secure_burden_raw[[1]],
      stringsAsFactors = FALSE
    )
  })

  block_df <- do.call(rbind, block_rows)
  block_csv_path <- file.path(context$cache_dir, "block_compare.csv")
  skat_plot_path <- file.path(context$cache_dir, "block_compare_skat_scatter.png")
  burden_plot_path <- file.path(context$cache_dir, "block_compare_burden_scatter.png")

  write.csv(block_df, block_csv_path, row.names = FALSE)
  has_skat_plot <- write_window_scatter_plot(
    block_df,
    "plain_skat_q",
    "secure_skat_q",
    skat_plot_path,
    "Block-Level SKAT Comparison"
  )
  has_burden_plot <- write_window_scatter_plot(
    block_df,
    "plain_burden_q",
    "secure_burden_q",
    burden_plot_path,
    "Block-Level Burden Comparison"
  )

  skat_corr <- suppressWarnings(cor(block_df$plain_skat_q, block_df$secure_skat_q, use = "complete.obs"))
  burden_corr <- suppressWarnings(cor(block_df$plain_burden_q, block_df$secure_burden_q, use = "complete.obs"))

  cat(sprintf("\n--- Block Comparison ---\n"))
  cat(sprintf("Blocks retained: %d\n", nrow(block_df)))
  cat(sprintf("Block summary CSV: %s\n", block_csv_path))
  if (has_skat_plot) {
    cat(sprintf("Block-level SKAT scatter plot: %s\n", skat_plot_path))
  }
  if (has_burden_plot) {
    cat(sprintf("Block-level Burden scatter plot: %s\n", burden_plot_path))
  }
  cat(sprintf("Block SKAT corr: %.10f\n", skat_corr))
  cat(sprintf("Block Burden corr: %.10f\n", burden_corr))

  invisible(block_df)
}
