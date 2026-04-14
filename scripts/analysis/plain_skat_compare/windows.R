# Sliding-window summaries built from the already reconstructed block-level
# plain and secure intermediate values.

# Convert one block's sorted positions into sliding windows over base-pair coordinates.
make_bp_windows <- function(pos_vec, block_idx, window_bp, step_bp, min_variants) {
  if (length(pos_vec) == 0) {
    return(data.frame())
  }

  starts <- seq.int(pos_vec[[1]], pos_vec[[length(pos_vec)]], by = step_bp)
  left_idx <- findInterval(starts - 1L, pos_vec) + 1L
  right_idx <- findInterval(starts + window_bp - 1L, pos_vec)
  n_variants <- pmax(0L, right_idx - left_idx + 1L)

  keep <- left_idx <= length(pos_vec) & right_idx >= left_idx & n_variants >= min_variants
  if (!any(keep)) {
    return(data.frame())
  }

  data.frame(
    block = rep(block_idx, sum(keep)),
    block_window_index = seq_len(length(starts))[keep],
    window_start_bp = starts[keep],
    window_end_bp = starts[keep] + window_bp - 1L,
    start_index = left_idx[keep],
    end_index = right_idx[keep],
    n_variants = n_variants[keep],
    stringsAsFactors = FALSE
  )
}

# Compute a safe relative difference that does not explode when the denominator is tiny.
safe_rel_diff <- function(x, y) {
  abs(x - y) / pmax(abs(y), 1e-12)
}

# Write a log-scale scatter plot comparing plain and secure window summaries.
write_window_scatter_plot <- function(df, plain_col, secure_col, out_path, plot_title) {
  keep <- is.finite(df[[plain_col]]) & is.finite(df[[secure_col]])
  plot_df <- df[keep, , drop = FALSE]
  if (nrow(plot_df) == 0) {
    return(FALSE)
  }

  eps <- 1e-6
  x <- log10(pmax(plot_df[[plain_col]], eps))
  y <- log10(pmax(plot_df[[secure_col]], eps))
  corr_val <- suppressWarnings(cor(plot_df[[plain_col]], plot_df[[secure_col]], use = "complete.obs"))

  png(out_path, width = 1800, height = 1600, res = 200)
  on.exit(dev.off(), add = TRUE)

  par(mar = c(5.5, 5.5, 4.5, 1.5))
  plot(
    x,
    y,
    pch = 19,
    cex = 0.65,
    col = rgb(0.10, 0.35, 0.70, 0.45),
    xlab = sprintf("log10(plain + %.2e)", eps),
    ylab = sprintf("log10(secure + %.2e)", eps),
    main = plot_title,
    sub = sprintf("n = %d, corr = %.6f", nrow(plot_df), corr_val)
  )
  grid(col = "grey90")
  abline(a = 0, b = 1, col = "firebrick", lwd = 2, lty = 2)

  TRUE
}

# Aggregate already-computed per-variant scores into window-level SKAT and burden summaries.
run_window_comparison <- function(context, plain_blocks, secure_results) {
  if (is.na(context$window_bp)) {
    return(invisible(NULL))
  }

  secure_window_scale <- resolve_selected_secure_scale(
    secure_results$secure_scale_global,
    context$skat_package_q_scale
  )

  window_rows <- list()
  row_idx <- 1L

  for (block_idx in seq_along(plain_blocks)) {
    plain_block <- plain_blocks[[block_idx]]
    block_windows <- make_bp_windows(
      plain_block$positions,
      block_idx,
      context$window_bp,
      context$step_bp,
      context$min_window_variants
    )

    if (nrow(block_windows) == 0) {
      next
    }

    if (!is.na(context$window_limit)) {
      remaining <- context$window_limit - (row_idx - 1L)
      if (remaining <= 0L) {
        break
      }
      if (nrow(block_windows) > remaining) {
        block_windows <- block_windows[seq_len(remaining), , drop = FALSE]
      }
    }

    secure_w2s2 <- plain_block$secure_weighted_score_sq
    secure_ws <- plain_block$secure_weighted_score

    for (i in seq_len(nrow(block_windows))) {
      start_idx <- block_windows$start_index[[i]]
      end_idx <- block_windows$end_index[[i]]

      plain_skat_q <- sum(plain_block$weighted_score_sq[start_idx:end_idx]) * context$skat_package_q_scale
      plain_burden_sum <- sum(plain_block$weighted_score[start_idx:end_idx])
      plain_burden_q <- compute_burden_q(plain_burden_sum, context$skat_package_q_scale)

      secure_skat_q <- NA_real_
      if (!is.null(secure_w2s2) && length(secure_w2s2) >= end_idx) {
        secure_skat_q <- sum(secure_w2s2[start_idx:end_idx]) * secure_window_scale
      }

      secure_burden_sum <- NA_real_
      secure_burden_q <- NA_real_
      if (!is.null(secure_ws) && length(secure_ws) >= end_idx) {
        secure_burden_sum <- sum(secure_ws[start_idx:end_idx])
        secure_burden_q <- compute_burden_q(secure_burden_sum, secure_window_scale)
      }

      window_rows[[row_idx]] <- data.frame(
        window_index = row_idx,
        block = block_idx,
        block_window_index = block_windows$block_window_index[[i]],
        window_start_bp = block_windows$window_start_bp[[i]],
        window_end_bp = block_windows$window_end_bp[[i]],
        n_variants = block_windows$n_variants[[i]],
        start_variant_id = plain_block$variant_ids[[start_idx]],
        end_variant_id = plain_block$variant_ids[[end_idx]],
        plain_skat_q = plain_skat_q,
        secure_skat_q = secure_skat_q,
        skat_abs_diff = abs(plain_skat_q - secure_skat_q),
        skat_rel_diff = safe_rel_diff(plain_skat_q, secure_skat_q),
        plain_burden_q = plain_burden_q,
        secure_burden_q = secure_burden_q,
        burden_abs_diff = abs(plain_burden_q - secure_burden_q),
        burden_rel_diff = safe_rel_diff(plain_burden_q, secure_burden_q),
        plain_burden_sum = plain_burden_sum,
        secure_burden_sum = secure_burden_sum,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
    }
  }

  if (length(window_rows) == 0) {
    cat(sprintf(
      "\n--- Window Comparison ---\nNo windows met the requested definition (%d bp, step %d bp, min variants %d).\n",
      context$window_bp,
      context$step_bp,
      context$min_window_variants
    ))
    return(invisible(NULL))
  }

  window_df <- do.call(rbind, window_rows)
  default_window_tag <- sprintf(
    "window_bp%d_step%d_minv%d",
    context$window_bp,
    context$step_bp,
    context$min_window_variants
  )
  window_tag <- if (is.null(context$window_output_tag) || !nzchar(context$window_output_tag)) {
    default_window_tag
  } else {
    gsub("[^A-Za-z0-9._-]+", "_", context$window_output_tag)
  }

  window_csv_path <- file.path(context$cache_dir, sprintf("window_compare_%s.csv", window_tag))
  skat_plot_path <- file.path(context$cache_dir, sprintf("window_compare_%s_skat_scatter.png", window_tag))
  burden_plot_path <- file.path(context$cache_dir, sprintf("window_compare_%s_burden_scatter.png", window_tag))

  write.csv(window_df, window_csv_path, row.names = FALSE)
  has_skat_plot <- write_window_scatter_plot(
    window_df,
    "plain_skat_q",
    "secure_skat_q",
    skat_plot_path,
    sprintf("Windowed SKAT Comparison (%d bp, step %d bp)", context$window_bp, context$step_bp)
  )
  has_burden_plot <- write_window_scatter_plot(
    window_df,
    "plain_burden_q",
    "secure_burden_q",
    burden_plot_path,
    sprintf("Windowed Burden Comparison (%d bp, step %d bp)", context$window_bp, context$step_bp)
  )

  skat_corr <- suppressWarnings(cor(window_df$plain_skat_q, window_df$secure_skat_q, use = "complete.obs"))
  burden_corr <- suppressWarnings(cor(window_df$plain_burden_q, window_df$secure_burden_q, use = "complete.obs"))
  skat_keep <- is.finite(window_df$skat_abs_diff)
  burden_keep <- is.finite(window_df$burden_abs_diff)
  skat_max_abs <- if (any(skat_keep)) max(window_df$skat_abs_diff[skat_keep]) else NA_real_
  burden_max_abs <- if (any(burden_keep)) max(window_df$burden_abs_diff[burden_keep]) else NA_real_
  worst_skat_idx <- if (any(skat_keep)) {
    which(skat_keep)[which.max(window_df$skat_abs_diff[skat_keep])]
  } else {
    NA_integer_
  }
  worst_burden_idx <- if (any(burden_keep)) {
    which(burden_keep)[which.max(window_df$burden_abs_diff[burden_keep])]
  } else {
    NA_integer_
  }

  cat(sprintf("\n--- Window Comparison ---\n"))
  cat(sprintf(
    "Window definition: %d bp, step %d bp, min variants %d\n",
    context$window_bp,
    context$step_bp,
    context$min_window_variants
  ))
  cat(sprintf("Windows retained: %d\n", nrow(window_df)))
  cat(sprintf("Window summary CSV: %s\n", window_csv_path))
  if (has_skat_plot) {
    cat(sprintf("SKAT scatter plot: %s\n", skat_plot_path))
  }
  if (has_burden_plot) {
    cat(sprintf("Burden scatter plot: %s\n", burden_plot_path))
  }
  cat(sprintf("SKAT window corr: %.10f\n", skat_corr))
  cat(sprintf("SKAT window max abs diff: %.10e\n", skat_max_abs))
  if (!is.na(worst_skat_idx)) {
    cat(sprintf(
      "Worst SKAT window: block %d [%d, %d], n=%d\n",
      window_df$block[worst_skat_idx],
      window_df$window_start_bp[worst_skat_idx],
      window_df$window_end_bp[worst_skat_idx],
      window_df$n_variants[worst_skat_idx]
    ))
  }
  cat(sprintf("Burden window corr: %.10f\n", burden_corr))
  cat(sprintf("Burden window max abs diff: %.10e\n", burden_max_abs))
  if (!is.na(worst_burden_idx)) {
    cat(sprintf(
      "Worst Burden window: block %d [%d, %d], n=%d\n",
      window_df$block[worst_burden_idx],
      window_df$window_start_bp[worst_burden_idx],
      window_df$window_end_bp[worst_burden_idx],
      window_df$n_variants[worst_burden_idx]
    ))
  }

  invisible(window_df)
}
