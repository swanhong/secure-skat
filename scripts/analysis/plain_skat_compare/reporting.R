# Comparison helpers and formatted reporting output.

# Return the leading slice of a vector, or NULL when the secure output is missing/short.
trim_or_null <- function(vec, target_len) {
  if (is.null(vec)) {
    return(NULL)
  }
  if (length(vec) < target_len) {
    return(NULL)
  }
  vec[seq_len(target_len)]
}

# Measure how far a secure vector is from the plain reference vector.
vector_diff_stats <- function(secure_vec, plain_vec, alt_plain_vec = NULL) {
  secure_vec <- trim_or_null(secure_vec, length(plain_vec))
  if (is.null(secure_vec)) {
    return(NULL)
  }

  out <- list(
    length = length(plain_vec),
    max_abs = max(abs(secure_vec - plain_vec)),
    mean_abs = mean(abs(secure_vec - plain_vec)),
    corr = suppressWarnings(cor(secure_vec, plain_vec)),
    secure_head = secure_vec[seq_len(min(5, length(secure_vec)))],
    plain_head = plain_vec[seq_len(min(5, length(plain_vec)))]
  )

  if (!is.null(alt_plain_vec)) {
    out$alt_max_abs <- max(abs(secure_vec - alt_plain_vec))
    out$alt_mean_abs <- mean(abs(secure_vec - alt_plain_vec))
    out$alt_corr <- suppressWarnings(cor(secure_vec, alt_plain_vec))
    out$best_match <- if (out$alt_max_abs < out$max_abs) "alt" else "direct"
  }

  out
}

# Compute vector differences while remembering which plain vector is the primary reference.
vector_diff_stats_primary <- function(secure_vec, plain_vec, primary_label, secondary_vec = NULL, secondary_label = NULL) {
  secure_vec <- trim_or_null(secure_vec, length(plain_vec))
  if (is.null(secure_vec)) {
    return(NULL)
  }

  out <- list(
    length = length(plain_vec),
    primary_label = primary_label,
    max_abs = max(abs(secure_vec - plain_vec)),
    mean_abs = mean(abs(secure_vec - plain_vec)),
    corr = suppressWarnings(cor(secure_vec, plain_vec)),
    secure_head = secure_vec[seq_len(min(5, length(secure_vec)))],
    plain_head = plain_vec[seq_len(min(5, length(plain_vec)))]
  )

  if (!is.null(secondary_vec) && !is.null(secondary_label)) {
    out$secondary_label <- secondary_label
    out$secondary_max_abs <- max(abs(secure_vec - secondary_vec))
    out$secondary_mean_abs <- mean(abs(secure_vec - secondary_vec))
    out$secondary_corr <- suppressWarnings(cor(secure_vec, secondary_vec))
    out$secondary_head <- secondary_vec[seq_len(min(5, length(secondary_vec)))]
  }

  out
}

# Compute absolute and relative error for one scalar comparison.
scalar_diff_stats <- function(secure_val, plain_val) {
  list(
    secure = secure_val,
    plain = plain_val,
    abs_diff = abs(secure_val - plain_val),
    rel_diff = abs(secure_val - plain_val) / max(abs(plain_val), 1e-12)
  )
}

# Print magnitude diagnostics for secure complex-valued outputs.
print_complex_summary <- function(block_idx, label, real_vec, imag_vec) {
  if (is.null(real_vec) || is.null(imag_vec)) {
    cat(sprintf("Block %02d %-18s complex missing\n", block_idx, label))
    return(invisible(NULL))
  }

  denom <- pmax(abs(real_vec), 1e-18)
  ratio <- abs(imag_vec) / denom

  cat(sprintf("Block %02d %-18s complex\n", block_idx, label))
  cat(sprintf(
    "max |real| = %.3e, mean |real| = %.3e\n",
    max(abs(real_vec)),
    mean(abs(real_vec))
  ))
  cat(sprintf(
    "max |imag| = %.3e, mean |imag| = %.3e\n",
    max(abs(imag_vec)),
    mean(abs(imag_vec))
  ))
  cat(sprintf(
    "max |imag|/|real| = %.3e, mean |imag|/|real| = %.3e\n",
    max(ratio),
    mean(ratio)
  ))
  cat(sprintf(
    "secure %s real = %s\n",
    label,
    paste(sprintf("%.6e", real_vec[seq_len(min(5, length(real_vec)))]), collapse = ", ")
  ))
  cat(sprintf(
    "secure %s imag = %s\n",
    label,
    paste(sprintf("%.6e", imag_vec[seq_len(min(5, length(imag_vec)))]), collapse = ", ")
  ))
}

# Print a standard secure-vs-plain vector comparison block.
print_vector_comparison <- function(block_idx, label, stats, alt_label = NULL) {
  if (is.null(stats)) {
    cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, label))
    return(invisible(NULL))
  }

  cat(sprintf("Block %02d %-18s\n", block_idx, label))
  cat(sprintf("max = %.3e, mean = %.3e\n", stats$max_abs, stats$mean_abs))
  cat(sprintf(
    "plain %s = %s\n",
    label,
    paste(sprintf("%.6e", stats$plain_head), collapse = ", ")
  ))
  cat(sprintf(
    "secure %s = %s\n",
    label,
    paste(sprintf("%.6e", stats$secure_head), collapse = ", ")
  ))

  if (!is.null(alt_label) && !is.null(stats$alt_max_abs)) {
    cat(sprintf(
      "alt %s max = %.3e, mean = %.3e\n",
      alt_label,
      stats$alt_max_abs,
      stats$alt_mean_abs
    ))
  }
}

# Print a vector comparison where the plain reference label matters for interpretation.
print_vector_comparison_primary <- function(block_idx, label, stats) {
  if (is.null(stats)) {
    cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, label))
    return(invisible(NULL))
  }

  cat(sprintf(
    "Block %02d %-18s [%s]\n",
    block_idx,
    label,
    stats$primary_label
  ))
  cat(sprintf("max = %.3e, mean = %.3e\n", stats$max_abs, stats$mean_abs))
  cat(sprintf(
    "plain %s = %s\n",
    stats$primary_label,
    paste(sprintf("%.6e", stats$plain_head), collapse = ", ")
  ))
  cat(sprintf(
    "secure %s = %s\n",
    label,
    paste(sprintf("%.6e", stats$secure_head), collapse = ", ")
  ))

  if (!is.null(stats$secondary_label)) {
    cat(sprintf(
      "raw %s = %s\n",
      stats$secondary_label,
      paste(sprintf("%.6e", stats$secondary_head), collapse = ", ")
    ))
  }
}

# Print a standard secure-vs-plain scalar comparison block.
print_scalar_comparison <- function(block_idx, label, stats) {
  cat(sprintf("Block %02d %-18s\n", block_idx, label))
  cat(sprintf("abs = %.3e, rel = %.3e\n", stats$abs_diff, stats$rel_diff))
  cat(sprintf("plain %s = %.6e\n", label, stats$plain))
  cat(sprintf("secure %s = %.6e\n", label, stats$secure))
}

# Compare the secure Q basis against the plain QR basis used by the null model.
summarize_q_debug <- function(secure_q_rows, y_vec, design_mat) {
  qr_obj <- qr(design_mat)
  q_plain <- qr.Q(qr_obj, complete = FALSE)
  q_plain_scaled <- t(q_plain * sqrt(length(y_vec)))

  if (!all(dim(secure_q_rows) == dim(q_plain_scaled))) {
    cat(sprintf(
      "Qcomb dimension mismatch! Secure: %s, Plain: %s\n",
      paste(dim(secure_q_rows), collapse = "x"),
      paste(dim(q_plain_scaled), collapse = "x")
    ))
    return(invisible(NULL))
  }

  gram_secure <- (secure_q_rows %*% t(secure_q_rows)) / length(y_vec)
  gram_plain <- (q_plain_scaled %*% t(q_plain_scaled)) / length(y_vec)

  cat(sprintf("\n--- Qcomb Debug ---\n"))
  cat("Secure diag(QQ^T / N): ")
  cat(paste(sprintf("%.6f", diag(gram_secure)), collapse = ", "))
  cat("\n")
  cat("Plain diag(QQ^T / N): ")
  cat(paste(sprintf("%.6f", diag(gram_plain)), collapse = ", "))
  cat("\n")
  cat(sprintf(
    "Secure max |offdiag(QQ^T / N)|: %.10e\n",
    max(abs(gram_secure - diag(diag(gram_secure))))
  ))

  qty_secure <- as.numeric(secure_q_rows %*% y_vec)
  qty_plain <- as.numeric(q_plain_scaled %*% y_vec)
  proj_secure <- as.numeric(crossprod(secure_q_rows, qty_secure)) / length(y_vec)
  proj_plain <- as.numeric(crossprod(q_plain_scaled, qty_plain)) / length(y_vec)

  cat(sprintf("\n--- Projection Debug (QQTy / N) ---\n"))
  cat(sprintf("first 10 proj (secure): %s\n", paste0(head(proj_secure, 10), collapse = ", ")))
  cat(sprintf("first 10 proj (plain): %s\n", paste0(head(proj_plain, 10), collapse = ", ")))
  cat(sprintf("Projection max abs diff: %.10e\n", max(abs(proj_secure - proj_plain))))
  cat(sprintf("Projection mean abs diff: %.10e\n", mean(abs(proj_secure - proj_plain))))
}
