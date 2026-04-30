# Orchestrate the end-to-end comparison by composing the smaller domain modules.

# Print top-level diagnostics that compare secure null-model intermediates to the plain path.
print_intermediate_value_comparisons <- function(context, plain_result, skat_summary, burden_summary) {
  cat(sprintf("\n--- Intermediate Value Comparisons ---\n"))
  cat(sprintf("Selected blocks for detailed output: %s\n", paste(context$blocks_to_print, collapse = ", ")))

  plain_qty <- as.vector(crossprod(context$Q_matrix, context$y))
  plain_qty_scaled <- sqrt(context$n_total) * plain_qty

  qty_files <- vapply(seq_along(context$party_dirs), function(party_idx) {
    file.path(secure_party_dir(context$secure_run_root, party_idx), "qty.txt")
  }, character(1))
  sys_qty_list <- Filter(Negate(is.null), lapply(qty_files, read_comma_numeric_vector))
  if (length(sys_qty_list) > 0) {
    secure_qty <- rowSums(do.call(cbind, sys_qty_list))
    cat(sprintf("Q^T y (scaled) max abs diff: %.10e\n", max(abs(secure_qty - plain_qty_scaled))))
  }

  yproj_files <- vapply(seq_along(context$party_dirs), function(party_idx) {
    file.path(secure_party_dir(context$secure_run_root, party_idx), "y_proj_rescaled.txt")
  }, character(1))
  sys_yproj_list <- Filter(Negate(is.null), lapply(yproj_files, read_comma_numeric_vector))
  if (length(sys_yproj_list) == 2) {
    secure_yproj_p1 <- as.vector(sys_yproj_list[[1]])
    secure_yproj_p2 <- as.vector(sys_yproj_list[[2]])

    plain_yproj <- as.vector(context$Q_matrix %*% plain_qty)
    offset_p1 <- context$party_sample_counts[[1]]
    plain_yproj_p1 <- plain_yproj[seq_len(offset_p1)]
    plain_yproj_p2 <- plain_yproj[(offset_p1 + 1):length(plain_yproj)]

    cat(sprintf("y_proj max abs diff (p1): %.10e\n", max(abs(secure_yproj_p1 - plain_yproj_p1))))
    cat(sprintf("y_proj max abs diff (p2): %.10e\n", max(abs(secure_yproj_p2 - plain_yproj_p2))))
  }

  cat(sprintf("Samples: %d\n", context$n_total))
  cat(sprintf("Variants: %d\n", plain_result$variant_total))
  cat(sprintf("Null-model s2 (SKAT package): %.10e\n", context$null_model_s2))
  cat(sprintf("Q scaling factor 1/(2*s2): %.10e\n", context$skat_package_q_scale))
  cat(sprintf("Plain Q (SKAT package-compatible): %.10e\n", skat_summary$q_total_secure_style))
  cat(sprintf("Plain Q (legacy single beta weight reference): %.10e\n", skat_summary$q_total_standard_weight))
  cat(sprintf("Plain Burden (SKAT package-compatible): %.10e\n", burden_summary$burden_q_total_secure_style))

  secure_ynew_list <- lapply(context$party_dirs, function(party_dir) {
    read_secure_ynew(party_dir, context$party_dirs, context$secure_run_root)
  })
  if (!any(sapply(secure_ynew_list, is.null))) {
    secure_ynew <- unlist(secure_ynew_list, use.names = FALSE)
    if (length(secure_ynew) == length(context$y_resid)) {
      cat(sprintf("Residuals (ynew.txt) max abs diff: %.10e\n", max(abs(context$y_resid - secure_ynew))))
      cat(sprintf("Residuals (ynew.txt) mean abs diff: %.10e\n", mean(abs(context$y_resid - secure_ynew))))
    } else {
      cat(sprintf(
        "Residuals length mismatch! Plain: %d, Secure: %d\n",
        length(context$y_resid),
        length(secure_ynew)
      ))
    }
  } else {
    cat("Residuals (ynew.txt) not found in all out directories.\n")
  }

  secure_qcomb_list <- lapply(seq_along(context$party_dirs), function(party_idx) {
    read_secure_matrix(file.path(secure_party_dir(context$secure_run_root, party_idx), "Qcomb.txt"))
  })
  if (!any(sapply(secure_qcomb_list, is.null))) {
    secure_qcomb <- do.call(cbind, secure_qcomb_list)
    summarize_q_debug(secure_qcomb, context$y, context$design)
  } else {
    cat("Qcomb.txt not found in all out directories.\n")
  }
}

# Compare secure block-level artifacts to plain block-level intermediates and export debug CSVs.
compare_block_level_intermediates <- function(context, plain_blocks) {
  cat(sprintf("\n--- Block-Level Intermediate Comparisons ---\n"))
  all_block_csv <- vector("list", length(plain_blocks))

  for (block_idx in seq_along(plain_blocks)) {
    plain_block <- plain_blocks[[block_idx]]
    secure_block_idx <- block_idx - 1L
    n_variants <- plain_block$n_variants
    print_block <- block_idx %in% context$blocks_to_print

    secure_local_dosage <- lapply(seq_along(context$party_dirs), function(party_idx) {
      trim_or_null(
        read_secure_vector(file.path(
          secure_party_dir(context$secure_run_root, party_idx),
          sprintf("assoc_cache_dos_sum.skat.%d.txt", secure_block_idx)
        )),
        n_variants
      )
    })

    if (print_block) {
      for (party_idx in seq_along(context$party_dirs)) {
        plain_local <- plain_block$local_dosage_sum[[party_idx]]
        complement_local <- 2.0 * context$party_sample_counts[[party_idx]] - plain_local
        print_vector_comparison_primary(
          block_idx,
          sprintf("dosageSum_%s", basename(context$party_dirs[[party_idx]])),
          vector_diff_stats_primary(
            secure_local_dosage[[party_idx]],
            complement_local,
            "plain dosage_sum_bar",
            plain_local,
            "plain dosage_sum"
          )
        )
      }
    }

    secure_p <- trim_or_null(
      read_secure_vector_any(c(
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_block%d.txt", secure_block_idx)),
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_enc_block%d.txt", secure_block_idx))
      )),
      n_variants
    )
    if (print_block && !is.null(secure_p)) {
      print_vector_comparison_primary(
        block_idx,
        "p",
        vector_diff_stats_primary(
          secure_p,
          plain_block$p,
          "plain p",
          plain_block$p_bar,
          "plain p_bar"
        )
      )
    }

    secure_p_bar <- trim_or_null(
      read_secure_vector_any(c(
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_bar_block%d.txt", secure_block_idx))
      )),
      n_variants
    )
    if (print_block && !is.null(secure_p_bar)) {
      print_vector_comparison_primary(
        block_idx,
        "p_bar",
        vector_diff_stats_primary(
          secure_p_bar,
          plain_block$p_bar,
          "plain p_bar",
          plain_block$p,
          "plain p"
        )
      )
    }

    secure_score <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("S_vec_block%d.txt", secure_block_idx)
      )),
      n_variants
    )
    if (print_block) {
      if (!any(sapply(secure_local_dosage, is.null))) {
        secure_global_dosage <- Reduce(`+`, secure_local_dosage)
        print_vector_comparison_primary(
          block_idx,
          "dosageSum_global",
          vector_diff_stats_primary(
            secure_global_dosage,
            plain_block$dosage_sum_bar,
            "plain dosage_sum_bar",
            plain_block$dosage_sum,
            "plain dosage_sum"
          )
        )
      } else {
        cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, "dosageSum_global"))
      }

      print_vector_comparison_primary(
        block_idx,
        "score",
        vector_diff_stats_primary(
          secure_score,
          plain_block$score_negated,
          "plain -score",
          plain_block$score,
          "plain score"
        )
      )
    }

    secure_weight <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("w_enc_block%d.txt", secure_block_idx)
      )),
      n_variants
    )

    secure_p_csv <- trim_or_null(
      read_secure_vector_any(c(
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_block%d.txt", secure_block_idx)),
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_enc_block%d.txt", secure_block_idx))
      )),
      n_variants
    )
    secure_p_imag_csv <- trim_or_null(
      read_secure_vector_any(c(
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_block%d_imag.txt", secure_block_idx)),
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_enc_block%d_imag.txt", secure_block_idx))
      )),
      n_variants
    )
    secure_p_bar_csv <- trim_or_null(
      read_secure_vector_any(c(
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_bar_block%d.txt", secure_block_idx))
      )),
      n_variants
    )
    secure_p_bar_imag_csv <- trim_or_null(
      read_secure_vector_any(c(
        file.path(secure_party_dir(context$secure_run_root, 1L), sprintf("p_bar_block%d_imag.txt", secure_block_idx))
      )),
      n_variants
    )
    secure_weight_imag_csv <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("w_enc_block%d_imag.txt", secure_block_idx)
      )),
      n_variants
    )
    secure_global_dosage_csv <- if (!any(sapply(secure_local_dosage, is.null))) {
      Reduce(`+`, secure_local_dosage)
    } else {
      rep(NA_real_, n_variants)
    }

    all_block_csv[[block_idx]] <- data.frame(
      chr = rep(block_idx, n_variants),
      block = rep(block_idx, n_variants),
      variant_index = seq_len(n_variants),
      variant_id = plain_block$variant_ids,
      plain_dosage_sum = plain_block$dosage_sum,
      plain_dosage_sum_bar = plain_block$dosage_sum_bar,
      secure_dosage_sum = secure_global_dosage_csv,
      plain_p = plain_block$p,
      plain_p_bar = plain_block$p_bar,
      secure_p = if (is.null(secure_p_csv)) rep(NA_real_, n_variants) else secure_p_csv,
      secure_p_imag = if (is.null(secure_p_imag_csv)) rep(NA_real_, n_variants) else secure_p_imag_csv,
      secure_p_bar = if (is.null(secure_p_bar_csv)) rep(NA_real_, n_variants) else secure_p_bar_csv,
      secure_p_bar_imag = if (is.null(secure_p_bar_imag_csv)) rep(NA_real_, n_variants) else secure_p_bar_imag_csv,
      plain_weight = plain_block$weight,
      secure_weight = if (is.null(secure_weight)) rep(NA_real_, n_variants) else secure_weight,
      secure_weight_imag = if (is.null(secure_weight_imag_csv)) rep(NA_real_, n_variants) else secure_weight_imag_csv,
      stringsAsFactors = FALSE
    )

    if (context$debug_mode) {
      write.csv(
        all_block_csv[[block_idx]],
        file.path(context$csv_dir, sprintf("variant_debug_block%02d.csv", block_idx)),
        row.names = FALSE
      )
    }

    if (!print_block) {
      next
    }

    print_complex_summary(block_idx, "p", secure_p, secure_p_imag_csv)
    print_complex_summary(block_idx, "p_bar", secure_p_bar, secure_p_bar_imag_csv)
    print_vector_comparison(block_idx, "weight", vector_diff_stats(secure_weight, plain_block$weight))

    secure_score_sq <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("S2_block%d.txt", secure_block_idx)
      )),
      n_variants
    )
    print_vector_comparison(block_idx, "score_sq", vector_diff_stats(secure_score_sq, plain_block$score_sq))

    secure_weight_sq <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("w2_block%d.txt", secure_block_idx)
      )),
      n_variants
    )
    print_vector_comparison(block_idx, "weight_sq", vector_diff_stats(secure_weight_sq, plain_block$weight_sq))

    secure_weighted_score_sq <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("w2S2_block%d.txt", secure_block_idx)
      )),
      n_variants
    )
    print_vector_comparison(block_idx, "w2S2", vector_diff_stats(secure_weighted_score_sq, plain_block$weighted_score_sq))

    secure_weighted_score <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("wS_block%d.txt", secure_block_idx)
      )),
      n_variants
    )
    plain_blocks[[block_idx]]$secure_weighted_score_sq <- secure_weighted_score_sq
    plain_blocks[[block_idx]]$secure_weighted_score <- secure_weighted_score
    print_vector_comparison(block_idx, "wS", vector_diff_stats(secure_weighted_score, plain_block$weighted_score_negated))

    secure_q_block <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("qBlock_block%d.txt", secure_block_idx)
      )),
      1L
    )
    if (!is.null(secure_q_block)) {
      print_scalar_comparison(
        block_idx,
        "qBlock",
        scalar_diff_stats(secure_q_block[[1]], plain_block$q_block)
      )
    } else {
      cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, "qBlock"))
    }

    secure_burden_block <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("qBurdenBlock_block%d.txt", secure_block_idx)
      )),
      1L
    )
    if (!is.null(secure_burden_block)) {
      print_scalar_comparison(
        block_idx,
        "qBurdenBlock",
        scalar_diff_stats(secure_burden_block[[1]], plain_block$burden_block_negated)
      )
    } else {
      cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, "qBurdenBlock"))
    }

    print_complex_summary(block_idx, "weight", secure_weight, secure_weight_imag_csv)
    print_complex_summary(
      block_idx,
      "weight_sq",
      secure_weight_sq,
      trim_or_null(
        read_secure_vector(file.path(
          secure_party_dir(context$secure_run_root, 1L),
          sprintf("w2_block%d_imag.txt", secure_block_idx)
        )),
        n_variants
      )
    )
    print_complex_summary(
      block_idx,
      "w2S2",
      secure_weighted_score_sq,
      trim_or_null(
        read_secure_vector(file.path(
          secure_party_dir(context$secure_run_root, 1L),
          sprintf("w2S2_block%d_imag.txt", secure_block_idx)
        )),
        n_variants
      )
    )
  }

  if (context$debug_mode) {
    write.csv(
      do.call(rbind, all_block_csv),
      file.path(context$csv_dir, "variant_debug_all.csv"),
      row.names = FALSE
    )
  }

  plain_blocks
}

# Print selected-block comparisons against the optional SKAT package reference outputs.
print_package_function_comparison <- function(context, skat_summary, burden_summary, package_result, selected_secure_skat, selected_secure_burden) {
  cat(sprintf("\n--- SKAT Package Function Comparison ---\n"))
  cat(sprintf("Selected blocks: %s\n", paste(context$blocks_to_print, collapse = ", ")))
  cat(sprintf("Selected variants (plain path): %d\n", skat_summary$package_selected_variants))
  if (!is.na(package_result$package_marker_count)) {
    cat(sprintf("Selected variants tested by SKAT package: %d\n", package_result$package_marker_count))
  }

  if (context$skip_skat_package) {
    cat("SKAT package result: skipped by --skip-skat-package.\n")
    return(invisible(NULL))
  }
  if (is.na(package_result$package_skat_q)) {
    cat("SKAT package result: package 'SKAT' not installed, skipping direct function comparison.\n")
    return(invisible(NULL))
  }

  cat(sprintf("SKAT::SKAT() Q: %.10e\n", package_result$package_skat_q))
  cat(sprintf("Plain selected-block Q: %.10e\n", skat_summary$selected_plain_q_total))
  cat(sprintf(
    "Absolute difference (plain vs SKAT package): %.10e\n",
    abs(skat_summary$selected_plain_q_total - package_result$package_skat_q)
  ))
  cat(sprintf(
    "Relative difference (plain vs SKAT package): %.10e\n",
    abs(skat_summary$selected_plain_q_total - package_result$package_skat_q) /
      max(abs(package_result$package_skat_q), 1e-12)
  ))
  if (!is.na(selected_secure_skat$selected_secure_q)) {
    cat(sprintf("Secure selected-block Q: %.10e\n", selected_secure_skat$selected_secure_q))
    cat(sprintf(
      "Absolute difference (secure vs SKAT package): %.10e\n",
      abs(selected_secure_skat$selected_secure_q - package_result$package_skat_q)
    ))
    cat(sprintf(
      "Relative difference (secure vs SKAT package): %.10e\n",
      abs(selected_secure_skat$selected_secure_q - package_result$package_skat_q) /
        max(abs(package_result$package_skat_q), 1e-12)
    ))
  }

  cat(sprintf("SKAT::SKAT(r.corr=1) burden Q: %.10e\n", package_result$package_burden_q))
  cat(sprintf("Plain selected-block burden Q: %.10e\n", burden_summary$selected_plain_burden_q))
  cat(sprintf(
    "Absolute difference (plain burden vs SKAT package): %.10e\n",
    abs(burden_summary$selected_plain_burden_q - package_result$package_burden_q)
  ))
  cat(sprintf(
    "Relative difference (plain burden vs SKAT package): %.10e\n",
    abs(burden_summary$selected_plain_burden_q - package_result$package_burden_q) /
      max(abs(package_result$package_burden_q), 1e-12)
  ))
  if (!is.na(selected_secure_burden$selected_secure_burden_q)) {
    cat(sprintf("Secure selected-block burden Q: %.10e\n", selected_secure_burden$selected_secure_burden_q))
    cat(sprintf(
      "Absolute difference (secure burden vs SKAT package): %.10e\n",
      abs(selected_secure_burden$selected_secure_burden_q - package_result$package_burden_q)
    ))
    cat(sprintf(
      "Relative difference (secure burden vs SKAT package): %.10e\n",
      abs(selected_secure_burden$selected_secure_burden_q - package_result$package_burden_q) /
        max(abs(package_result$package_burden_q), 1e-12)
    ))
  }
}

# Print final full-run SKAT and burden summaries against the secure outputs.
print_run_level_summaries <- function(skat_summary, burden_summary, secure_results) {
  if (!is.na(secure_results$secure_q)) {
    cat(sprintf("\n--- SKAT Results ---\n"))
    cat(sprintf("Secure Q (%s): %.10e\n", secure_results$secure_q_path, secure_results$secure_q))
    if (!identical(secure_results$secure_q_for_compare, secure_results$secure_q)) {
      cat(sprintf(
        "Secure Q after inferred SKAT-package scaling: %.10e\n",
        secure_results$secure_q_for_compare
      ))
    }
    cat(sprintf("Plain Q (SKAT package-compatible): %.10e\n", skat_summary$q_total_secure_style))
    cat(sprintf(
      "Absolute difference (plain): %.10e\n",
      abs(skat_summary$q_total_secure_style - secure_results$secure_q_for_compare)
    ))
    cat(sprintf(
      "Relative difference (|secure - plain| / secure): %.10e\n",
      abs(skat_summary$q_total_secure_style - secure_results$secure_q_for_compare) /
        max(abs(secure_results$secure_q_for_compare), 1e-12)
    ))
  }

  if (!is.na(secure_results$secure_burden)) {
    cat(sprintf("\n--- Burden Results ---\n"))
    cat(sprintf(
      "Secure Burden (%s): %.10e\n",
      secure_results$secure_burden_path,
      secure_results$secure_burden
    ))
    if (!identical(secure_results$secure_burden_for_compare, secure_results$secure_burden)) {
      cat(sprintf(
        "Secure Burden after inferred SKAT-package scaling: %.10e\n",
        secure_results$secure_burden_for_compare
      ))
    }
    cat(sprintf(
      "Plain Burden (SKAT package-compatible): %.10e\n",
      burden_summary$burden_q_total_secure_style
    ))
    cat(sprintf(
      "Absolute difference (Burden): %.10e\n",
      abs(burden_summary$burden_q_total_secure_style - secure_results$secure_burden_for_compare)
    ))
    cat(sprintf(
      "Relative difference (|secure - plain| / secure): %.10e\n",
      abs(burden_summary$burden_q_total_secure_style - secure_results$secure_burden_for_compare) /
        max(abs(secure_results$secure_burden_for_compare), 1e-12)
    ))
  }
}

# Run the full plain-vs-secure comparison pipeline from CLI args to final summaries.
run_plain_skat_compare <- function(args) {
  context <- load_compare_context(resolve_compare_paths(parse_compare_args(args)))
  plain_result <- compute_plain_skat_blocks(context)
  skat_summary <- summarize_plain_skat_totals(
    plain_result,
    context$skat_package_q_scale,
    context$blocks_to_print
  )
  burden_summary <- summarize_plain_burden_totals(
    plain_result,
    context$skat_package_q_scale,
    context$blocks_to_print
  )
  package_result <- run_package_skat_reference(context, plain_result)
  secure_results <- load_secure_scalar_results(context)

  print_intermediate_value_comparisons(context, plain_result, skat_summary, burden_summary)
  plain_result$plain_blocks <- compare_block_level_intermediates(context, plain_result$plain_blocks)
  run_block_comparison(context, plain_result$plain_blocks, secure_results)
  run_window_comparison(context, plain_result$plain_blocks, secure_results)

  selected_secure_blocks <- read_selected_secure_block_sums(context)
  selected_secure_scale <- resolve_selected_secure_scale(
    secure_results$secure_scale_global,
    context$skat_package_q_scale
  )
  selected_secure_skat <- summarize_selected_secure_skat(selected_secure_blocks, selected_secure_scale)
  selected_secure_burden <- summarize_selected_secure_burden(selected_secure_blocks, selected_secure_scale)

  print_package_function_comparison(
    context,
    skat_summary,
    burden_summary,
    package_result,
    selected_secure_skat,
    selected_secure_burden
  )
  print_run_level_summaries(skat_summary, burden_summary, secure_results)

  invisible(list(
    context = context,
    plain_result = plain_result,
    skat_summary = skat_summary,
    burden_summary = burden_summary,
    package_result = package_result,
    secure_results = secure_results
  ))
}
