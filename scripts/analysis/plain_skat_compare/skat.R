# SKAT-specific block reconstruction and package-reference comparisons.

# Compute the SKAT beta weight used for one vector of minor allele frequencies.
compute_beta_weight <- function(maf) {
  25.0 * ((1.0 - maf)^24)
}

# Rebuild all per-block plain SKAT ingredients and retain intermediates for debugging.
compute_plain_skat_blocks <- function(context) {
  q_total_secure_style_raw <- 0.0
  q_total_standard_weight_raw <- 0.0
  burden_total_secure_style_raw <- 0.0
  variant_total <- 0L
  plain_blocks <- vector("list", context$n_blocks)
  package_geno_blocks <- vector("list", context$n_blocks)
  package_weight_blocks <- vector("list", context$n_blocks)

  for (chr_index in seq_len(context$n_blocks)) {
    party_exports <- lapply(
      context$party_dirs,
      export_chr_matrix,
      chr_index = chr_index,
      cache_dir = context$cache_dir,
      plink2 = context$plink2
    )

    variant_ids <- party_exports[[1]]$variant_ids
    if (!identical(variant_ids, party_exports[[2]]$variant_ids)) {
      stop(sprintf("Variant order mismatch at chromosome %d", chr_index))
    }

    local_geno <- lapply(party_exports, `[[`, "geno")
    block_start <- if (chr_index == 1L) 1L else sum(context$chrom_sizes[seq_len(chr_index - 1L)]) + 1L
    block_end <- sum(context$chrom_sizes[seq_len(chr_index)])
    block_positions <- context$all_positions[block_start:block_end]

    if (!is.null(context$secure_qc_filter)) {
      block_filter <- context$secure_qc_filter[block_start:block_end]
      variant_ids <- variant_ids[block_filter]
      local_geno <- lapply(local_geno, function(mat) mat[, block_filter, drop = FALSE])
      block_positions <- block_positions[block_filter]
    }

    local_dosage_sum <- lapply(local_geno, function(mat) {
      mat[is.na(mat)] <- 0.0
      colSums(mat)
    })

    geno <- do.call(rbind, local_geno)
    if (nrow(geno) != context$n_total) {
      stop(sprintf("Sample count mismatch at chromosome %d", chr_index))
    }

    geno[is.na(geno)] <- 0.0

    score <- as.numeric(crossprod(geno, context$y_resid))
    dosage_sum <- colSums(geno)
    p <- dosage_sum / (2.0 * context$n_total)
    p_bar <- 1.0 - p
    maf <- pmin(p, p_bar)
    beta_weight <- compute_beta_weight(maf)

    q_block <- sum((beta_weight^2) * (score^2))
    burden_block <- sum(beta_weight * score)

    q_total_secure_style_raw <- q_total_secure_style_raw + q_block
    q_total_standard_weight_raw <- q_total_standard_weight_raw + sum(beta_weight * (score^2))
    burden_total_secure_style_raw <- burden_total_secure_style_raw + burden_block
    variant_total <- variant_total + ncol(geno)
    package_geno_blocks[[chr_index]] <- geno
    package_weight_blocks[[chr_index]] <- beta_weight

    plain_blocks[[chr_index]] <- list(
      chr_index = chr_index,
      n_variants = ncol(geno),
      variant_ids = variant_ids,
      positions = block_positions,
      local_dosage_sum = local_dosage_sum,
      dosage_sum = dosage_sum,
      dosage_sum_bar = 2.0 * context$n_total - dosage_sum,
      p = p,
      p_bar = p_bar,
      maf = maf,
      score = score,
      score_negated = -score,
      weight = beta_weight,
      score_sq = score^2,
      weight_sq = beta_weight^2,
      weighted_score_sq = (beta_weight^2) * (score^2),
      weighted_score = beta_weight * score,
      weighted_score_negated = -(beta_weight * score),
      q_block = q_block,
      burden_block = burden_block,
      burden_block_negated = -burden_block
    )
  }

  list(
    q_total_secure_style_raw = q_total_secure_style_raw,
    q_total_standard_weight_raw = q_total_standard_weight_raw,
    burden_total_secure_style_raw = burden_total_secure_style_raw,
    variant_total = variant_total,
    plain_blocks = plain_blocks,
    package_geno_blocks = package_geno_blocks,
    package_weight_blocks = package_weight_blocks
  )
}

# Scale raw plain SKAT sums into the final package-compatible totals.
summarize_plain_skat_totals <- function(plain_result, scale_factor, blocks_to_print) {
  list(
    q_total_secure_style = plain_result$q_total_secure_style_raw * scale_factor,
    q_total_standard_weight = plain_result$q_total_standard_weight_raw * scale_factor,
    selected_plain_q_total = sum(vapply(blocks_to_print, function(block_idx) {
      plain_result$plain_blocks[[block_idx]]$q_block
    }, numeric(1))) * scale_factor,
    package_selected_variants = sum(vapply(blocks_to_print, function(block_idx) {
      plain_result$plain_blocks[[block_idx]]$n_variants
    }, integer(1)))
  )
}

# Run the public R SKAT package on the selected blocks for a direct reference result.
run_package_skat_reference <- function(context, plain_result) {
  package_skat_q <- NA_real_
  package_burden_q <- NA_real_
  package_marker_count <- NA_integer_

  if (!context$skip_skat_package && requireNamespace("SKAT", quietly = TRUE)) {
    package_geno <- do.call(cbind, plain_result$package_geno_blocks[context$blocks_to_print])
    package_weights <- unlist(plain_result$package_weight_blocks[context$blocks_to_print], use.names = FALSE)

    cov_df <- as.data.frame(context$X)
    names(cov_df) <- sprintf("X%d", seq_len(ncol(cov_df)))
    null_df <- cbind(data.frame(y = context$y), cov_df)
    null_formula <- if (ncol(cov_df) > 0) {
      as.formula(paste("y ~", paste(names(cov_df), collapse = " + ")))
    } else {
      y ~ 1
    }

    null_obj <- SKAT::SKAT_Null_Model(null_formula, out_type = "C", data = null_df)

    skat_pkg_res <- SKAT::SKAT(
      package_geno,
      null_obj,
      kernel = "linear.weighted",
      method = "davies",
      weights = package_weights,
      is_check_genotype = FALSE,
      is_dosage = TRUE,
      impute.method = "fixed",
      missing_cutoff = 1,
      max_maf = 1,
      estimate_MAF = 1,
      r.corr = 0
    )
    package_skat_q <- as.numeric(skat_pkg_res$Q)
    if (!is.null(skat_pkg_res$param$n.marker.test)) {
      package_marker_count <- as.integer(skat_pkg_res$param$n.marker.test)
    }

    burden_pkg_res <- SKAT::SKAT(
      package_geno,
      null_obj,
      kernel = "linear.weighted",
      method = "davies",
      weights = package_weights,
      is_check_genotype = FALSE,
      is_dosage = TRUE,
      impute.method = "fixed",
      missing_cutoff = 1,
      max_maf = 1,
      estimate_MAF = 1,
      r.corr = 1
    )
    package_burden_q <- as.numeric(burden_pkg_res$Q)
  }

  list(
    package_skat_q = package_skat_q,
    package_burden_q = package_burden_q,
    package_marker_count = package_marker_count
  )
}

# Scale the selected secure per-block SKAT sums into the final comparable Q statistic.
summarize_selected_secure_skat <- function(selected_secure_blocks, secure_scale) {
  if (!isTRUE(selected_secure_blocks$available)) {
    return(list(selected_secure_q = NA_real_))
  }

  list(
    selected_secure_q = selected_secure_blocks$selected_secure_q_raw * secure_scale
  )
}
