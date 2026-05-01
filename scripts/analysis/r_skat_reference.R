#!/usr/bin/env Rscript

# Tiny R bridge for invoking the public SKAT package from the Python compare
# pipeline. It expects pre-exported PLINK `.raw` block files plus pheno/cov/weight
# tables already written by Python.

pop_flag_value <- function(args, flag) {
  flag_idx <- match(flag, args)
  if (is.na(flag_idx)) {
    return(list(args = args, value = NULL))
  }
  if (flag_idx == length(args)) {
    stop(sprintf("Missing value for %s", flag))
  }

  list(
    args = args[-c(flag_idx, flag_idx + 1L)],
    value = args[[flag_idx + 1L]]
  )
}

parse_args <- function(args) {
  manifest_res <- pop_flag_value(args, "--manifest")
  args <- manifest_res$args
  pheno_res <- pop_flag_value(args, "--pheno")
  args <- pheno_res$args
  cov_res <- pop_flag_value(args, "--cov")
  args <- cov_res$args
  weights_res <- pop_flag_value(args, "--weights")
  args <- weights_res$args
  out_res <- pop_flag_value(args, "--out")
  args <- out_res$args

  if (length(args) != 0) {
    stop(sprintf("Unexpected trailing arguments: %s", paste(args, collapse = " ")))
  }

  required <- list(
    manifest = manifest_res$value,
    pheno = pheno_res$value,
    cov = cov_res$value,
    weights = weights_res$value,
    out = out_res$value
  )
  missing <- names(required)[vapply(required, is.null, logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
  }

  required
}

read_raw_matrix <- function(path, keep_ids_path) {
  raw_df <- read.table(path, header = TRUE, check.names = FALSE)
  if (ncol(raw_df) <= 6L) {
    stop(sprintf("No genotype columns found in %s", path))
  }

  keep_ids <- readLines(keep_ids_path, warn = FALSE)
  keep_ids <- trimws(keep_ids)
  keep_ids <- keep_ids[nzchar(keep_ids)]
  geno_cols <- colnames(raw_df)[-(1:6)]
  keep_idx <- match(keep_ids, geno_cols)
  if (any(is.na(keep_idx))) {
    missing_ids <- keep_ids[is.na(keep_idx)]
    stop(sprintf(
      "Failed to match %d reference variant ids in %s (example: %s)",
      length(missing_ids),
      path,
      missing_ids[[1]]
    ))
  }

  geno <- as.matrix(raw_df[, keep_idx + 6L, drop = FALSE])
  storage.mode(geno) <- "numeric"
  geno
}

main <- function() {
  if (!requireNamespace("SKAT", quietly = TRUE)) {
    message("Missing R package 'SKAT'. Install it or rerun the Python pipeline with --skip-reference.")
    quit(save = "no", status = 2L)
  }

  args <- parse_args(commandArgs(trailingOnly = TRUE))

  manifest_df <- read.delim(args$manifest, stringsAsFactors = FALSE)
  if (nrow(manifest_df) == 0L) {
    stop("Reference manifest is empty")
  }

  geno_blocks <- lapply(seq_len(nrow(manifest_df)), function(idx) {
    keep_ids_path <- manifest_df$variant_ids_path[[idx]]
    p1 <- read_raw_matrix(manifest_df$raw_path_party1[[idx]], keep_ids_path)
    p2 <- read_raw_matrix(manifest_df$raw_path_party2[[idx]], keep_ids_path)
    if (ncol(p1) != ncol(p2)) {
      stop(sprintf("Variant count mismatch between party raw exports for block %s", manifest_df$block[[idx]]))
    }
    rbind(p1, p2)
  })

  geno <- do.call(cbind, geno_blocks)
  geno[is.na(geno)] <- 0.0

  pheno_df <- read.delim(args$pheno, stringsAsFactors = FALSE)
  cov_df <- read.delim(args$cov, stringsAsFactors = FALSE)
  weights_df <- read.delim(args$weights, stringsAsFactors = FALSE)

  y <- pheno_df[[1]]
  weights <- as.numeric(weights_df[[1]])
  if (nrow(geno) != length(y)) {
    stop(sprintf("Sample count mismatch: geno has %d rows but phenotype has %d rows", nrow(geno), length(y)))
  }
  if (ncol(geno) != length(weights)) {
    stop(sprintf("Marker count mismatch: geno has %d columns but weight vector has %d rows", ncol(geno), length(weights)))
  }

  null_df <- cbind(data.frame(y = y), cov_df)
  cov_names <- names(cov_df)
  null_formula <- if (length(cov_names) > 0L) {
    as.formula(paste("y ~", paste(cov_names, collapse = " + ")))
  } else {
    y ~ 1
  }

  null_obj <- SKAT::SKAT_Null_Model(null_formula, out_type = "C", data = null_df)

  skat_res <- SKAT::SKAT(
    geno,
    null_obj,
    kernel = "linear.weighted",
    method = "davies",
    weights = weights,
    is_check_genotype = FALSE,
    is_dosage = TRUE,
    impute.method = "fixed",
    missing_cutoff = 1,
    max_maf = 1,
    estimate_MAF = 1,
    r.corr = 0
  )
  burden_res <- SKAT::SKAT(
    geno,
    null_obj,
    kernel = "linear.weighted",
    method = "davies",
    weights = weights,
    is_check_genotype = FALSE,
    is_dosage = TRUE,
    impute.method = "fixed",
    missing_cutoff = 1,
    max_maf = 1,
    estimate_MAF = 1,
    r.corr = 1
  )

  n_markers <- if (!is.null(skat_res$param$n.marker.test)) {
    as.integer(skat_res$param$n.marker.test)
  } else {
    as.integer(ncol(geno))
  }

  out_df <- data.frame(
    skat_q = as.numeric(skat_res$Q),
    burden_q = as.numeric(burden_res$Q),
    n_markers = n_markers
  )
  write.table(out_df, file = args$out, sep = "\t", row.names = FALSE, quote = FALSE)
}

main()
