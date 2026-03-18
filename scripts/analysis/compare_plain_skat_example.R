args <- commandArgs(trailingOnly = TRUE)
debug_mode <- "--debug" %in% args
args <- args[args != "--debug"]

dataset_root <- "example_data"
dataset_flag_idx <- match("--dataset", args)
if (!is.na(dataset_flag_idx)) {
  if (dataset_flag_idx == length(args)) {
    stop("Missing value for --dataset")
  }
  dataset_root <- args[[dataset_flag_idx + 1]]
  args <- args[-c(dataset_flag_idx, dataset_flag_idx + 1)]
}

repo_root <- if (length(args) >= 1) args[[1]] else "."
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
dataset_root <- if (grepl("^/", dataset_root)) {
  normalizePath(dataset_root, winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(repo_root, dataset_root), winslash = "/", mustWork = TRUE)
}
run_id <- if (length(args) >= 2) args[[2]] else ""
blocks_arg_index <- if (length(args) >= 3) 3 else NA_integer_
blocks_to_print <- if (!is.na(blocks_arg_index)) {
  as.integer(strsplit(args[[blocks_arg_index]], ",", fixed = TRUE)[[1]])
} else {
  c(1L, 22L)
}
blocks_to_print <- sort(unique(blocks_to_print[!is.na(blocks_to_print) & blocks_to_print >= 1L & blocks_to_print <= 22L]))
if (length(blocks_to_print) == 0) {
  blocks_to_print <- c(1L, 22L)
}

party_dirs <- c(
  file.path(dataset_root, "party1"),
  file.path(dataset_root, "party2")
)

for (dir_path in party_dirs) {
  if (!dir.exists(dir_path)) {
    stop(sprintf("Missing data directory: %s", dir_path))
  }
}

plink2 <- Sys.which("plink2")
if (plink2 == "") {
  stop("plink2 not found in PATH")
}

cache_dir <- file.path(repo_root, ".local", "tmp", "plain_skat_compare")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
csv_dir <- file.path(cache_dir, "variant_debug_csv")
if (debug_mode) {
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
}

secure_party_dir <- function(party_idx) {
  file.path(secure_run_root, sprintf("party%d", party_idx))
}

resolve_secure_run_root <- function() {
  if (!nzchar(run_id)) {
    return(file.path(repo_root, "out"))
  }

  candidates <- Sys.glob(file.path(repo_root, "out", sprintf("output_*_%s", run_id)))
  candidates <- candidates[dir.exists(candidates)]

  if (length(candidates) == 0) {
    stop(sprintf("No secure output directory found for run id: %s", run_id))
  }

  if (length(candidates) > 1) {
    info <- file.info(candidates)
    candidates <- candidates[order(info$mtime, decreasing = TRUE)]
  }

  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

secure_run_root <- resolve_secure_run_root()

read_pheno <- function(party_dir) {
  as.numeric(read.table(file.path(party_dir, "pheno.txt"), header = FALSE)[, 1])
}

read_cov <- function(party_dir) {
  as.matrix(read.table(file.path(party_dir, "cov.txt"), header = FALSE, sep = "\t"))
}

export_chr_matrix <- function(party_dir, chr_index) {
  out_prefix <- file.path(cache_dir, sprintf("%s_chr%02d", basename(party_dir), chr_index))
  raw_path <- sprintf("%s.raw", out_prefix)

  if (!file.exists(raw_path)) {
    pfile_prefix <- file.path(party_dir, "geno", sprintf("chr%d", chr_index))
    sample_keep <- file.path(party_dir, "sample_keep.txt")

    cmd <- sprintf(
      "\"%s\" --pfile \"%s\" --keep \"%s\" --export A --out \"%s\"",
      plink2, pfile_prefix, sample_keep, out_prefix
    )

    status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (status != 0) {
      stop(sprintf("plink2 export failed for %s chr%d", basename(party_dir), chr_index))
    }
  }

  raw_df <- read.table(raw_path, header = TRUE, check.names = FALSE)
  if (ncol(raw_df) <= 6) {
    stop(sprintf("No genotype columns found in %s", raw_path))
  }

  geno <- as.matrix(raw_df[, -(1:6), drop = FALSE])
  storage.mode(geno) <- "numeric"

  list(
    geno = geno,
    variant_ids = colnames(raw_df)[-(1:6)]
  )
}

party_pheno <- lapply(party_dirs, read_pheno)
party_cov <- lapply(party_dirs, read_cov)
party_sample_counts <- vapply(party_pheno, length, integer(1))

y <- unlist(party_pheno, use.names = FALSE)

cov_num_cols <- 4
X <- do.call(
  rbind,
  lapply(party_cov, function(mat) {
    if (ncol(mat) < cov_num_cols) {
      stop("Covariate file has fewer than 4 columns")
    }
    mat[, seq_len(cov_num_cols), drop = FALSE]
  })
)

design <- cbind("(Intercept)" = 1.0, X)
fit <- lm.fit(x = design, y = y)
y_resid <- fit$residuals
Q_matrix <- qr.Q(fit$qr)

if (any(is.na(y_resid))) {
  stop("Residual calculation produced NA values")
}


n_total <- length(y_resid)
q_total_secure_style <- 0.0
q_total_standard_weight <- 0.0
burden_total_secure_style <- 0.0
variant_total <- 0L
plain_blocks <- vector("list", 22)
all_block_csv <- vector("list", 22)

for (chr_index in seq_len(22)) {
  party_exports <- lapply(party_dirs, export_chr_matrix, chr_index = chr_index)

  variant_ids <- party_exports[[1]]$variant_ids
  if (!identical(variant_ids, party_exports[[2]]$variant_ids)) {
    stop(sprintf("Variant order mismatch at chromosome %d", chr_index))
  }

  local_geno <- lapply(party_exports, `[[`, "geno")
  local_dosage_sum <- lapply(local_geno, function(mat) {
    mat[is.na(mat)] <- 0.0
    colSums(mat)
  })

  geno <- do.call(rbind, local_geno)
  if (nrow(geno) != n_total) {
    stop(sprintf("Sample count mismatch at chromosome %d", chr_index))
  }

  geno[is.na(geno)] <- 0.0

  score <- as.numeric(crossprod(geno, y_resid))
  dosage_sum <- colSums(geno)
  p <- dosage_sum / (2.0 * n_total)
  p_bar <- 1.0 - p
  maf <- pmin(p, p_bar)

  # Keep p and p_bar explicit throughout the script.
  # The SKAT beta weight is Beta(MAF; 1, 25) = 25 * (1-MAF)^24.
  beta_weight <- 25.0 * ((1.0 - maf)^24)

  q_block <- sum((beta_weight^2) * (score^2))
  burden_block <- sum(beta_weight * score)

  q_total_secure_style <- q_total_secure_style + q_block
  q_total_standard_weight <- q_total_standard_weight + sum(beta_weight * (score^2))
  burden_total_secure_style <- burden_total_secure_style + burden_block
  variant_total <- variant_total + ncol(geno)

  plain_blocks[[chr_index]] <- list(
    chr_index = chr_index,
    n_variants = ncol(geno),
    variant_ids = variant_ids,
    local_dosage_sum = local_dosage_sum,
    dosage_sum = dosage_sum,
    dosage_sum_bar = 2.0 * n_total - dosage_sum,
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
    q_block = q_block,
    burden_block = burden_block
  )
}

burden_q_total_secure_style <- burden_total_secure_style^2

secure_q_path <- file.path(secure_party_dir(1), "skat_out.txt")
secure_q <- if (file.exists(secure_q_path)) {
  as.numeric(scan(secure_q_path, quiet = TRUE, nmax = 1))
} else {
  NA_real_
}

secure_burden_path <- file.path(secure_party_dir(1), "burden_out.txt")
secure_burden <- if (file.exists(secure_burden_path)) {
  as.numeric(scan(secure_burden_path, quiet = TRUE, nmax = 1))
} else {
  NA_real_
}

read_secure_matrix <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) {
    return(NULL)
  }

  do.call(
    rbind,
    lapply(lines, function(line) {
      vals <- unlist(strsplit(line, ",", fixed = TRUE))
      vals <- vals[vals != ""]
      as.numeric(vals)
    })
  )
}

read_secure_vector <- function(path) {
  mat <- read_secure_matrix(path)
  if (is.null(mat)) {
    return(NULL)
  }

  as.numeric(t(mat))
}

read_secure_vector_any <- function(paths) {
  for (path in paths) {
    vec <- read_secure_vector(path)
    if (!is.null(vec)) {
      return(vec)
    }
  }
  NULL
}

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

trim_or_null <- function(vec, target_len) {
  if (is.null(vec)) {
    return(NULL)
  }
  if (length(vec) < target_len) {
    return(NULL)
  }
  vec[seq_len(target_len)]
}

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

scalar_diff_stats <- function(secure_val, plain_val) {
  list(
    secure = secure_val,
    plain = plain_val,
    abs_diff = abs(secure_val - plain_val),
    rel_diff = abs(secure_val - plain_val) / max(abs(plain_val), 1e-12)
  )
}

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

print_scalar_comparison <- function(block_idx, label, stats) {
  cat(sprintf("Block %02d %-18s\n", block_idx, label))
  cat(sprintf("abs = %.3e, rel = %.3e\n", stats$abs_diff, stats$rel_diff))
  cat(sprintf("plain %s = %.6e\n", label, stats$plain))
  cat(sprintf("secure %s = %.6e\n", label, stats$secure))
}

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


cat(sprintf("\n--- Intermediate Value Comparisons ---\n"))
cat(sprintf("Selected blocks for detailed output: %s\n", paste(blocks_to_print, collapse = ", ")))
plain_qty <- as.vector(crossprod(Q_matrix, y))

# Q^T y check
qty_files <- vapply(seq_along(party_dirs), function(party_idx) {
  file.path(secure_party_dir(party_idx), "qty.txt")
}, character(1))
sys_qty_list <- lapply(qty_files, function(f) {
  if (file.exists(f)) {
    vals <- as.numeric(strsplit(readLines(f)[1], ",\\s*")[[1]])
    vals[!is.na(vals)]
  } else {
    NULL
  }
})
sys_qty_list <- Filter(Negate(is.null), sys_qty_list)

if (length(sys_qty_list) > 0) {
  secure_qty <- rowSums(do.call(cbind, sys_qty_list))

  cat(sprintf("Q^T y max abs diff: %.10e\n", max(abs(secure_qty - plain_qty))))
}

# y_proj_rescaled check (before subtracting from y)
yproj_files <- vapply(seq_along(party_dirs), function(party_idx) {
  file.path(secure_party_dir(party_idx), "y_proj_rescaled.txt")
}, character(1))
sys_yproj_list <- lapply(yproj_files, function(f) {
  if (file.exists(f)) {
    vals <- as.numeric(strsplit(readLines(f)[1], ",\\s*")[[1]])
    vals[!is.na(vals)]
  } else {
    NULL
  }
})
sys_yproj_list <- Filter(Negate(is.null), sys_yproj_list)

if (length(sys_yproj_list) == 2) {
  secure_yproj_p1 <- as.vector(sys_yproj_list[[1]])
  secure_yproj_p2 <- as.vector(sys_yproj_list[[2]])

  plain_yproj <- as.vector(Q_matrix %*% plain_qty)
  plain_yproj_p1 <- plain_yproj[1:1000]
  plain_yproj_p2 <- plain_yproj[1001:2000]

  cat(sprintf("y_proj max abs diff (p1): %.10e\n", max(abs(secure_yproj_p1 - plain_yproj_p1))))
  cat(sprintf("y_proj max abs diff (p2): %.10e\n", max(abs(secure_yproj_p2 - plain_yproj_p2))))
}

# Compare residuals (ynew.txt)
read_secure_ynew <- function(party_dir) {
  party_idx <- match(basename(party_dir), basename(party_dirs))
  ynew_path <- file.path(secure_party_dir(party_idx), "ynew.txt")
  if (file.exists(ynew_path)) {
    lines <- readLines(ynew_path, warn = FALSE)
    if (length(lines) > 0) {
      # Values are comma-separated in the file
      vals_str <- unlist(strsplit(lines, ",", fixed = TRUE))
      # Remove empty strings from trailing commas
      vals_str <- vals_str[vals_str != ""]
      as.numeric(vals_str)
    } else {
      NULL
    }
  } else {
    NULL
  }
}

cat(sprintf("Samples: %d\n", n_total))
cat(sprintf("Variants: %d\n", variant_total))
cat(sprintf("Plain Q (secure-style weights^2): %.10e\n", q_total_secure_style))
cat(sprintf("Plain Q (single beta weight): %.10e\n", q_total_standard_weight))
cat(sprintf("Plain Burden (secure-style): %.10e\n", burden_q_total_secure_style))

secure_ynew_list <- lapply(party_dirs, read_secure_ynew)
if (!any(sapply(secure_ynew_list, is.null))) {
  secure_ynew <- unlist(secure_ynew_list, use.names = FALSE)
  if (length(secure_ynew) == length(y_resid)) {
    max_resid_diff <- max(abs(y_resid - secure_ynew))
    mean_resid_diff <- mean(abs(y_resid - secure_ynew))

    cat(sprintf("Residuals (ynew.txt) max abs diff: %.10e\n", max_resid_diff))
    cat(sprintf("Residuals (ynew.txt) mean abs diff: %.10e\n", mean_resid_diff))
  } else {
    cat(sprintf("Residuals length mismatch! Plain: %d, Secure: %d\n", length(y_resid), length(secure_ynew)))
  }
} else {
  cat("Residuals (ynew.txt) not found in all out directories.\n")
}

secure_qcomb_list <- lapply(
  seq_along(party_dirs),
  function(party_idx) read_secure_matrix(file.path(secure_party_dir(party_idx), "Qcomb.txt"))
)
if (!any(sapply(secure_qcomb_list, is.null))) {
  secure_qcomb <- do.call(cbind, secure_qcomb_list)
  summarize_q_debug(secure_qcomb, y, design)
} else {
  cat("Qcomb.txt not found in all out directories.\n")
}

cat(sprintf("\n--- Block-Level Intermediate Comparisons ---\n"))
for (block_idx in seq_along(plain_blocks)) {
  if (!(block_idx %in% blocks_to_print)) {
    # still export CSV for all blocks even if they are not printed
  }

  plain_block <- plain_blocks[[block_idx]]
  secure_block_idx <- block_idx - 1L
  n_variants <- plain_block$n_variants
  print_block <- block_idx %in% blocks_to_print

  secure_local_dosage <- lapply(
    seq_along(party_dirs),
    function(party_idx) {
      trim_or_null(
        read_secure_vector(file.path(
          secure_party_dir(party_idx),
          sprintf("assoc_cache_dos_sum.skat.%d.txt", secure_block_idx)
        )),
        n_variants
      )
    }
  )

  if (print_block) {
    for (party_idx in seq_along(party_dirs)) {
      plain_local <- plain_block$local_dosage_sum[[party_idx]]
      complement_local <- 2.0 * party_sample_counts[[party_idx]] - plain_local
      print_vector_comparison_primary(
        block_idx,
        sprintf("dosageSum_%s", basename(party_dirs[[party_idx]])),
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
      file.path(secure_party_dir(1), sprintf("p_block%d.txt", secure_block_idx)),
      file.path(secure_party_dir(1), sprintf("p_enc_block%d.txt", secure_block_idx))
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
      file.path(secure_party_dir(1), sprintf("p_bar_block%d.txt", secure_block_idx))
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
    read_secure_vector(file.path(secure_party_dir(1), sprintf("S_vec_block%d.txt", secure_block_idx))),
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
    read_secure_vector(file.path(secure_party_dir(1), sprintf("w_enc_block%d.txt", secure_block_idx))),
    n_variants
  )

  secure_p_csv <- trim_or_null(
    read_secure_vector_any(c(
      file.path(secure_party_dir(1), sprintf("p_block%d.txt", secure_block_idx)),
      file.path(secure_party_dir(1), sprintf("p_enc_block%d.txt", secure_block_idx))
    )),
    n_variants
  )
  secure_p_imag_csv <- trim_or_null(
    read_secure_vector_any(c(
      file.path(secure_party_dir(1), sprintf("p_block%d_imag.txt", secure_block_idx)),
      file.path(secure_party_dir(1), sprintf("p_enc_block%d_imag.txt", secure_block_idx))
    )),
    n_variants
  )
  secure_p_bar_csv <- trim_or_null(
    read_secure_vector_any(c(
      file.path(secure_party_dir(1), sprintf("p_bar_block%d.txt", secure_block_idx))
    )),
    n_variants
  )
  secure_p_bar_imag_csv <- trim_or_null(
    read_secure_vector_any(c(
      file.path(secure_party_dir(1), sprintf("p_bar_block%d_imag.txt", secure_block_idx))
    )),
    n_variants
  )
  secure_weight_csv <- secure_weight
  secure_weight_imag_csv <- trim_or_null(
    read_secure_vector(file.path(secure_party_dir(1), sprintf("w_enc_block%d_imag.txt", secure_block_idx))),
    n_variants
  )
  secure_global_dosage_csv <- if (!any(sapply(secure_local_dosage, is.null))) {
    Reduce(`+`, secure_local_dosage)
  } else {
    rep(NA_real_, n_variants)
  }

  block_csv <- data.frame(
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
    secure_weight = if (is.null(secure_weight_csv)) rep(NA_real_, n_variants) else secure_weight_csv,
    secure_weight_imag = if (is.null(secure_weight_imag_csv)) rep(NA_real_, n_variants) else secure_weight_imag_csv,
    stringsAsFactors = FALSE
  )
  all_block_csv[[block_idx]] <- block_csv
  if (debug_mode) {
    write.csv(
      block_csv,
      file.path(csv_dir, sprintf("variant_debug_block%02d.csv", block_idx)),
      row.names = FALSE
    )
  }

  if (!print_block) {
    next
  }
  print_complex_summary(
    block_idx,
    "p",
    secure_p,
    secure_p_imag_csv
  )
  print_complex_summary(
    block_idx,
    "p_bar",
    secure_p_bar,
    secure_p_bar_imag_csv
  )
  print_vector_comparison(
    block_idx,
    "weight",
    vector_diff_stats(secure_weight, plain_block$weight)
  )

  secure_score_sq <- trim_or_null(
    read_secure_vector(file.path(secure_party_dir(1), sprintf("S2_block%d.txt", secure_block_idx))),
    n_variants
  )
  print_vector_comparison(
    block_idx,
    "score_sq",
    vector_diff_stats(secure_score_sq, plain_block$score_sq)
  )

  secure_weight_sq <- trim_or_null(
    read_secure_vector(file.path(secure_party_dir(1), sprintf("w2_block%d.txt", secure_block_idx))),
    n_variants
  )
  print_vector_comparison(
    block_idx,
    "weight_sq",
    vector_diff_stats(secure_weight_sq, plain_block$weight_sq)
  )

  secure_weighted_score_sq <- trim_or_null(
    read_secure_vector(file.path(secure_party_dir(1), sprintf("w2S2_block%d.txt", secure_block_idx))),
    n_variants
  )
  print_vector_comparison(
    block_idx,
    "w2S2",
    vector_diff_stats(secure_weighted_score_sq, plain_block$weighted_score_sq)
  )

  secure_weighted_score <- trim_or_null(
    read_secure_vector(file.path(secure_party_dir(1), sprintf("wS_block%d.txt", secure_block_idx))),
    n_variants
  )
  print_vector_comparison(
    block_idx,
    "wS",
    vector_diff_stats(secure_weighted_score, plain_block$weighted_score)
  )

  secure_q_block <- read_secure_vector(file.path(secure_party_dir(1), sprintf("qBlock_block%d.txt", secure_block_idx)))
  secure_q_block <- trim_or_null(secure_q_block, 1)
  if (!is.null(secure_q_block)) {
    print_scalar_comparison(
      block_idx,
      "qBlock",
      scalar_diff_stats(secure_q_block[[1]], plain_block$q_block)
    )
  } else {
    cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, "qBlock"))
  }

  secure_burden_block <- read_secure_vector(file.path(secure_party_dir(1), sprintf("qBurdenBlock_block%d.txt", secure_block_idx)))
  secure_burden_block <- trim_or_null(secure_burden_block, 1)
  if (!is.null(secure_burden_block)) {
    print_scalar_comparison(
      block_idx,
      "qBurdenBlock",
      scalar_diff_stats(secure_burden_block[[1]], plain_block$burden_block)
    )
  } else {
    cat(sprintf("Block %02d %-18s missing secure output\n", block_idx, "qBurdenBlock"))
  }

  print_complex_summary(
    block_idx,
    "weight",
    secure_weight,
    trim_or_null(
      read_secure_vector(file.path(secure_party_dir(1), sprintf("w_enc_block%d_imag.txt", secure_block_idx))),
      n_variants
    )
  )
  print_complex_summary(
    block_idx,
    "weight_sq",
    secure_weight_sq,
    trim_or_null(
      read_secure_vector(file.path(secure_party_dir(1), sprintf("w2_block%d_imag.txt", secure_block_idx))),
      n_variants
    )
  )
  print_complex_summary(
    block_idx,
    "w2S2",
    secure_weighted_score_sq,
    trim_or_null(
      read_secure_vector(file.path(secure_party_dir(1), sprintf("w2S2_block%d_imag.txt", secure_block_idx))),
      n_variants
    )
  )
}

if (debug_mode) {
  write.csv(
    do.call(rbind, all_block_csv),
    file.path(csv_dir, "variant_debug_all.csv"),
    row.names = FALSE
  )
}

if (!is.na(secure_q)) {
  cat(sprintf("\n--- SKAT Results ---\n"))
  cat(sprintf("Secure Q (%s): %.10e\n", secure_q_path, secure_q))
  cat(sprintf("Plain Q (secure-style weights^2): %.10e\n", q_total_secure_style))
  cat(sprintf("Absolute difference (secure-style): %.10e\n", abs(q_total_secure_style - secure_q)))
  cat(sprintf("Relative difference (|secure - plain| / secure): %.10e\n", abs(q_total_secure_style - secure_q) / max(abs(secure_q), 1e-12)))
}

if (!is.na(secure_burden)) {
  cat(sprintf("\n--- Burden Results ---\n"))
  cat(sprintf("Secure Burden (%s): %.10e\n", secure_burden_path, secure_burden))
  cat(sprintf("Plain Burden (secure-style): %.10e\n", burden_q_total_secure_style))
  cat(sprintf("Absolute difference (Burden): %.10e\n", abs(burden_q_total_secure_style - secure_burden)))
  cat(sprintf("Relative difference (|secure - plain| / secure): %.10e\n", abs(burden_q_total_secure_style - secure_burden) / max(abs(secure_burden), 1e-12)))
}
