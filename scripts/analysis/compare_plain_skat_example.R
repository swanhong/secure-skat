args <- commandArgs(trailingOnly = TRUE)

repo_root <- if (length(args) >= 1) args[[1]] else "."
repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)

party_dirs <- c(
  file.path(repo_root, "example_data", "party1"),
  file.path(repo_root, "example_data", "party2")
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

if (any(is.na(y_resid))) {
  stop("Residual calculation produced NA values")
}


n_total <- length(y_resid)
q_total_secure_style <- 0.0
q_total_standard_weight <- 0.0
burden_total_secure_style <- 0.0
variant_total <- 0L

for (chr_index in seq_len(22)) {
  party_exports <- lapply(party_dirs, export_chr_matrix, chr_index = chr_index)

  variant_ids <- party_exports[[1]]$variant_ids
  if (!identical(variant_ids, party_exports[[2]]$variant_ids)) {
    stop(sprintf("Variant order mismatch at chromosome %d", chr_index))
  }

  geno <- do.call(rbind, lapply(party_exports, `[[`, "geno"))
  if (nrow(geno) != n_total) {
    stop(sprintf("Sample count mismatch at chromosome %d", chr_index))
  }

  geno[is.na(geno)] <- 0.0

  score <- as.numeric(crossprod(geno, y_resid))
  dosage_sum <- colSums(geno)
  alt_af <- dosage_sum / (2.0 * n_total)
  maf <- ifelse(alt_af > 0.5, 1.0 - alt_af, alt_af)

  beta_weight <- 25.0 * (1.0 - maf)^24

  q_block <- sum((beta_weight^2) * (score^2))
  burden_block <- sum(beta_weight * score)

  q_total_secure_style <- q_total_secure_style + q_block
  q_total_standard_weight <- q_total_standard_weight + sum(beta_weight * (score^2))
  burden_total_secure_style <- burden_total_secure_style + burden_block
  variant_total <- variant_total + ncol(geno)
}

burden_q_total_secure_style <- burden_total_secure_style^2

secure_q_path <- file.path(repo_root, "out", "party1", "skat_out.txt")
secure_q <- if (file.exists(secure_q_path)) {
  as.numeric(scan(secure_q_path, quiet = TRUE, nmax = 1))
} else {
  NA_real_
}

secure_burden_path <- file.path(repo_root, "out", "party1", "burden_out.txt")
secure_burden <- if (file.exists(secure_burden_path)) {
  as.numeric(scan(secure_burden_path, quiet = TRUE, nmax = 1))
} else {
  NA_real_
}

cat(sprintf("Samples: %d\n", n_total))
cat(sprintf("Variants: %d\n", variant_total))
cat(sprintf("Plain Q (secure-style weights^2): %.10e\n", q_total_secure_style))
cat(sprintf("Plain Q (single beta weight): %.10e\n", q_total_standard_weight))
cat(sprintf("Plain Burden (secure-style): %.10e\n", burden_q_total_secure_style))

if (!is.na(secure_q)) {
  cat(sprintf("\n--- SKAT Results ---\n"))
  cat(sprintf("Secure Q (out/party1/skat_out.txt): %.10e\n", secure_q))
  cat(sprintf("Absolute difference (secure-style): %.10e\n", abs(q_total_secure_style - secure_q)))
  cat(sprintf("Relative difference (secure-style): %.10e\n", abs(q_total_secure_style - secure_q) / max(abs(secure_q), 1e-12)))
}

if (!is.na(secure_burden)) {
  cat(sprintf("\n--- Burden Results ---\n"))
  cat(sprintf("Secure Burden (out/party1/burden_out.txt): %.10e\n", secure_burden))
  cat(sprintf("Absolute difference (Burden): %.10e\n", abs(burden_q_total_secure_style - secure_burden)))
  cat(sprintf("Relative difference (Burden): %.10e\n", abs(burden_q_total_secure_style - secure_burden) / max(abs(secure_burden), 1e-12)))
}
