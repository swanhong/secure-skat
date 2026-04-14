# Dataset loading and null-model preparation.

# Read the global SNP position table and verify it matches the expected variant count.
read_positions <- function(party_dir, total_dataset_variants) {
  pos_table <- read.table(
    file.path(party_dir, "snp_pos.txt"),
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  if (ncol(pos_table) < 1) {
    stop(sprintf("No columns found in %s", file.path(party_dir, "snp_pos.txt")))
  }

  pos_col <- if (ncol(pos_table) >= 2) 2L else 1L
  all_positions <- as.integer(pos_table[[pos_col]])
  if (length(all_positions) != total_dataset_variants) {
    stop(sprintf(
      "Position vector length mismatch: expected %d variants but found %d in %s",
      total_dataset_variants,
      length(all_positions),
      file.path(party_dir, "snp_pos.txt")
    ))
  }

  all_positions
}

# Load the secure QC keep mask when the secure pipeline emitted one.
read_qc_filter <- function(secure_run_root, total_dataset_variants) {
  secure_qc_filter_path <- file.path(secure_run_root, "cache", "party1", "gkeep.txt")
  if (!file.exists(secure_qc_filter_path)) {
    return(NULL)
  }

  secure_qc_filter <- readLines(secure_qc_filter_path, warn = FALSE)
  secure_qc_filter <- secure_qc_filter[nzchar(secure_qc_filter)]
  secure_qc_filter <- secure_qc_filter == "1"
  if (length(secure_qc_filter) != total_dataset_variants) {
    stop(sprintf(
      "QC filter length mismatch: expected %d variants but found %d in %s",
      total_dataset_variants,
      length(secure_qc_filter),
      secure_qc_filter_path
    ))
  }

  secure_qc_filter
}

# Read the phenotype vector for one party.
read_pheno <- function(party_dir) {
  as.numeric(read.table(file.path(party_dir, "pheno.txt"), header = FALSE)[, 1])
}

# Read the covariate matrix for one party.
read_cov <- function(party_dir) {
  as.matrix(read.table(file.path(party_dir, "cov.txt"), header = FALSE, sep = "\t"))
}

# Export one chromosome into a dense dosage matrix using PLINK's `.raw` format.
export_chr_matrix <- function(party_dir, chr_index, cache_dir, plink2) {
  out_prefix <- file.path(cache_dir, sprintf("%s_chr%02d", basename(party_dir), chr_index))
  raw_path <- sprintf("%s.raw", out_prefix)

  if (!file.exists(raw_path)) {
    pfile_prefix <- file.path(party_dir, "geno", sprintf("chr%d", chr_index))
    sample_keep <- file.path(party_dir, "sample_keep.txt")
    pfile_modifier <- if (file.exists(sprintf("%s.pvar.zst", pfile_prefix))) {
      " vzs"
    } else {
      ""
    }

    cmd <- sprintf(
      "\"%s\" --pfile \"%s\"%s --keep \"%s\" --export A --out \"%s\"",
      plink2, pfile_prefix, pfile_modifier, sample_keep, out_prefix
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

# Assemble phenotypes, covariates, and null-model residuals shared by all blocks.
load_model_inputs <- function(party_dirs, cov_num_cols = 4L) {
  party_pheno <- lapply(party_dirs, read_pheno)
  party_cov <- lapply(party_dirs, read_cov)
  party_sample_counts <- vapply(party_pheno, length, integer(1))

  y <- unlist(party_pheno, use.names = FALSE)
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

  null_model_s2 <- sum(y_resid^2) / (nrow(design) - ncol(design))

  list(
    party_pheno = party_pheno,
    party_cov = party_cov,
    party_sample_counts = party_sample_counts,
    y = y,
    X = X,
    design = design,
    fit = fit,
    y_resid = y_resid,
    Q_matrix = qr.Q(fit$qr),
    null_model_s2 = null_model_s2,
    skat_package_q_scale = 1.0 / (2.0 * null_model_s2),
    n_total = length(y_resid)
  )
}

# Build the full comparison context by combining path, metadata, and model inputs.
load_compare_context <- function(config) {
  total_dataset_variants <- sum(config$chrom_sizes)

  c(
    config,
    list(
      total_dataset_variants = total_dataset_variants,
      all_positions = read_positions(config$party_dirs[[1]], total_dataset_variants),
      secure_qc_filter = read_qc_filter(config$secure_run_root, total_dataset_variants)
    ),
    load_model_inputs(config$party_dirs)
  )
}
