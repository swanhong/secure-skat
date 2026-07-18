#!/usr/bin/env Rscript
# R::SKAT ground-truth reference on fed_prep.py / run_fed.sh output (federated-private layout).
# Sibling of skat_compare.R, adapted to the FED_OUT directory layout so it runs unchanged on the
# AoU workbench (blocks format only; no plink2 needed).
#
# For each gene it assembles the POOLED full matrix Z = public list ∪ private variants, with cohort A
# contributing zero columns for the private (B-only) variants, then runs R::SKAT to get the
# ground-truth SKAT and Burden Q/p. This mirrors what the secure pipeline computes.
#
# fed_prep layout (the FED_OUT dir; see scripts/preprocessing/fed_prep.py):
#   A/geno.<g>.bin   n_A x pub_m    public list                     (int8, row-major buf[i*m+j])
#   B/geno.<g>.bin   n_B x pub_m    public list, aligned (public_only cols = 0)
#   B/priv.<g>.bin   n_B x priv_m   B-only private variants         (empty file if priv_m = 0)
#   A/pheno.txt A/cov.txt  B/pheno.txt B/cov.txt                    (cov = N_PCS PCs, tab-separated)
#   block_sizes.txt  pub_m per gene (one per line)
#
# ORIENTATION NOTE: R::SKAT computes the pooled MAF/orientation (the ground truth). The secure
# pipeline uses BLIND LOCAL orientation, which equals the pooled orientation only when the two cohorts
# are concordant for a variant. A per-gene mismatch therefore flags a straddle (see .local/warning.md).
# Private variants use MAF = count/(2N) with cohort A treated as all-reference (the fed_prep spec),
# which is exactly the zero-column assembly below.
#
# Usage:
#   Rscript fed_skat_compare.R <FED_OUT_dir> [secure_csv] [out_dir]
#     <FED_OUT_dir>  run_fed.sh output dir (has A/, B/, block_sizes.txt)
#     [secure_csv]   optional fed_results.csv (per-gene secure p-values) to compare against
#     [out_dir]      default <FED_OUT_dir>/plain_R

suppressMessages(library(SKAT))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 1) stop("usage: fed_skat_compare.R <FED_OUT_dir> [secure_csv] [out_dir]")
fed_out    <- normalizePath(a[[1]], mustWork = TRUE)
secure_csv <- if (length(a) >= 2 && nzchar(a[[2]])) a[[2]] else ""
out_dir    <- if (length(a) >= 3 && nzchar(a[[3]])) a[[3]] else file.path(fed_out, "plain_R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pub_m   <- as.integer(readLines(file.path(fed_out, "block_sizes.txt")))
nblocks <- length(pub_m)

read_pheno <- function(d) as.numeric(read.table(file.path(fed_out, d, "pheno.txt"))[, 1])
read_cov   <- function(d) as.matrix(read.table(file.path(fed_out, d, "cov.txt"), sep = "\t"))
yA <- read_pheno("A"); yB <- read_pheno("B")
nA <- length(yA); nB <- length(yB)

# int8 block loader: row-major n x m (buf[i*m + j]); size 0 / missing file -> n x 0.
# m inferred from file size (1 byte/cell) so no manifest/jsonlite dependency.
load_bin <- function(path, n) {
  if (!file.exists(path) || file.info(path)$size == 0) return(matrix(0, n, 0))
  m <- as.integer(file.info(path)$size / n)
  raw <- readBin(path, "integer", n = n * m, size = 1, signed = TRUE)
  matrix(raw, nrow = n, ncol = m, byrow = TRUE)
}

# Pooled null model (y ~ PCs); SKAT adds its own intercept, drop constant covariate columns.
y <- c(yA, yB)
X <- rbind(read_cov("A"), read_cov("B"))
X <- X[, apply(X, 2, function(col) length(unique(col)) > 1), drop = FALSE]
obj <- SKAT_Null_Model(y ~ X, out_type = "C")

res <- data.frame(gene = integer(0), pub_m = integer(0), priv_m = integer(0),
                  Q_skat = numeric(0), p_skat = numeric(0),
                  Q_burden = numeric(0), p_burden = numeric(0))
for (g in seq_len(nblocks) - 1L) {
  A_pub  <- load_bin(file.path(fed_out, "A", sprintf("geno.%d.bin", g)), nA)
  B_pub  <- load_bin(file.path(fed_out, "B", sprintf("geno.%d.bin", g)), nB)
  B_priv <- load_bin(file.path(fed_out, "B", sprintf("priv.%d.bin", g)), nB)
  vm <- ncol(B_priv)
  # cohort A carries no private variants -> zero columns; pooled full matrix = public list + private.
  Z <- rbind(cbind(A_pub, matrix(0, nA, vm)), cbind(B_pub, B_priv))
  Z[is.na(Z)] <- 0
  if (ncol(Z) == 0) { cat(sprintf("gene %d: empty, skipped\n", g)); next }
  sk <- SKAT(Z, obj, weights.beta = c(1, 25))
  bu <- SKAT(Z, obj, weights.beta = c(1, 25), r.corr = 1)
  res <- rbind(res, data.frame(gene = g, pub_m = ncol(A_pub), priv_m = vm,
                               Q_skat = sk$Q, p_skat = sk$p.value,
                               Q_burden = bu$Q, p_burden = bu$p.value))
  cat(sprintf("gene %d/%d  m=%d(+%d priv)  Q_skat=%.6e p=%.3e | Q_burden=%.6e p=%.3e\n",
              g, nblocks - 1L, ncol(A_pub), vm, sk$Q, sk$p.value, bu$Q, bu$p.value))
}
csv_out <- file.path(out_dir, "fed_plain_R.csv")
write.csv(res, csv_out, row.names = FALSE)
cat(sprintf("wrote %s  (%d genes)   totals: SKAT Q=%.6e  Burden Q=%.6e\n",
            csv_out, nrow(res), sum(res$Q_skat), sum(res$Q_burden)))

# --- optional: compare secure fed_results.csv against this R::SKAT reference ---
if (nzchar(secure_csv) && file.exists(secure_csv)) {
  sec <- read.csv(secure_csv, check.names = FALSE)
  cat(sprintf("\nsecure csv: %s  columns: %s\n", secure_csv, paste(names(sec), collapse = ", ")))
  # best-effort column match; adjust the patterns if fed_results.csv uses other names.
  pick <- function(pat) { i <- grep(pat, names(sec), ignore.case = TRUE); if (length(i)) sec[[i[1]]] else NULL }
  compare1 <- function(label, secval, refval) {
    if (is.null(secval) || length(secval) != nrow(res)) {
      cat(sprintf("  [%s] no matching secure column of length %d -> skip\n", label, nrow(res))); return()
    }
    ss_tot <- sum((refval - mean(refval))^2)
    r2  <- if (ss_tot > 0) 1 - sum((secval - refval)^2) / ss_tot else NA
    rat <- median(secval / refval, na.rm = TRUE)
    cat(sprintf("  [%s vs R::SKAT]  R^2=%.6f  median(secure/R)=%.6g  max|rel|=%.2e\n",
                label, r2, rat, max(abs(secval - refval) / pmax(abs(refval), 1e-12))))
  }
  compare1("SKAT p",   pick("skat.*p|p.*skat"),     res$p_skat)
  compare1("Burden p", pick("burden.*p|p.*burden"), res$p_burden)
}
