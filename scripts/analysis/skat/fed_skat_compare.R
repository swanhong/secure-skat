#!/usr/bin/env Rscript
# R::SKAT ground-truth reference on fed_prep.py / run_fed.sh output (federated-private layout),
# with an optional 3-way compare (secure vs Python-plain vs R::SKAT) against fed_results.csv.
# Sibling of skat_compare.R, adapted to the FED_OUT layout so it runs unchanged on the AoU workbench
# (blocks format only; no plink2 needed).
#
# Per gene it assembles the POOLED full matrix Z = public list ∪ private variants (cohort A contributes
# zero columns for the B-only private variants), then runs R::SKAT. For SKAT it reports the p-value by
# BOTH methods -- Davies (exact) and Liu (moment-matching, closest to our Wilson-Hilferty). For Burden
# (r.corr=1) the statistic is 1-dof so there is a single exact p-value.
#
# fed_prep layout (FED_OUT dir; see scripts/preprocessing/fed_prep.py):
#   A/geno.<g>.bin  n_A x pub_m   B/geno.<g>.bin  n_B x pub_m (aligned)   B/priv.<g>.bin  n_B x priv_m
#   {A,B}/pheno.txt {A,B}/cov.txt   block_sizes.txt (pub_m per gene)
# fed_results.csv (from run_fed.sh FED_CSV=1) columns:
#   gene,chrom,pos,burden_p_secure,burden_p_plain[,skat_p_secure,skat_p_plain]
#
# ORIENTATION NOTE: R::SKAT computes the pooled MAF/orientation (ground truth); the secure/plain
# pipeline uses BLIND LOCAL orientation, which equals pooled only when the two cohorts are concordant.
# A per-gene divergence therefore flags a straddle (see .local/warning.md). SKAT p by Davies is the
# exact reference; our secure SKAT p is Wilson-Hilferty, so it should track Liu more closely than Davies
# in the tail (a Davies mismatch there is the approximation, not a bug). Burden p is 1-dof, so
# secure/plain/R should agree directly.
#
# Usage: Rscript fed_skat_compare.R <FED_OUT_dir> [fed_results.csv] [out_dir]

suppressMessages(library(SKAT))

a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 1) stop("usage: fed_skat_compare.R <FED_OUT_dir> [fed_results.csv] [out_dir]")
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

# int8 block loader: row-major n x m (buf[i*m + j]); m inferred from file size (1 byte/cell).
load_bin <- function(path, n) {
  if (!file.exists(path) || file.info(path)$size == 0) return(matrix(0, n, 0))
  m <- as.integer(file.info(path)$size / n)
  matrix(readBin(path, "integer", n = n * m, size = 1, signed = TRUE), nrow = n, ncol = m, byrow = TRUE)
}

y <- c(yA, yB)
X <- rbind(read_cov("A"), read_cov("B"))
X <- X[, apply(X, 2, function(col) length(unique(col)) > 1), drop = FALSE]
obj <- SKAT_Null_Model(y ~ X, out_type = "C")

res <- data.frame()
for (g in seq_len(nblocks) - 1L) {
  A_pub  <- load_bin(file.path(fed_out, "A", sprintf("geno.%d.bin", g)), nA)
  B_pub  <- load_bin(file.path(fed_out, "B", sprintf("geno.%d.bin", g)), nB)
  B_priv <- load_bin(file.path(fed_out, "B", sprintf("priv.%d.bin", g)), nB)
  vm <- ncol(B_priv)
  Z <- rbind(cbind(A_pub, matrix(0, nA, vm)), cbind(B_pub, B_priv))
  Z[is.na(Z)] <- 0
  if (ncol(Z) == 0) { cat(sprintf("gene %d: empty, skipped\n", g)); next }
  sk_d <- SKAT(Z, obj, weights.beta = c(1, 25), method = "davies")
  sk_l <- SKAT(Z, obj, weights.beta = c(1, 25), method = "liu")
  bu   <- SKAT(Z, obj, weights.beta = c(1, 25), r.corr = 1)   # 1-dof burden, exact p
  res <- rbind(res, data.frame(gene = g, pub_m = ncol(A_pub), priv_m = vm,
                               Q_skat = sk_d$Q,
                               skat_p_R_davies = sk_d$p.value, skat_p_R_liu = sk_l$p.value,
                               burden_p_R = bu$p.value))
  cat(sprintf("gene %d/%d  m=%d(+%d)  SKAT p: davies=%.3e liu=%.3e | Burden p=%.3e\n",
              g, nblocks - 1L, ncol(A_pub), vm, sk_d$p.value, sk_l$p.value, bu$p.value))
}
csv_out <- file.path(out_dir, "fed_plain_R.csv")
write.csv(res, csv_out, row.names = FALSE)
cat(sprintf("wrote %s  (%d genes)\n", csv_out, nrow(res)))

# --- 3-way compare against fed_results.csv (secure + Python plain) ---
if (nzchar(secure_csv) && file.exists(secure_csv)) {
  sec <- read.csv(secure_csv, check.names = FALSE)
  m <- merge(sec, res, by = "gene")
  m <- m[order(m$gene), ]
  merged_out <- file.path(out_dir, "fed_compare_all.csv")
  write.csv(m, merged_out, row.names = FALSE)
  cat(sprintf("\nmerged %d genes -> %s\n", nrow(m), merged_out))

  agree <- function(x, ref) {
    ok <- is.finite(x) & is.finite(ref)
    x <- x[ok]; ref <- ref[ok]; ss <- sum((ref - mean(ref))^2)
    r2 <- if (ss > 0) 1 - sum((x - ref)^2) / ss else NA_real_
    sprintf("R^2=%.5f  max|dp|=%.2e  max|d log10p|=%.3f",
            r2, max(abs(x - ref)), max(abs(-log10(pmax(x,1e-300)) + log10(pmax(ref,1e-300)))))
  }
  cat("\n=== Burden p (direct: secure/plain vs R::SKAT) ===\n")
  if ("burden_p_secure" %in% names(m)) cat("  secure vs R:", agree(m$burden_p_secure, m$burden_p_R), "\n")
  if ("burden_p_plain"  %in% names(m)) cat("  plain  vs R:", agree(m$burden_p_plain,  m$burden_p_R), "\n")

  if ("skat_p_secure" %in% names(m)) {
    cat("\n=== SKAT p (secure/plain vs R::SKAT davies & liu) ===\n")
    cat("  secure vs davies:", agree(m$skat_p_secure, m$skat_p_R_davies), "\n")
    cat("  secure vs liu   :", agree(m$skat_p_secure, m$skat_p_R_liu),    "\n")
    if ("skat_p_plain" %in% names(m)) {
      cat("  plain  vs davies:", agree(m$skat_p_plain, m$skat_p_R_davies), "\n")
      cat("  plain  vs liu   :", agree(m$skat_p_plain, m$skat_p_R_liu),    "\n")
    }
  } else {
    cat("\n(fed_results.csv has no skat_p columns -> run with FED_PROBES>0 for the SKAT p compare)\n")
  }
}
