#!/usr/bin/env Rscript
# Plain R::SKAT reference + (when secure outputs exist) per-block compare & PNG.
# Handles two on-disk geno formats: "pgen" (plink2 --export A) and "blocks" (raw
# int8 .bin, row-major buf[i*bs + j], one file per block).
#
# NOTE: R::SKAT's Q is sigma^2-normalized; the v2 secure Q is NOT yet normalized,
# so the two differ by a (per-null-model) constant. Expected at v2 -> it goes away
# once the secure pipeline normalizes (v6+).
#
# Usage: Rscript skat_compare.R <repo_root> <data_dir> <pgen|blocks> [out_dir] [secure_dir]

suppressMessages(library(SKAT))

a <- commandArgs(trailingOnly = TRUE)
repo_root   <- normalizePath(a[[1]], mustWork = TRUE)
data_dir    <- a[[2]]
geno_format <- a[[3]]
out_dir     <- if (length(a) >= 4 && nzchar(a[[4]])) a[[4]] else file.path(repo_root, "out", "plain")
secure_dir  <- if (length(a) >= 5 && nzchar(a[[5]])) a[[5]] else ""
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

party_dirs  <- file.path(repo_root, data_dir, c("party1", "party2"))
block_sizes <- as.integer(readLines(file.path(party_dirs[1], "chrom_sizes.txt")))
nblocks     <- length(block_sizes)

read_pheno <- function(d) as.numeric(read.table(file.path(d, "pheno.txt"))[, 1])
read_cov   <- function(d) as.matrix(read.table(file.path(d, "cov.txt"), sep = "\t"))

# --- per-party geno block loaders -------------------------------------------
load_pgen <- local({
  plink2 <- Sys.which("plink2"); if (plink2 == "") stop("plink2 not in PATH")
  cache <- file.path(repo_root, ".local", "tmp", "plain_skat")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  function(party_dir, b) {
    pre <- file.path(cache, sprintf("%s_chr%02d", basename(party_dir), b + 1))
    raw <- sprintf("%s.raw", pre)
    if (!file.exists(raw)) {
      cmd <- sprintf('"%s" --pfile "%s" --keep "%s" --export A --out "%s"',
                     plink2, file.path(party_dir, "geno", sprintf("chr%d", b + 1)),
                     file.path(party_dir, "sample_keep.txt"), pre)
      if (system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) != 0)
        stop(sprintf("plink2 export failed: %s chr%d", basename(party_dir), b + 1))
    }
    df <- read.table(raw, header = TRUE, check.names = FALSE)
    as.matrix(df[, -(1:6), drop = FALSE])
  }
})
load_blocks <- function(party_dir, b, n) {
  bs <- block_sizes[b + 1]
  raw <- readBin(file.path(party_dir, sprintf("geno.%d.bin", b)), "integer",
                 n = n * bs, size = 1, signed = TRUE)
  matrix(raw, nrow = n, ncol = bs, byrow = TRUE)  # row-major buf[i*bs + j]
}

# --- null model (pooled pheno + covariates; SKAT adds its own intercept) -----
y <- unlist(lapply(party_dirs, read_pheno), use.names = FALSE)
X <- do.call(rbind, lapply(party_dirs, read_cov))
X <- X[, apply(X, 2, function(col) length(unique(col)) > 1), drop = FALSE]
obj <- SKAT_Null_Model(y ~ X, out_type = "C")
n_party <- vapply(party_dirs, function(d) length(read_pheno(d)), integer(1))

# --- per-block R::SKAT -------------------------------------------------------
q_skat <- numeric(nblocks); q_burden <- numeric(nblocks)
for (b in seq_len(nblocks) - 1L) {
  Z <- if (geno_format == "pgen")
    do.call(rbind, lapply(party_dirs, load_pgen, b = b))
  else
    do.call(rbind, Map(function(d, n) load_blocks(d, b, n), party_dirs, n_party))
  Z[is.na(Z)] <- 0
  q_skat[b + 1]   <- SKAT(Z, obj, weights.beta = c(1, 25))$Q
  q_burden[b + 1] <- SKAT(Z, obj, weights.beta = c(1, 25), r.corr = 1)$Q
  writeLines(format(q_skat[b + 1], digits = 10), file.path(out_dir, sprintf("plain_qBlock_block%d.txt", b)))
  cat(sprintf("block %d/%d  Q_skat=%.6e  Q_burden=%.6e\n", b, nblocks - 1, q_skat[b + 1], q_burden[b + 1]))
}
writeLines(format(sum(q_skat), digits = 10),   file.path(out_dir, "plain_skat_out.txt"))
writeLines(format(sum(q_burden), digits = 10), file.path(out_dir, "plain_burden_out.txt"))
cat(sprintf("Plain R::SKAT totals: SKAT=%.6e Burden=%.6e\n", sum(q_skat), sum(q_burden)))

# --- compare to secure per-block (if present) --------------------------------
read1 <- function(p) if (file.exists(p)) as.numeric(readLines(p)[1]) else NA_real_
secure <- if (nzchar(secure_dir)) {
  vapply(seq_len(nblocks) - 1L,
         function(b) read1(file.path(secure_dir, sprintf("qBlock_block%d.txt", b))), numeric(1))
} else {
  rep(NA_real_, nblocks)
}

if (all(is.na(secure))) {
  cat("(no secure per-block files -> compare skipped; run secure with SFGWAS_DEBUG=1)\n")
} else {
  ratio <- secure / q_skat
  print(data.frame(block = seq_len(nblocks) - 1L, secure = secure, plain = q_skat,
                   ratio = ratio, rel_err = abs(secure - q_skat) / pmax(abs(q_skat), 1e-12)), digits = 6)
  cat(sprintf("secure/plain ratio: median=%.6g (flat => only the missing sigma^2 normalization)\n",
              median(ratio, na.rm = TRUE)))
  png_out <- file.path(out_dir, "skat_block_compare.png")
  png(png_out, width = 900, height = 450); par(mfrow = c(1, 2))
  plot(q_skat, secure, pch = 19, xlab = "plain R::SKAT Q", ylab = "secure Q", main = "per-block Q: secure vs plain")
  abline(0, median(ratio, na.rm = TRUE), col = "red", lty = 2)
  plot(seq_len(nblocks) - 1L, ratio, type = "b", pch = 19, xlab = "block", ylab = "secure/plain", main = "per-block ratio")
  invisible(dev.off())
  cat(sprintf("PNG: %s\n", png_out))
}
