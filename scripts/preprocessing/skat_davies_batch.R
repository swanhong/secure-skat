#!/usr/bin/env Rscript

# Batched, fail-closed wrapper around the Davies/qfc routine bundled with R::SKAT.
# Low-tail Ruben and closed-form cases are handled in Python before this helper is called.

suppressPackageStartupMessages(library(SKAT))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: skat_davies_batch.R CASES.csv OUTPUT.csv TOL")
}

input_path <- args[[1]]
output_path <- args[[2]]
requested_tol <- as.numeric(args[[3]])
if (!is.finite(requested_tol) || requested_tol <= 0 || requested_tol >= 1) {
  stop("TOL must lie in (0,1)")
}

cases <- read.csv(input_path, stringsAsFactors = FALSE)
required <- c("case_id", "eigen_index", "q", "lambda")
if (!all(required %in% names(cases))) {
  stop("input CSV is missing required columns")
}

DAVIES_LIM <- 10000000L
BASE_ACC <- max(requested_tol, 1e-10)
DEEP_ACC <- max(min(BASE_ACC * 1e-3, 1e-13), .Machine$double.eps * 16)
davies_fn <- get0("SKAT_davies", envir = asNamespace("SKAT"), inherits = FALSE)
if (is.null(davies_fn)) {
  stop("installed SKAT package does not provide the required SKAT_davies backend")
}

run_davies <- function(q, lambda, acc) {
  fit <- davies_fn(q, lambda, acc = acc, lim = DAVIES_LIM)
  if (!is.list(fit) || !all(c("Qq", "ifault") %in% names(fit))) {
    stop("SKAT_davies returned an unsupported result shape")
  }
  fit
}

groups <- split(cases, cases$case_id)
ids <- sort(as.integer(names(groups)))
result <- vector("list", length(ids))
for (i in seq_along(ids)) {
  id <- ids[[i]]
  rows <- groups[[as.character(id)]]
  rows <- rows[order(rows$eigen_index), ]
  q <- unique(rows$q)
  lambda <- rows$lambda
  if (length(q) != 1L || !is.finite(q) || q <= 0 ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop(sprintf("invalid normalized mixture case %d", id))
  }

  fit <- run_davies(q, lambda, BASE_ACC)
  pvalue <- fit$Qq
  ifault <- fit$ifault
  acc_used <- BASE_ACC

  # qfc's requested accuracy is absolute. Re-evaluate small tails more tightly; a value is accepted
  # only when it is separated from zero by at least ten requested error units.
  needs_deep <- (
    ifault != 0L || length(pvalue) != 1L || !is.finite(pvalue) ||
    pvalue <= 10 * BASE_ACC || pvalue > 1
  )
  if (needs_deep) {
    strict_fit <- run_davies(q, lambda, DEEP_ACC)
    strict_p <- strict_fit$Qq
    if (strict_fit$ifault == 0L && length(strict_p) == 1L && is.finite(strict_p) &&
        strict_p > 10 * DEEP_ACC && strict_p <= 1) {
      fit <- strict_fit
      pvalue <- strict_p
      ifault <- strict_fit$ifault
      acc_used <- DEEP_ACC
    } else {
      pvalue <- strict_p
      ifault <- if (strict_fit$ifault == 0L) 9L else strict_fit$ifault
      acc_used <- DEEP_ACC
    }
  }

  result[[i]] <- data.frame(
    case_id = id,
    pvalue = pvalue,
    ifault = ifault,
    acc_used = acc_used,
    stringsAsFactors = FALSE
  )
}

options(digits = 17)
write.csv(do.call(rbind, result), output_path, row.names = FALSE, quote = FALSE)
