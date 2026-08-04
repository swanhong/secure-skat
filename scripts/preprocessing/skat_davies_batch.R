#!/usr/bin/env Rscript

# Batched Davies and Liu p-values from the routines bundled with R::SKAT.

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
BASE_ACC <- max(requested_tol, 1e-6)
TAIL_ACC <- max(min(requested_tol, 1e-10), .Machine$double.eps * 16)
davies_fn <- get0("SKAT_davies", envir = asNamespace("SKAT"), inherits = FALSE)
liu_fn <- get0("SKAT_liu", envir = asNamespace("SKAT"), inherits = FALSE)
if (is.null(davies_fn) || is.null(liu_fn)) {
  stop("installed SKAT package does not provide the required p-value backends")
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
  if (length(q) == 1L && q == 0 && length(lambda) == 1L && lambda == 0) {
    result[[i]] <- data.frame(
      case_id = id, davies = 1, liu = 1, davies_ifault = 0L,
      acc_used = BASE_ACC, stringsAsFactors = FALSE
    )
    next
  }
  if (length(q) != 1L || !is.finite(q) || q <= 0 ||
      any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop(sprintf("invalid normalized mixture case %d", id))
  }

  fit <- run_davies(q, lambda, BASE_ACC)
  pvalue <- fit$Qq
  ifault <- fit$ifault
  acc_used <- BASE_ACC

  # Match R::SKAT's regular accuracy, and only ask for more precision in a resolved small tail.
  if (ifault == 0L && length(pvalue) == 1L && is.finite(pvalue) &&
      pvalue > 0 && pvalue <= 10 * BASE_ACC) {
    tail_fit <- run_davies(q, lambda, TAIL_ACC)
    tail_p <- tail_fit$Qq
    if (tail_fit$ifault == 0L && length(tail_p) == 1L && is.finite(tail_p) &&
        tail_p > 10 * TAIL_ACC && tail_p <= 1) {
      pvalue <- tail_p
      ifault <- tail_fit$ifault
      acc_used <- TAIL_ACC
    } else {
      pvalue <- tail_p
      ifault <- if (tail_fit$ifault == 0L) 9L else tail_fit$ifault
      acc_used <- TAIL_ACC
    }
  }

  result[[i]] <- data.frame(
    case_id = id,
    davies = pvalue,
    liu = liu_fn(q, lambda),
    davies_ifault = ifault,
    acc_used = acc_used,
    stringsAsFactors = FALSE
  )
}

options(digits = 17)
write.csv(do.call(rbind, result), output_path, row.names = FALSE, quote = FALSE)
