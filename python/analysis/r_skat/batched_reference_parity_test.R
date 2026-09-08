args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript batched_reference_parity_test.R <fixture-dir>")
}

script_argument <- grep(
  "^--file=",
  commandArgs(),
  value = TRUE
)
script_path <- sub("^--file=", "", script_argument[[1]])
script_dir <- dirname(normalizePath(script_path))

source(file.path(script_dir, "load_preprocessed.R"))
source(file.path(script_dir, "skat_reference.R"))

top_level_pvalues <- function(genotype, null_model) {
  if (ncol(genotype) == 0) {
    return(list(
      burden_p = 1,
      skat_davies_p = 1,
      skat_davies_converged = NA_integer_,
      skat_liu_p = 1
    ))
  }

  burden <- SKAT::SKAT(
    genotype,
    null_model,
    kernel = "linear.weighted",
    method = "davies",
    weights.beta = c(1, 25),
    r.corr = 1
  )
  skat_davies <- SKAT::SKAT(
    genotype,
    null_model,
    kernel = "linear.weighted",
    method = "davies",
    weights.beta = c(1, 25),
    r.corr = 0
  )
  skat_liu <- SKAT::SKAT(
    genotype,
    null_model,
    kernel = "linear.weighted",
    method = "liu",
    weights.beta = c(1, 25),
    r.corr = 0
  )

  davies_converged <- skat_davies$param$Is_Converged
  if (length(davies_converged) == 0) {
    davies_converged <- NA_integer_
  }

  list(
    burden_p = burden$p.value,
    skat_davies_p = skat_davies$p.value,
    skat_davies_converged = davies_converged,
    skat_liu_p = skat_liu$p.value
  )
}

assert_pvalues_equal <- function(actual, expected, label) {
  for (field in c("burden_p", "skat_davies_p", "skat_liu_p")) {
    actual_logp <- -log10(actual[[field]])
    expected_logp <- -log10(expected[[field]])
    if (!isTRUE(all.equal(
      actual_logp,
      expected_logp,
      tolerance = 1e-12
    ))) {
      stop(sprintf("%s field %s differs on -log10(p)", label, field))
    }
  }
  if (!identical(
    actual$skat_davies_converged,
    expected$skat_davies_converged
  )) {
    stop(sprintf("%s Davies convergence differs", label))
  }
}

input <- read_preprocessed_input(args[[1]])
null_models <- fit_null_models(input)
genotypes <- lapply(
  seq_along(input$genes),
  function(gene) read_gene_genotype(input, gene)
)
genotypes[[length(genotypes) + 1]] <- matrix(
  0L,
  nrow = nrow(input$covariates),
  ncol = 2
)

for (gene in seq_along(genotypes)) {
  genotype <- genotypes[[gene]]
  actual <- suppressWarnings(
    compute_skat_pvalues(genotype, null_models)
  )

  for (phenotype in seq_along(null_models)) {
    expected <- suppressWarnings(
      top_level_pvalues(
        genotype,
        null_models[[phenotype]]
      )
    )
    assert_pvalues_equal(
      actual[[phenotype]],
      expected,
      sprintf("gene %d phenotype %d", gene - 1, phenotype - 1)
    )
  }
}

# Exercise the non-converged Davies fallback near p = 1e-50.
set.seed(20260831)
sample_count <- 20000
variant_count <- 20
deep_tail_genotype <- matrix(
  rbinom(sample_count * variant_count, 2, 0.01),
  nrow = sample_count
)
deep_tail_data <- data.frame(
  covariate1 = rnorm(sample_count),
  covariate2 = rnorm(sample_count)
)
genetic_score <- as.numeric(
  scale(rowSums(deep_tail_genotype))
)
deep_tail_data$y <-
  0.125 * genetic_score + rnorm(sample_count)
deep_tail_model <- SKAT::SKAT_Null_Model(
  y ~ .,
  data = deep_tail_data,
  out_type = "C",
  Adjustment = FALSE
)
deep_tail_expected <- top_level_pvalues(
  deep_tail_genotype,
  deep_tail_model
)
deep_tail_actual <- compute_skat_pvalues(
  deep_tail_genotype,
  list(deep_tail_model)
)[[1]]
deep_tail_logp <- -log10(deep_tail_expected$skat_davies_p)

if (deep_tail_logp < 45 || deep_tail_logp > 55) {
  stop(sprintf("expected deep-tail log-p near 50, got %.6f", deep_tail_logp))
}
if (deep_tail_expected$skat_davies_converged != 0) {
  stop("expected deep-tail Davies fallback")
}
assert_pvalues_equal(
  deep_tail_actual,
  deep_tail_expected,
  "deep tail"
)
