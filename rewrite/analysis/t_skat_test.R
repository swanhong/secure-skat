source("rewrite/analysis/skat.R")

sample_count <- 20
null_data <- data.frame(
  y = sin(seq_len(sample_count)),
  covariate = seq_len(sample_count)
)
null_model <- SKAT::SKAT_Null_Model(
  y ~ covariate,
  data = null_data,
  out_type = "C",
  Adjustment = FALSE
)

genotype <- matrix(
  numeric(),
  nrow = sample_count,
  ncol = 0
)
pvalues <- compute_skat_pvalues(genotype, null_model)

stopifnot(
  pvalues$burden_p == 1,
  pvalues$skat_davies_p == 1,
  is.na(pvalues$skat_davies_converged),
  pvalues$skat_liu_p == 1
)
