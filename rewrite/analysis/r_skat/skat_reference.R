fit_null_models <- function(input) {
  phenotype_count <- ncol(input$phenotypes)
  null_models <- vector("list", phenotype_count)

  for (phenotype in seq_len(phenotype_count)) {
    # 1. Build the pooled regression table for this phenotype.
    null_data <- data.frame(
      y = input$phenotypes[, phenotype],
      input$covariates,
      check.names = FALSE
    )
    names(null_data) <- c(
      "y",
      paste0("covariate", seq_len(ncol(input$covariates)))
    )

    # 2. Fit the continuous pooled null model with an intercept.
    null_models[[phenotype]] <- SKAT::SKAT_Null_Model(
      y ~ .,
      data = null_data,
      out_type = "C",
      Adjustment = FALSE
    )
  }

  null_models
}

compute_skat_pvalues <- function(genotype, null_models) {
  phenotype_count <- length(null_models)
  constant_result <- list(
    burden_p = 1,
    skat_davies_p = 1,
    skat_davies_converged = NA_integer_,
    skat_liu_p = 1
  )

  # Preserve empty genes in output order with their zero contribution.
  if (ncol(genotype) == 0) {
    return(rep(list(constant_result), phenotype_count))
  }

  reference_model <- null_models[[1]]
  for (null_model in null_models[-1]) {
    if (
      !identical(null_model$id_include, reference_model$id_include) ||
      !identical(null_model$X1, reference_model$X1)
    ) {
      stop("SKAT null models do not share rows and covariates")
    }
  }

  sample_count <- reference_model$n.all
  if (is.null(sample_count)) {
    sample_count <- nrow(genotype)
  }
  checked_genotype <- SKAT:::SKAT_MAIN_Check_Z(
    Z = genotype,
    n = sample_count,
    id_include = reference_model$id_include,
    SetID = NULL,
    weights = NULL,
    weights.beta = c(1, 25),
    impute.method = "fixed",
    is_check_genotype = TRUE,
    is_dosage = FALSE,
    missing_cutoff = 0.15,
    max_maf = 1,
    estimate_MAF = 1
  )
  if (checked_genotype$return == 1) {
    constant_result$burden_p <- checked_genotype$p.value
    constant_result$skat_davies_p <- checked_genotype$p.value
    constant_result$skat_liu_p <- checked_genotype$p.value
    return(rep(list(constant_result), phenotype_count))
  }

  # Reuse the weighted genotype and projected kernels across phenotypes.
  weighted_genotype <- t(
    t(checked_genotype$Z.test) * checked_genotype$weights
  )
  covariates <- reference_model$X1
  covariate_inverse <- solve(t(covariates) %*% covariates)

  genotype_covariates <- t(weighted_genotype) %*% covariates
  skat_kernel <-
    t(weighted_genotype) %*% weighted_genotype -
    genotype_covariates %*% covariate_inverse %*%
      t(genotype_covariates)

  burden_genotype <- cbind(rowSums(weighted_genotype))
  burden_covariates <- t(burden_genotype) %*% covariates
  burden_kernel <-
    t(burden_genotype) %*% burden_genotype -
    burden_covariates %*% covariate_inverse %*%
      t(burden_covariates)

  lapply(null_models, function(null_model) {
    skat_score <- t(null_model$res) %*% weighted_genotype
    skat_q <-
      skat_score %*% t(skat_score) / null_model$s2 / 2
    skat_davies <- SKAT:::Get_Davies_PVal(
      skat_q,
      skat_kernel
    )
    skat_liu <- SKAT:::Get_Liu_PVal(skat_q, skat_kernel)

    burden_score <- t(null_model$res) %*% burden_genotype
    burden_q <-
      burden_score %*% t(burden_score) / null_model$s2 / 2
    burden <- SKAT:::Get_Davies_PVal(
      burden_q,
      burden_kernel
    )

    list(
      burden_p = burden$p.value,
      skat_davies_p = skat_davies$p.value,
      skat_davies_converged =
        skat_davies$param$Is_Converged,
      skat_liu_p = skat_liu$p.value
    )
  })
}

run_skat_analysis <- function(input) {
  phenotype_count <- ncol(input$phenotypes)
  gene_count <- length(input$genes)
  progress_every <- as.integer(
    Sys.getenv("REFERENCE_PROGRESS_EVERY", "250")
  )
  if (is.na(progress_every) || progress_every < 0) {
    stop("REFERENCE_PROGRESS_EVERY must be a non-negative integer")
  }
  progress_label <- sprintf(
    "%s %s",
    basename(dirname(input$input_dir)),
    basename(input$input_dir)
  )
  started_at <- proc.time()[["elapsed"]]

  # 1. Fit one pooled null model per phenotype.
  null_models <- fit_null_models(input)

  # 2. Allocate the gene-major output rows.
  results <- vector(
    "list",
    gene_count * phenotype_count
  )
  result_index <- 1

  for (gene in seq_along(input$genes)) {
    # 3. Load the pooled public and B-private genotypes for this gene.
    genotype <- read_gene_genotype(input, gene)
    pvalues_by_phenotype <- compute_skat_pvalues(
      genotype,
      null_models
    )

    for (phenotype in seq_len(phenotype_count)) {
      # 4. Select the three p-values for this phenotype.
      pvalues <- pvalues_by_phenotype[[phenotype]]

      # 5. Store the result in g*q+t order.
      results[[result_index]] <- data.frame(
        gene_index = gene - 1,
        gene_id = input$genes[[gene]],
        phenotype_index = phenotype - 1,
        burden_p = pvalues$burden_p,
        skat_davies_p = pvalues$skat_davies_p,
        skat_davies_converged =
          pvalues$skat_davies_converged,
        skat_liu_p = pvalues$skat_liu_p
      )
      result_index <- result_index + 1
    }

    if (
      progress_every > 0 &&
      (gene %% progress_every == 0 || gene == gene_count)
    ) {
      elapsed_seconds <- proc.time()[["elapsed"]] - started_at
      remaining_seconds <- elapsed_seconds * (gene_count - gene) / gene
      message(sprintf(
        "[%s] %d/%d genes (%.1f%%), elapsed=%.1fm, ETA=%.1fm",
        progress_label,
        gene,
        gene_count,
        100 * gene / gene_count,
        elapsed_seconds / 60,
        remaining_seconds / 60
      ))
    }
  }

  do.call(rbind, results)
}
