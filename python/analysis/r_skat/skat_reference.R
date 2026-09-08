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

  lapply(null_models, function(null_model) {
    skat_davies <- SKAT::SKAT(
      genotype,
      null_model,
      method = "davies",
      r.corr = 0,
      kernel = "linear.weighted",
      weights.beta = c(1, 25)
    )
    davies_converged <- skat_davies$param$Is_Converged
    if (length(davies_converged) == 0) {
      davies_converged <- NA_integer_
    }

    list(
      burden_p = SKAT::SKAT(
        genotype,
        null_model,
        method = "davies",
        r.corr = 1,
        kernel = "linear.weighted",
        weights.beta = c(1, 25)
      )$p.value,
      skat_davies_p = skat_davies$p.value,
      skat_davies_converged = davies_converged,
      skat_liu_p = SKAT::SKAT(
        genotype,
        null_model,
        method = "liu",
        r.corr = 0,
        kernel = "linear.weighted",
        weights.beta = c(1, 25)
      )$p.value
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
        skat_davies_converged = pvalues$skat_davies_converged,
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
