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

compute_skat_pvalues <- function(genotype, null_model) {
  # Preserve empty genes in output order with their zero contribution.
  if (ncol(genotype) == 0) {
    return(list(
      burden_p = 1,
      skat_davies_p = 1,
      skat_davies_converged = NA_integer_,
      skat_liu_p = 1
    ))
  }

  # 1. Compute the Burden p-value using r.corr = 1.
  burden <- SKAT::SKAT(
    genotype,
    null_model,
    kernel = "linear.weighted",
    method = "davies",
    weights.beta = c(1, 25),
    r.corr = 1
  )

  # 2. Compute the SKAT p-value using Davies' method.
  skat_davies <- SKAT::SKAT(
    genotype,
    null_model,
    kernel = "linear.weighted",
    method = "davies",
    weights.beta = c(1, 25),
    r.corr = 0
  )

  # 3. Compute the SKAT p-value using Liu's method.
  skat_liu <- SKAT::SKAT(
    genotype,
    null_model,
    kernel = "linear.weighted",
    method = "liu",
    weights.beta = c(1, 25),
    r.corr = 0
  )

  # 4. Preserve missing convergence for degenerate genes
  davies_converged <- skat_davies$param$Is_Converged
  if (length(davies_converged) == 0) {
    davies_converged <- NA_integer_
  }

  # 5. Return only the comparison outputs.
  list(
    burden_p = burden$p.value,
    skat_davies_p = skat_davies$p.value,
    skat_davies_converged = davies_converged,
    skat_liu_p = skat_liu$p.value
  )
}

run_skat_analysis <- function(input) {
  phenotype_count <- ncol(input$phenotypes)

  # 1. Fit one pooled null model per phenotype.
  null_models <- fit_null_models(input)

  # 2. Allocate the gene-major output rows.
  results <- vector(
    "list",
    length(input$genes) * phenotype_count
  )
  result_index <- 1

  for (gene in seq_along(input$genes)) {
    # 3. Load the pooled public and B-private genotypes for this gene.
    genotype <- read_gene_genotype(input, gene)

    for (phenotype in seq_len(phenotype_count)) {
      # 4. Compute the three p-values for this gene and phenotype.
      pvalues <- compute_skat_pvalues(
        genotype,
        null_models[[phenotype]]
      )

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
  }

  do.call(rbind, results)
}
