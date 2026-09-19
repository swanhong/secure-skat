# Invoked by run_r.py; reads existing prepared blocks without modifying them.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L)
  stop("usage: run_genes.R loader.R prepared genes.tsv phenotypes.tsv chromosome output.csv")
if (!requireNamespace("SKAT", quietly = TRUE)) stop("Install the R package SKAT first")
source(args[[1]])
options(digits = 17)

input <- read_preprocessed_input(args[[2]])
genes <- read.delim(args[[3]], colClasses = "character", check.names = FALSE)
phenotypes <- read.delim(args[[4]], colClasses = "character", check.names = FALSE)
indices <- as.integer(phenotypes$index)
if (anyNA(indices) || any(indices < 1L | indices > ncol(input$phenotypes)))
  stop("Phenotype indices are outside the prepared matrix")
if (nrow(input$covariates) != nrow(input$phenotypes) ||
    any(!is.finite(input$covariates)) || any(!is.finite(input$phenotypes)))
  stop("Invalid or mismatched covariate/phenotype matrices")

# Fit once per phenotype within this prepared directory, then reuse for all genes.
null_models <- lapply(indices, function(index) {
  data <- data.frame(y = input$phenotypes[, index], input$covariates)
  names(data) <- c("y", paste0("cov", seq_len(ncol(input$covariates))))
  SKAT::SKAT_Null_Model(y ~ ., data = data, out_type = "C", Adjustment = FALSE)
})

run_test <- function(genotype, null_model, rho) {
  if (ncol(genotype) == 0L)
    return(list(p = NA_real_, flag = NA_integer_, markers = 0L,
                status = "not_tested", message = "No prepared variants"))
  messages <- character()
  result <- tryCatch(withCallingHandlers(
    SKAT::SKAT(genotype, null_model, method = "davies", r.corr = rho,
               kernel = "linear.weighted", weights.beta = c(1, 25)),
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  ), error = identity)
  if (inherits(result, "error"))
    return(list(p = NA_real_, flag = NA_integer_, markers = NA_integer_,
                status = "error", message = conditionMessage(result)))
  flag <- result$param$Is_Converged
  markers <- result$param$n.marker.test
  if (length(flag) != 1L) flag <- NA_integer_
  if (length(markers) != 1L) markers <- NA_integer_
  p <- result$p.value
  status <- "ok"
  if (!is.na(markers) && markers == 0L) {
    status <- "not_tested"
    p <- NA_real_
  } else if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
    status <- "invalid_p"
    p <- NA_real_
  } else if (!is.na(flag) && flag == 0L) {
    status <- "nonconverged"
  } else if (is.na(flag) || length(messages)) {
    status <- "check_diagnostics"
  }
  list(p = p, flag = flag, markers = markers, status = status,
       message = paste(unique(messages), collapse = " | "))
}

rows <- list()
for (gene in seq_len(nrow(genes))) {
  index <- as.integer(genes$index[[gene]])
  if (is.na(index) || index < 1L || index > length(input$genes) ||
      input$genes[[index]] != genes$prepared_gene_id[[gene]]) {
    stop("Gene index does not match genes.txt")
  }
  genotype <- read_gene_genotype(input, index)
  if (nrow(genotype) != nrow(input$phenotypes)) stop("Genotype row count mismatch")
  for (pheno in seq_len(nrow(phenotypes))) {
    burden <- run_test(genotype, null_models[[pheno]], 1)
    skat <- run_test(genotype, null_models[[pheno]], 0)
    rows[[length(rows) + 1L]] <- data.frame(
      gene_id = genes$gene_id[[gene]], gene_symbol = genes$gene_symbol[[gene]],
      chromosome = args[[5]], phenotype = phenotypes$name[[pheno]],
      phenotype_column = phenotypes$column[[pheno]],
      axa_phenotype_id = phenotypes$axa_id[[pheno]],
      R_burden = burden$p, R_skat = skat$p,
      R_burden_status = burden$status, R_skat_status = skat$status,
      R_skat_converged = skat$flag,
      n_samples = nrow(genotype), n_variants_input = ncol(genotype),
      n_variants_burden = burden$markers, n_variants_skat = skat$markers,
      R_burden_message = burden$message, R_skat_message = skat$message,
      R_version = R.version.string, SKAT_version = as.character(packageVersion("SKAT"))
    )
  }
  message(sprintf("chr%s %s (%d/%d genes): five phenotypes finished",
                  args[[5]], genes$gene_id[[gene]], gene, nrow(genes)))
}
write.csv(do.call(rbind, rows), file = args[[6]], row.names = FALSE, na = "NA")
