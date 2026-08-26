args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript run_reference_skat.R <preprocessed-dir> <output.csv>")
}

# 1. Load the shared preprocessing-block reader.
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_argument[[1]])
workbench_dir <- dirname(normalizePath(script_path))
source(file.path(workbench_dir, "..", "analysis", "input.R"))
source(file.path(workbench_dir, "..", "analysis", "skat.R"))

input <- read_preprocessed_input(args[[1]])
phenotype_count <- ncol(input$phenotypes)
gene_metadata <- read.delim(
  file.path(args[[1]], "gene_metadata.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
phenotype_names <- readLines(file.path(args[[1]], "phenotypes.txt"))

# 2. Fit one pooled continuous null model per phenotype.
null_models <- fit_null_models(input)

# 3. Compute pooled Burden and SKAT-Davies for every gene and phenotype.
results <- vector("list", length(input$genes) * phenotype_count)
result_index <- 1
for (gene in seq_along(input$genes)) {
  genotype <- read_gene_genotype(input, gene)

  for (phenotype in seq_len(phenotype_count)) {
    null_model <- null_models[[phenotype]]
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

    results[[result_index]] <- data.frame(
      chromosome = gene_metadata$chromosome[[gene]],
      gene_index = gene - 1,
      gene_id = input$genes[[gene]],
      gene_symbol = gene_metadata$gene_symbol[[gene]],
      gene_order = gene_metadata$gene_order[[gene]],
      phenotype_index = phenotype - 1,
      phenotype_name = phenotype_names[[phenotype]],
      r_burden_p = burden$p.value,
      r_skat_davies_p = skat_davies$p.value,
      r_skat_davies_converged = skat_davies$param$Is_Converged
    )
    result_index <- result_index + 1
  }
}

# 4. Write gene-major CSV output for the Workbench join step.
output <- do.call(rbind, results)
for (column in c("r_burden_p", "r_skat_davies_p")) {
  output[[column]] <- sprintf("%.17g", output[[column]])
}
write.table(
  output,
  file = args[[2]],
  sep = ",",
  quote = TRUE,
  qmethod = "double",
  row.names = FALSE,
  na = "NA"
)
