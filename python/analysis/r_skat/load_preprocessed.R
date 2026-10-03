read_text_matrix <- function(path) {
  as.matrix(read.table(path, header = FALSE, check.names = FALSE))
}

read_dosage_matrix <- function(path, row_count, column_count = NULL) {
  # 1. Infer the number of columns when reading a private block.
  byte_count <- file.info(path)$size

  if (is.null(column_count)) {
    if (byte_count %% row_count != 0) {
      stop(sprintf(
        "%s has %d bytes, not divisible by %d rows",
        path, byte_count, row_count
      ))
    }
    column_count <- byte_count / row_count
  }

  # 2. Read the signed int8 dosages in preprocessing row order.
  values <- readBin(
    path,
    integer(),
    n = row_count * column_count,
    size = 1,
    signed = TRUE
  )

  # 3. Restore the row-major dosage matrix.
  matrix(
    values,
    nrow = row_count,
    ncol = column_count,
    byrow = TRUE
  )
}

read_preprocessed_input <- function(input_dir_a, input_dir_b = input_dir_a) {
  input_dir_a <- normalizePath(input_dir_a, mustWork = TRUE)

  input_dir_b <- normalizePath(input_dir_b, mustWork = TRUE)
  for (name in c("genes.txt", "block_sizes.txt", "public_variants.tsv")) {
    if (input_dir_a != input_dir_b && !identical(
      readLines(file.path(input_dir_a, name)), readLines(file.path(input_dir_b, name))
    )) stop(sprintf("A/B input mismatch: %s", name))
  }

  # 1. Read the ordered gene list and public variant counts.
  genes <- readLines(file.path(input_dir_a, "genes.txt"))
  public_variant_counts <- scan(
    file.path(input_dir_a, "block_sizes.txt"),
    what = integer(),
    quiet = TRUE
  )

  if (length(genes) != length(public_variant_counts)) {
    stop("genes.txt and block_sizes.txt have different lengths")
  }

  # 2. Read cohort-local covariates and phenotypes.
  covariates_a <- read_text_matrix(
    file.path(input_dir_a, "A", "cov.txt")
  )
  covariates_b <- read_text_matrix(
    file.path(input_dir_b, "B", "cov.txt")
  )
  phenotypes_a <- read_text_matrix(
    file.path(input_dir_a, "A", "pheno.txt")
  )
  phenotypes_b <- read_text_matrix(
    file.path(input_dir_b, "B", "pheno.txt")
  )

  # 3. Assemble the pooled X and Y in A-then-B sample order.
  list(
    input_dir_a = input_dir_a,
    input_dir_b = input_dir_b,
    genes = genes,
    public_variant_counts = public_variant_counts,
    sample_count_a = nrow(covariates_a),
    sample_count_b = nrow(covariates_b),
    covariates = rbind(covariates_a, covariates_b),
    phenotypes = rbind(phenotypes_a, phenotypes_b)
  )
}

read_gene_genotype <- function(input, gene) {
  block_name <- sprintf("block.%d.bin", gene - 1)
  public_variant_count <- input$public_variant_counts[[gene]]

  # 1. Read the public genotype block from A and B.
  public_a <- read_dosage_matrix(
    file.path(input$input_dir_a, "A", "geno", block_name),
    input$sample_count_a,
    public_variant_count
  )
  public_b <- read_dosage_matrix(
    file.path(input$input_dir_b, "B", "geno", block_name),
    input$sample_count_b,
    public_variant_count
  )

  # 2. Read the B-private genotype block and infer its width.
  private_b <- read_dosage_matrix(
    file.path(input$input_dir_b, "B", "private", block_name),
    input$sample_count_b
  )
  private_variant_count <- ncol(private_b)

  # 3. Build pooled G = [Gp, Gv] with zeros for A-private entries.
  rbind(
    cbind(
      public_a,
      matrix(0, input$sample_count_a, private_variant_count)
    ),
    cbind(
      public_b,
      private_b
    )
  )
}
