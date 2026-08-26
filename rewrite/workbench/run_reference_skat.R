args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0(
    "usage: Rscript run_reference_skat.R ",
    "<preprocessed-dir> <output.csv> <timing.csv>"
  ))
}

timing_fields <- c(
  "component",
  "scope",
  "chromosome",
  "party",
  "lane",
  "phase",
  "parent_phase",
  "event_index",
  "batch_index",
  "batch_width",
  "batch_gene_count",
  "gene_index",
  "gene_id",
  "phenotype_index",
  "phenotype_name",
  "trace_mode",
  "measurement_kind",
  "elapsed_seconds",
  "sample_count_a",
  "sample_count_b",
  "public_variant_count",
  "private_variant_count",
  "ckks",
  "data_bits",
  "frac_bits",
  "probes",
  "status"
)

timing_rows <- list()
timing_chromosome <- ""
timing_sample_count_a <- ""
timing_sample_count_b <- ""

add_timing <- function(
    phase,
    scope,
    elapsed,
    status = "success",
    parent_phase = "",
    gene_index = "",
    gene_id = "",
    phenotype_index = "",
    phenotype_name = "",
    public_variant_count = "",
    private_variant_count = "") {
  row <- as.list(setNames(rep("", length(timing_fields)), timing_fields))
  row$component <- "r_reference"
  row$scope <- scope
  row$chromosome <- timing_chromosome
  row$phase <- phase
  row$parent_phase <- parent_phase
  row$event_index <- length(timing_rows)
  row$gene_index <- gene_index
  row$gene_id <- gene_id
  row$phenotype_index <- phenotype_index
  row$phenotype_name <- phenotype_name
  row$measurement_kind <- "actual"
  row$elapsed_seconds <- sprintf("%.9f", elapsed)
  row$sample_count_a <- timing_sample_count_a
  row$sample_count_b <- timing_sample_count_b
  row$public_variant_count <- public_variant_count
  row$private_variant_count <- private_variant_count
  row$status <- status
  timing_rows[[length(timing_rows) + 1]] <<- row
}

measure_value <- function(
    operation,
    phase,
    scope,
    parent_phase = "",
    gene_index = "",
    gene_id = "",
    phenotype_index = "",
    phenotype_name = "",
    public_variant_count = "",
    private_variant_count = "") {
  started <- proc.time()[["elapsed"]]
  status <- "success"
  tryCatch(
    force(operation),
    error = function(error) {
      status <<- "failure"
      stop(error)
    },
    finally = add_timing(
      phase = phase,
      scope = scope,
      elapsed = proc.time()[["elapsed"]] - started,
      status = status,
      parent_phase = parent_phase,
      gene_index = gene_index,
      gene_id = gene_id,
      phenotype_index = phenotype_index,
      phenotype_name = phenotype_name,
      public_variant_count = public_variant_count,
      private_variant_count = private_variant_count
    )
  )
}

write_timing <- function(path) {
  if (length(timing_rows) == 0) {
    output <- as.data.frame(
      setNames(rep(list(character()), length(timing_fields)), timing_fields)
    )
  } else {
    output <- do.call(
      rbind.data.frame,
      c(timing_rows, stringsAsFactors = FALSE)
    )
    output <- output[, timing_fields, drop = FALSE]
    output$chromosome[output$chromosome == ""] <- timing_chromosome
    output$sample_count_a[output$sample_count_a == ""] <- timing_sample_count_a
    output$sample_count_b[output$sample_count_b == ""] <- timing_sample_count_b
  }

  temporary <- paste0(path, ".tmp")
  write.table(
    output,
    file = temporary,
    sep = ",",
    quote = TRUE,
    qmethod = "double",
    row.names = FALSE,
    na = ""
  )
  if (!file.rename(temporary, path)) {
    stop("failed to install timing CSV: ", path)
  }
}

run_reference <- function() {
  run_started <- proc.time()[["elapsed"]]
  run_status <- "failure"
  on.exit({
    add_timing(
      phase = "r_total",
      scope = "chromosome",
      elapsed = proc.time()[["elapsed"]] - run_started,
      status = run_status
    )
    write_timing(args[[3]])
  }, add = TRUE)

  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  script_path <- sub("^--file=", "", script_argument[[1]])
  workbench_dir <- dirname(normalizePath(script_path))
  source(file.path(workbench_dir, "..", "analysis", "input.R"))
  source(file.path(workbench_dir, "..", "analysis", "skat.R"))

  input <- measure_value(
    read_preprocessed_input(args[[1]]),
    phase = "input_load",
    scope = "chromosome",
    parent_phase = "r_total"
  )
  phenotype_count <- ncol(input$phenotypes)
  gene_metadata <- read.delim(
    file.path(args[[1]], "gene_metadata.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  phenotype_names <- readLines(file.path(args[[1]], "phenotypes.txt"))
  timing_chromosome <<- paste0("chr", gene_metadata$chromosome[[1]])
  timing_sample_count_a <<- input$sample_count_a
  timing_sample_count_b <<- input$sample_count_b

  null_models <- vector("list", phenotype_count)
  for (phenotype in seq_len(phenotype_count)) {
    null_models[[phenotype]] <- measure_value(
      fit_null_model(input, phenotype),
      phase = "null_model_fit",
      scope = "phenotype",
      parent_phase = "r_total",
      phenotype_index = phenotype - 1,
      phenotype_name = phenotype_names[[phenotype]]
    )
  }

  results <- vector("list", length(input$genes) * phenotype_count)
  result_index <- 1
  for (gene in seq_along(input$genes)) {
    gene_started <- proc.time()[["elapsed"]]
    gene_status <- "failure"
    gene_id <- input$genes[[gene]]
    public_variant_count <- input$public_variant_counts[[gene]]
    private_variant_count <- ""

    tryCatch({
      genotype <- measure_value(
        read_gene_genotype(input, gene),
        phase = "genotype_load",
        scope = "gene",
        parent_phase = "r_gene_total",
        gene_index = gene - 1,
        gene_id = gene_id,
        public_variant_count = public_variant_count
      )
      private_variant_count <- ncol(genotype) - public_variant_count

      for (phenotype in seq_len(phenotype_count)) {
        null_model <- null_models[[phenotype]]
        phenotype_name <- phenotype_names[[phenotype]]
        burden <- measure_value(
          SKAT::SKAT(
            genotype,
            null_model,
            kernel = "linear.weighted",
            method = "davies",
            weights.beta = c(1, 25),
            r.corr = 1
          ),
          phase = "burden",
          scope = "gene_phenotype",
          parent_phase = "r_gene_total",
          gene_index = gene - 1,
          gene_id = gene_id,
          phenotype_index = phenotype - 1,
          phenotype_name = phenotype_name,
          public_variant_count = public_variant_count,
          private_variant_count = private_variant_count
        )
        skat_davies <- measure_value(
          SKAT::SKAT(
            genotype,
            null_model,
            kernel = "linear.weighted",
            method = "davies",
            weights.beta = c(1, 25),
            r.corr = 0
          ),
          phase = "skat_davies",
          scope = "gene_phenotype",
          parent_phase = "r_gene_total",
          gene_index = gene - 1,
          gene_id = gene_id,
          phenotype_index = phenotype - 1,
          phenotype_name = phenotype_name,
          public_variant_count = public_variant_count,
          private_variant_count = private_variant_count
        )

        results[[result_index]] <- data.frame(
          chromosome = gene_metadata$chromosome[[gene]],
          gene_index = gene - 1,
          gene_id = gene_id,
          gene_symbol = gene_metadata$gene_symbol[[gene]],
          gene_order = gene_metadata$gene_order[[gene]],
          phenotype_index = phenotype - 1,
          phenotype_name = phenotype_name,
          r_burden_p = burden$p.value,
          r_skat_davies_p = skat_davies$p.value,
          r_skat_davies_converged = skat_davies$param$Is_Converged
        )
        result_index <- result_index + 1
      }
      gene_status <- "success"
    }, finally = add_timing(
      phase = "r_gene_total",
      scope = "gene",
      elapsed = proc.time()[["elapsed"]] - gene_started,
      status = gene_status,
      parent_phase = "r_total",
      gene_index = gene - 1,
      gene_id = gene_id,
      public_variant_count = public_variant_count,
      private_variant_count = private_variant_count
    ))
  }

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
  run_status <- "success"
}

run_reference()
