args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript compute_reference.R <preprocessed-dir>")
}

# 1. Locate and load the analysis modules.
script_argument <- grep(
  "^--file=",
  commandArgs(),
  value = TRUE
)
script_path <- sub("^--file=", "", script_argument[[1]])
script_dir <- dirname(normalizePath(script_path))

source(file.path(script_dir, "load_preprocessed.R"))
source(file.path(script_dir, "skat_reference.R"))

# 2. Read the preprocessing output.
input <- read_preprocessed_input(args[[1]])

# 3. Compute Burden, SKAT-Davies, and SKAT-Liu.
results <- run_skat_analysis(input)

# 4. Write gene-major TSV output.
write.table(
  results,
  file = stdout(),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
