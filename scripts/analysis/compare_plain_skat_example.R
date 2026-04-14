# Thin entrypoint: load the responsibility-specific modules and hand off the
# entire workflow to run_plain_skat_compare().
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(file.path(getwd(), "scripts", "analysis"), winslash = "/", mustWork = TRUE)
}

module_dir <- file.path(script_dir, "plain_skat_compare")
module_files <- c(
  "args.R",
  "data_io.R",
  "burden.R",
  "secure_io.R",
  "reporting.R",
  "skat.R",
  "windows.R",
  "blocks.R",
  "workflow.R"
)

for (module_name in module_files) {
  source(file.path(module_dir, module_name), local = TRUE)
}

# Execute the modular plain-vs-secure SKAT workflow using the current CLI arguments.
main <- function() {
  run_plain_skat_compare(commandArgs(trailingOnly = TRUE))
}

if (sys.nframe() == 0) {
  main()
}
