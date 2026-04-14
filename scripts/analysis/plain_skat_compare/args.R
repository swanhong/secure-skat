# Argument parsing and path resolution for the plain-vs-secure SKAT comparison.

# Pop a `--flag value` pair from argv and return both the value and remaining args.
pop_flag_value <- function(args, flag, cast = identity) {
  flag_idx <- match(flag, args)
  if (is.na(flag_idx)) {
    return(list(args = args, value = NULL))
  }
  if (flag_idx == length(args)) {
    stop(sprintf("Missing value for %s", flag))
  }

  list(
    args = args[-c(flag_idx, flag_idx + 1L)],
    value = cast(args[[flag_idx + 1L]])
  )
}

# Parse CLI flags into a normalized config list before any filesystem access.
parse_compare_args <- function(args) {
  debug_mode <- "--debug" %in% args
  args <- args[args != "--debug"]
  skip_skat_package <- "--skip-skat-package" %in% args
  args <- args[args != "--skip-skat-package"]

  dataset_root <- "example_data"
  dataset_res <- pop_flag_value(args, "--dataset", identity)
  args <- dataset_res$args
  if (!is.null(dataset_res$value)) {
    dataset_root <- dataset_res$value
  }

  window_bp_res <- pop_flag_value(args, "--window-bp", as.integer)
  args <- window_bp_res$args
  window_bp <- if (is.null(window_bp_res$value)) NA_integer_ else window_bp_res$value

  step_bp_res <- pop_flag_value(args, "--step-bp", as.integer)
  args <- step_bp_res$args
  step_bp <- if (is.null(step_bp_res$value)) NA_integer_ else step_bp_res$value

  min_window_variants_res <- pop_flag_value(args, "--min-window-variants", as.integer)
  args <- min_window_variants_res$args
  min_window_variants <- if (is.null(min_window_variants_res$value)) 1L else min_window_variants_res$value

  window_limit_res <- pop_flag_value(args, "--window-limit", as.integer)
  args <- window_limit_res$args
  window_limit <- if (is.null(window_limit_res$value)) NA_integer_ else window_limit_res$value

  window_output_tag_res <- pop_flag_value(args, "--window-output-tag", identity)
  args <- window_output_tag_res$args
  window_output_tag <- window_output_tag_res$value

  if (!is.na(window_bp) && window_bp <= 0L) {
    stop("--window-bp must be a positive integer")
  }
  if (is.na(step_bp) && !is.na(window_bp)) {
    step_bp <- window_bp
  }
  if (!is.na(window_bp) && (is.na(step_bp) || step_bp <= 0L)) {
    stop("--step-bp must be a positive integer")
  }
  if (is.na(min_window_variants) || min_window_variants <= 0L) {
    stop("--min-window-variants must be a positive integer")
  }
  if (!is.na(window_limit) && window_limit <= 0L) {
    stop("--window-limit must be a positive integer")
  }

  list(
    debug_mode = debug_mode,
    skip_skat_package = skip_skat_package,
    dataset_root = dataset_root,
    window_bp = window_bp,
    step_bp = step_bp,
    min_window_variants = min_window_variants,
    window_limit = window_limit,
    window_output_tag = window_output_tag,
    repo_root = if (length(args) >= 1) args[[1]] else ".",
    run_id = if (length(args) >= 2) args[[2]] else "",
    blocks_arg = if (length(args) >= 3) args[[3]] else NULL
  )
}

# Resolve the secure run directory, preferring the newest match for a given run id.
resolve_secure_run_root <- function(repo_root, run_id) {
  if (!nzchar(run_id)) {
    return(file.path(repo_root, "out"))
  }

  candidates <- Sys.glob(file.path(repo_root, "out", sprintf("output_*_%s", run_id)))
  candidates <- candidates[dir.exists(candidates)]

  if (length(candidates) == 0) {
    stop(sprintf("No secure output directory found for run id: %s", run_id))
  }

  if (length(candidates) > 1) {
    info <- file.info(candidates)
    candidates <- candidates[order(info$mtime, decreasing = TRUE)]
  }

  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

# Interpret the optional block list argument and clamp it to valid block indices.
resolve_blocks_to_print <- function(blocks_arg, n_blocks) {
  blocks_to_print <- if (!is.null(blocks_arg)) {
    as.integer(strsplit(blocks_arg, ",", fixed = TRUE)[[1]])
  } else {
    c(1L, n_blocks)
  }

  blocks_to_print <- sort(unique(
    blocks_to_print[
      !is.na(blocks_to_print) &
      blocks_to_print >= 1L &
      blocks_to_print <= n_blocks
    ]
  ))

  if (length(blocks_to_print) == 0) {
    blocks_to_print <- c(1L, n_blocks)
  }

  blocks_to_print
}

# Build the path to one party's secure-output directory.
secure_party_dir <- function(secure_run_root, party_idx) {
  file.path(secure_run_root, sprintf("party%d", party_idx))
}

# Expand repo-relative paths, validate required inputs, and prepare cache locations.
resolve_compare_paths <- function(config) {
  repo_root <- normalizePath(config$repo_root, winslash = "/", mustWork = TRUE)
  dataset_root <- if (grepl("^/", config$dataset_root)) {
    normalizePath(config$dataset_root, winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(file.path(repo_root, config$dataset_root), winslash = "/", mustWork = TRUE)
  }

  party_dirs <- c(
    file.path(dataset_root, "party1"),
    file.path(dataset_root, "party2")
  )
  for (dir_path in party_dirs) {
    if (!dir.exists(dir_path)) {
      stop(sprintf("Missing data directory: %s", dir_path))
    }
  }

  chrom_sizes <- as.integer(readLines(file.path(party_dirs[[1]], "chrom_sizes.txt"), warn = FALSE))
  chrom_sizes <- chrom_sizes[!is.na(chrom_sizes)]
  n_blocks <- length(chrom_sizes)
  if (n_blocks == 0) {
    stop("No blocks found in chrom_sizes.txt")
  }

  blocks_to_print <- resolve_blocks_to_print(config$blocks_arg, n_blocks)

  plink2 <- Sys.which("plink2")
  if (plink2 == "") {
    stop("plink2 not found in PATH")
  }

  dataset_cache_tag <- gsub("[^A-Za-z0-9._-]+", "_", dataset_root)
  cache_dir <- file.path(repo_root, ".local", "tmp", "plain_skat_compare", dataset_cache_tag)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  csv_dir <- file.path(cache_dir, "variant_debug_csv")
  if (config$debug_mode) {
    dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  }

  c(
    config,
    list(
      repo_root = repo_root,
      dataset_root = dataset_root,
      party_dirs = party_dirs,
      chrom_sizes = chrom_sizes,
      n_blocks = n_blocks,
      blocks_to_print = blocks_to_print,
      plink2 = plink2,
      cache_dir = cache_dir,
      csv_dir = csv_dir,
      secure_run_root = resolve_secure_run_root(repo_root, config$run_id)
    )
  )
}
