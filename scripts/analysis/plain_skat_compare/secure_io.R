# Readers for secure pipeline outputs and scale-normalization helpers.

# Read a secure output file that stores a numeric matrix as comma-delimited rows.
read_secure_matrix <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) {
    return(NULL)
  }

  do.call(
    rbind,
    lapply(lines, function(line) {
      vals <- unlist(strsplit(line, ",", fixed = TRUE))
      vals <- trimws(vals)
      vals <- vals[nzchar(vals)]
      as.numeric(vals)
    })
  )
}

# Flatten a secure matrix file into a numeric vector in the script's expected order.
read_secure_vector <- function(path) {
  mat <- read_secure_matrix(path)
  if (is.null(mat)) {
    return(NULL)
  }

  as.numeric(t(mat))
}

# Try multiple candidate secure paths and return the first vector that exists.
read_secure_vector_any <- function(paths) {
  for (path in paths) {
    vec <- read_secure_vector(path)
    if (!is.null(vec)) {
      return(vec)
    }
  }
  NULL
}

# Read a simple comma-delimited numeric file into one vector.
read_comma_numeric_vector <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }

  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) {
    return(NULL)
  }

  vals <- unlist(strsplit(lines, ",", fixed = TRUE))
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  as.numeric(vals)
}

# Read one party's secure residual output (`ynew.txt`) when it is available.
read_secure_ynew <- function(party_dir, party_dirs, secure_run_root) {
  party_idx <- match(basename(party_dir), basename(party_dirs))
  ynew_path <- file.path(secure_party_dir(secure_run_root, party_idx), "ynew.txt")
  read_comma_numeric_vector(ynew_path)
}

# Pick the secure scale when present, otherwise fall back to the plain reference scale.
resolve_selected_secure_scale <- function(secure_scale_global, default_scale) {
  if (is.na(secure_scale_global)) {
    default_scale
  } else {
    secure_scale_global
  }
}

# Load the top-level secure SKAT and burden statistics.
# These files already contain the final secure outputs for direct comparison.
load_secure_scalar_results <- function(context) {
  secure_q_path <- file.path(secure_party_dir(context$secure_run_root, 1L), "skat_out.txt")
  secure_q <- if (file.exists(secure_q_path)) {
    as.numeric(scan(secure_q_path, quiet = TRUE, nmax = 1))
  } else {
    NA_real_
  }

  secure_scale_global <- NA_real_

  secure_burden_path <- file.path(secure_party_dir(context$secure_run_root, 1L), "burden_out.txt")
  secure_burden <- if (file.exists(secure_burden_path)) {
    as.numeric(scan(secure_burden_path, quiet = TRUE, nmax = 1))
  } else {
    NA_real_
  }

  list(
    secure_q_path = secure_q_path,
    secure_q = secure_q,
    secure_scale_global = secure_scale_global,
    secure_q_for_compare = secure_q,
    secure_burden_path = secure_burden_path,
    secure_burden = secure_burden,
    secure_burden_for_compare = secure_burden
  )
}

# Sum secure per-block SKAT and burden outputs for the selected block subset.
read_selected_secure_block_sums <- function(context) {
  selected_secure_q_raw <- NA_real_
  selected_secure_burden_sum_raw <- 0.0
  available <- TRUE

  for (block_idx in context$blocks_to_print) {
    secure_block_idx <- block_idx - 1L
    q_block_vec <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("qBlock_block%d.txt", secure_block_idx)
      )),
      1L
    )
    q_burden_block_vec <- trim_or_null(
      read_secure_vector(file.path(
        secure_party_dir(context$secure_run_root, 1L),
        sprintf("qBurdenBlock_block%d.txt", secure_block_idx)
      )),
      1L
    )

    if (is.null(q_block_vec) || is.null(q_burden_block_vec)) {
      available <- FALSE
      break
    }

    if (is.na(selected_secure_q_raw)) {
      selected_secure_q_raw <- q_block_vec[[1]]
    } else {
      selected_secure_q_raw <- selected_secure_q_raw + q_block_vec[[1]]
    }
    selected_secure_burden_sum_raw <- selected_secure_burden_sum_raw + q_burden_block_vec[[1]]
  }

  list(
    available = available,
    selected_secure_q_raw = selected_secure_q_raw,
    selected_secure_burden_sum_raw = selected_secure_burden_sum_raw
  )
}
