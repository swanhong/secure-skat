# Burden-specific summaries are kept separate from SKAT Q summaries so the
# main workflow can call them independently.

# Convert a burden sum into the SKAT-style burden Q statistic with the chosen scale.
compute_burden_q <- function(burden_sum, scale_factor) {
  (burden_sum^2) * scale_factor
}

# Summarize full-run and selected-block burden values from the plain reconstruction.
summarize_plain_burden_totals <- function(plain_result, scale_factor, blocks_to_print) {
  selected_plain_burden_sum <- sum(vapply(blocks_to_print, function(block_idx) {
    plain_result$plain_blocks[[block_idx]]$burden_block
  }, numeric(1)))

  list(
    burden_total_secure_style = plain_result$burden_total_secure_style_raw,
    burden_q_total_secure_style = compute_burden_q(plain_result$burden_total_secure_style_raw, scale_factor),
    selected_plain_burden_sum = selected_plain_burden_sum,
    selected_plain_burden_q = compute_burden_q(selected_plain_burden_sum, scale_factor)
  )
}

# Scale selected secure burden block sums into the final comparable burden statistic.
summarize_selected_secure_burden <- function(selected_secure_blocks, secure_scale) {
  if (!isTRUE(selected_secure_blocks$available)) {
    return(list(selected_secure_burden_q = NA_real_))
  }

  list(
    selected_secure_burden_q = compute_burden_q(
      selected_secure_blocks$selected_secure_burden_sum_raw,
      secure_scale
    )
  )
}
