#' Extract flagged trials
#'
#' Returns a data frame of all trials flagged by the detection pipeline,
#' suitable for downstream filtering via \code{dplyr::anti_join} or similar.
#'
#' @param x An \code{inough_detected} object.
#' @return Data frame with columns \code{id}, \code{trial_idx},
#'   \code{flag_type} (\code{"bailout"} or \code{"chunk"}), \code{chunk_id}
#'   (integer, \code{NA} for bailout), \code{p_adj} (numeric, \code{NA} for
#'   bailout).
#' @export
flags <- function(x) UseMethod("flags")

#' @export
flags.inough_detected <- function(x) {
  trial <- x$signals$trial
  rows  <- list()

  # Bail-out: all trials for bailed-out participants
  if (nrow(x$bailout) > 0L) {
    for (bid in x$bailout$id) {
      bt <- trial[trial$id == bid, ]
      if (nrow(bt) > 0L) {
        rows[[length(rows) + 1L]] <- data.frame(
          id        = bt$id,
          trial_idx = bt$trial_idx,
          flag_type = "bailout",
          chunk_id  = NA_integer_,
          p_adj     = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Significant chunks
  sig <- x$results[x$results$significant, , drop = FALSE]
  if (nrow(sig) > 0L) {
    for (j in seq_len(nrow(sig))) {
      pid <- sig$id[j]
      pt  <- trial[trial$id == pid, ]
      in_chunk <- pt$trial_idx >= sig$start[j] & pt$trial_idx <= sig$end[j]
      ct  <- pt[in_chunk, ]
      if (nrow(ct) > 0L) {
        rows[[length(rows) + 1L]] <- data.frame(
          id        = ct$id,
          trial_idx = ct$trial_idx,
          flag_type = "chunk",
          chunk_id  = j,
          p_adj     = sig$p_adj[j],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      id = character(0L), trial_idx = integer(0L),
      flag_type = character(0L), chunk_id = integer(0L),
      p_adj = numeric(0L), stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
