#' Stage 1: Candidate Screening
#'
#' Flags regions where the absolute rolling mean of the lag-1 response
#' exceeds the per-participant threshold: \code{|roll_resp| > threshold}.
#' Catches both repetition (positive) and switching (negative) stereotypy
#' with a single track. Overlapping regions are merged; chunks shorter than
#' \code{min_chunk} are discarded. Boundary extension mode is controlled by
#' \code{heuristics$boundary_mode}.
#'
#' @param trial Trial data frame from \code{inough_signals}.
#' @param participant Participant data frame from \code{inough_signals}.
#' @param bail_ids Character vector of bailed-out participant IDs.
#' @param control \code{inough_control} object.
#' @param heuristics \code{inough_heuristics} object.
#' @return A list with \code{$candidates} (raw threshold crossings) and
#'   \code{$chunks} (after merge + min-length filter + boundary extension).
#' @keywords internal
screen <- function(trial, participant, bail_ids, control, heuristics) {
  w          <- control$window_size
  min_trials <- 3L * (2L * w + 1L)
  thresh_val <- control$screening_threshold
  min_chunk  <- control$min_chunk

  # Build weight vector once
  positions <- seq_len(2L * w + 1L) - (w + 1L)
  weights   <- switch(control$window_weight,
    uniform    = rep(1, 2L * w + 1L),
    triangular = (w + 1L) - abs(positions)
  )

  # Extension amount: triangular extends less (detected boundary closer to true onset)
  ext <- if (control$window_weight == "triangular") {
    as.integer(round(2 / 3 * w))
  } else {
    w
  }

  boundary_mode <- heuristics$boundary_mode

  ids <- setdiff(unique(as.character(trial$id)), bail_ids)

  all_candidates <- list()
  all_chunks     <- list()

  for (i in seq_along(ids)) {
    pid  <- ids[i]
    rows <- which(trial$id == pid)
    n    <- length(rows)
    if (n < min_trials) next

    resp_lag1 <- trial$resp_lag1[rows]
    tidx      <- trial$trial_idx[rows]

    roll_resp <- rolling_wmean(resp_lag1, weights)
    abs_resp  <- abs(roll_resp)

    cands <- rle_regions(abs_resp > thresh_val, tidx)
    if (nrow(cands) == 0L) next

    cands$id <- pid
    all_candidates[[length(all_candidates) + 1L]] <- cands

    # Merge overlapping + discard short
    merged <- merge_regions(cands[, c("start", "end")])
    merged <- merged[merged$end - merged$start + 1L >= min_chunk, , drop = FALSE]

    if (nrow(merged) > 0L) {
      # End boundary: always fixed extension
      merged$end <- pmin(tidx[length(tidx)], merged$end + ext)

      # Start boundary: depends on mode
      if (boundary_mode == "heuristic") {
        for (j in seq_len(nrow(merged))) {
          # Determine stereotype direction from rolling mean inside detected chunk
          chunk_mask <- tidx >= merged$start[j] & tidx <= merged$end[j]
          direction  <- sign(mean(roll_resp[chunk_mask]))

          if (direction == 0) {
            # Ambiguous — fall back to fixed
            merged$start[j] <- pmax(tidx[1L], merged$start[j] - ext)
          } else {
            # Walk backwards contiguously while resp_lag1 matches stereotype
            start_pos <- which(tidx == merged$start[j])
            candidate_start <- merged$start[j]
            for (step in seq_len(ext)) {
              check_pos <- start_pos - step
              if (check_pos < 1L) break
              if (resp_lag1[check_pos] == direction) {
                candidate_start <- tidx[check_pos]
              } else {
                break
              }
            }
            merged$start[j] <- candidate_start
          }
        }
      } else {
        # Fixed: symmetric extension
        merged$start <- pmax(tidx[1L], merged$start - ext)
      }

      # Clip to trial range and re-merge (expanded chunks may overlap)
      merged$start <- pmax(tidx[1L], merged$start)
      merged <- merge_regions(merged)
      merged$id <- pid
      all_chunks[[length(all_chunks) + 1L]] <- merged
    }
  }

  candidates <- if (length(all_candidates) > 0L) {
    do.call(rbind, all_candidates)
  } else {
    data.frame(id = character(0L), start = integer(0L), end = integer(0L),
               stringsAsFactors = FALSE)
  }

  chunks <- if (length(all_chunks) > 0L) {
    out <- do.call(rbind, all_chunks)
    out[, c("id", "start", "end")]
  } else {
    data.frame(id = character(0L), start = integer(0L), end = integer(0L),
               stringsAsFactors = FALSE)
  }

  list(candidates = candidates, chunks = chunks)
}


#' Apply spurious-start accuracy heuristic
#'
#' For each chunk, checks if the first \code{n} trials have suspiciously
#' high accuracy (\code{>= k} correct). If so, marks the chunk and shifts
#' the test start past those trials. Chunks with fewer than \code{min_left}
#' trials remaining are dropped.
#'
#' @param chunks Data frame with \code{id}, \code{start}, \code{end}.
#' @param trial Trial data frame (needs \code{id}, \code{trial_idx}, \code{correct}).
#' @param heuristics \code{inough_heuristics} object.
#' @return Updated chunks data frame with added columns: \code{test_start},
#'   \code{spurious_start}.
#' @keywords internal
apply_spurious_heuristic <- function(chunks, trial, heuristics) {
  if (nrow(chunks) == 0L) {
    chunks$test_start     <- integer(0L)
    chunks$spurious_start <- logical(0L)
    return(chunks)
  }

  n_check  <- heuristics$spurious_n
  k_thresh <- heuristics$spurious_k
  min_left <- heuristics$min_left

  chunks$test_start     <- chunks$start
  chunks$spurious_start <- FALSE

  for (i in seq_len(nrow(chunks))) {
    pid <- chunks$id[i]
    pt  <- trial[trial$id == pid, ]
    ch_trials <- pt[pt$trial_idx >= chunks$start[i] &
                    pt$trial_idx <= chunks$end[i], ]

    if (nrow(ch_trials) < n_check) next

    n_correct <- sum(ch_trials$correct[seq_len(n_check)])
    if (n_correct >= k_thresh) {
      chunks$spurious_start[i] <- TRUE
      chunks$test_start[i]     <- ch_trials$trial_idx[n_check + 1L]
    }
  }

  # Drop chunks where remaining length < min_left
  remaining <- chunks$end - chunks$test_start + 1L
  chunks <- chunks[remaining >= min_left, , drop = FALSE]
  rownames(chunks) <- NULL

  chunks
}


# --- Helpers (package-internal) ---

#' Centered weighted rolling mean
#'
#' \code{weights} is the full weight vector of length \code{2*k+1} (centre at
#' position \code{k+1}). At the edges the weight vector is cropped and
#' renormalized so the SD formula stays valid.
#' @keywords internal
rolling_wmean <- function(x, weights) {
  n  <- length(x)
  k  <- (length(weights) - 1L) %/% 2L
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo  <- max(1L, i - k)
    hi  <- min(n, i + k)
    wi  <- weights[(lo - i + k + 1L):(hi - i + k + 1L)]
    out[i] <- sum(wi * x[lo:hi]) / sum(wi)
  }
  out
}

#' Robust z-score (median / MAD); returns NULL if MAD = 0
#' @keywords internal
robust_zscore <- function(x) {
  med <- stats::median(x)
  mad_val <- stats::mad(x)
  if (mad_val == 0) return(NULL)
  (x - med) / mad_val
}

#' Contiguous TRUE regions to start/end data.frame
#' @keywords internal
rle_regions <- function(above, trial_idx) {
  r <- rle(above)
  if (!any(r$values)) {
    return(data.frame(start = integer(0L), end = integer(0L),
                      stringsAsFactors = FALSE))
  }
  ends   <- cumsum(r$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  ok     <- r$values
  data.frame(
    start = trial_idx[starts[ok]],
    end   = trial_idx[ends[ok]],
    stringsAsFactors = FALSE
  )
}

#' Sort and merge overlapping regions
#' @keywords internal
merge_regions <- function(regions) {
  if (nrow(regions) == 0L) return(regions)
  regions <- regions[order(regions$start), , drop = FALSE]

  merged_start <- regions$start[1L]
  merged_end   <- regions$end[1L]
  out_start <- integer(0L)
  out_end   <- integer(0L)

  for (i in seq_len(nrow(regions))[-1L]) {
    if (regions$start[i] <= merged_end + 1L) {
      # Overlapping or adjacent — extend
      merged_end <- max(merged_end, regions$end[i])
    } else {
      out_start <- c(out_start, merged_start)
      out_end   <- c(out_end,   merged_end)
      merged_start <- regions$start[i]
      merged_end   <- regions$end[i]
    }
  }
  out_start <- c(out_start, merged_start)
  out_end   <- c(out_end,   merged_end)

  data.frame(start = out_start, end = out_end, stringsAsFactors = FALSE)
}
