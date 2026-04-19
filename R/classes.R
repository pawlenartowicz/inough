# S3 constructors and methods for inough objects

# --- inough_signals ---

new_inough_signals <- function(trial, participant, model, call, formula) {
  stopifnot(is.data.frame(trial), is.data.frame(participant))
  structure(
    list(
      trial       = trial,
      participant = participant,
      model       = model,
      call        = call,
      formula     = formula
    ),
    class = "inough_signals"
  )
}

#' @export
print.inough_signals <- function(x, ...) {
  n_part  <- nrow(x$participant)
  n_trial <- nrow(x$trial)

  conv <- if (inherits(x$model, "glmerMod")) {
    if (is.null(x$model@optinfo$conv$lme4$messages)) "converged" else "converged (with warnings)"
  } else {
    "glm (no random effects)"
  }

  cat(sprintf("inough_signals: %d participants, %s trials\n",
              n_part, format(n_trial, big.mark = ",")))
  cat(sprintf("  Model: %s\n", conv))
  cat(sprintf("  LZ range: [%.3f, %.3f]\n",
              min(x$participant$lz), max(x$participant$lz)))
  cat(sprintf("  Accuracy range: [%.3f, %.3f]\n",
              min(x$participant$mean_accuracy), max(x$participant$mean_accuracy)))
  invisible(x)
}

#' @export
summary.inough_signals <- function(object, ...) {
  object$participant
}

# --- inough_detected ---

new_inough_detected <- function(signals, bailout, candidates, chunks,
                                results, control, heuristics, fdr_alpha,
                                call) {
  structure(
    list(
      signals    = signals,
      bailout    = bailout,
      candidates = candidates,
      chunks     = chunks,
      results    = results,
      control    = control,
      heuristics = heuristics,
      fdr_alpha  = fdr_alpha,
      call       = call
    ),
    class = "inough_detected"
  )
}

#' @export
print.inough_detected <- function(x, ...) {
  n_part  <- nrow(x$signals$participant)
  n_trial <- nrow(x$signals$trial)
  n_bail  <- nrow(x$bailout)

  sig <- x$results[x$results$significant, , drop = FALSE]
  n_sig_chunks    <- nrow(sig)
  n_flagged_parts <- length(unique(sig$id))

  fl <- flags(x)
  n_flagged_trials <- nrow(fl)

  # Bail-out reasons
  bail_str <- if (n_bail > 0L) {
    tbl <- table(x$bailout$reason)
    paste(vapply(names(tbl), function(r) sprintf("%d %s", tbl[[r]], r),
                 character(1L)), collapse = ", ")
  } else ""

  thresh_str <- if (is.character(x$control$screening_threshold)) {
    x$control$screening_threshold
  } else {
    format(x$control$screening_threshold)
  }

  cat(sprintf("inough detection: %d participants, %s trials\n",
              n_part, format(n_trial, big.mark = ",")))
  if (n_bail > 0L) {
    cat(sprintf("  Bailed out: %d (%.1f%%) \u2014 %s\n",
                n_bail, 100 * n_bail / n_part, bail_str))
  } else {
    cat("  Bailed out: 0\n")
  }
  cat(sprintf("  Flagged chunks: %d across %d participants (%.1f%%)\n",
              n_sig_chunks, n_flagged_parts, 100 * n_flagged_parts / n_part))
  cat(sprintf("  Flagged trials: %s (%.1f%%)\n",
              format(n_flagged_trials, big.mark = ","),
              100 * n_flagged_trials / n_trial))
  cat(sprintf("  Control: window=%d, screening=%s, fdr=%.2f\n",
              x$control$window_size, thresh_str, x$fdr_alpha))
  cat(sprintf("  Heuristics: boundary=%s, spurious=%s\n",
              x$heuristics$boundary_mode,
              if (x$heuristics$spurious) sprintf("yes (%d/%d)",
                x$heuristics$spurious_k, x$heuristics$spurious_n) else "no"))
  invisible(x)
}

#' @export
summary.inough_detected <- function(object, ...) {
  parts <- object$signals$participant
  fl <- flags(object)

  # Per-participant flagged trial count
  if (nrow(fl) > 0L) {
    flag_agg <- stats::aggregate(fl$trial_idx, by = list(id = fl$id), FUN = length)
    names(flag_agg)[2L] <- "n_flagged"
  } else {
    flag_agg <- data.frame(id = character(0L), n_flagged = integer(0L),
                           stringsAsFactors = FALSE)
  }

  # Significant chunk count per participant
  sig <- object$results[object$results$significant, , drop = FALSE]
  if (nrow(sig) > 0L) {
    chunk_agg <- stats::aggregate(sig$significant, by = list(id = sig$id), FUN = sum)
    names(chunk_agg)[2L] <- "n_chunks"
  } else {
    chunk_agg <- data.frame(id = character(0L), n_chunks = integer(0L),
                            stringsAsFactors = FALSE)
  }

  out <- parts[, c("id", "lz", "n_trials", "mean_accuracy")]
  out <- merge(out, flag_agg, by = "id", all.x = TRUE)
  out <- merge(out, chunk_agg, by = "id", all.x = TRUE)
  out$n_flagged[is.na(out$n_flagged)] <- 0L
  out$n_chunks[is.na(out$n_chunks)]   <- 0L
  out$pct_flagged <- round(100 * out$n_flagged / out$n_trials, 2)

  bail_ids <- object$bailout$id
  out$status <- ifelse(out$id %in% bail_ids, "bailout",
                ifelse(out$n_chunks > 0L, "flagged", "clean"))

  out <- out[order(-out$pct_flagged), ]
  out <- out[, c("id", "status", "n_chunks", "pct_flagged", "lz", "mean_accuracy")]
  rownames(out) <- NULL
  out
}
