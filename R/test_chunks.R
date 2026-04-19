#' Stage 3: Formal Test
#'
#' Welch t-test on accuracy residuals (chunk vs comparison set) with
#' LZ-informed local FDR. When \code{test_start} column is present (from
#' spurious-accuracy trimming), uses it for inside trials while still
#' excluding the full chunk from the comparison set.
#'
#' Local FDR is computed per chunk:
#' \deqn{lfdr = \pi_0 f_0(t) / [\pi_0 f_0(t) + \pi_1 f_1(t)]}
#' where \eqn{\pi_0 = \min(1, LZ)} is the LZ-informed null prior,
#' \eqn{f_0} is the central t density (null), and \eqn{f_1} is
#' a non-central t density with noncentrality parameter reflecting the
#' expected accuracy drop to chance.
#'
#' @param chunks Data frame from Stage 2 (\code{id}, \code{start}, \code{end},
#'   and optionally \code{test_start}, \code{spurious_start}).
#' @param trial Trial data frame from \code{inough_signals}.
#' @param participant Participant data frame from \code{inough_signals}.
#' @param fdr_alpha Significance threshold for local FDR (default 0.05).
#'   A chunk is flagged when \code{lfdr < fdr_alpha}.
#' @param comparison \code{"clean"} or \code{"rest"}.
#' @return Data frame with \code{id}, \code{start}, \code{end}, \code{t_stat},
#'   \code{df}, \code{p_raw}, \code{p_adj}, \code{lfdr}, \code{effect_size},
#'   \code{significant}.
#' @keywords internal
test_chunks <- function(chunks, trial, participant, fdr_alpha, comparison) {
  if (nrow(chunks) == 0L) {
    return(data.frame(
      id = character(0L), start = integer(0L), end = integer(0L),
      t_stat = numeric(0L), df = numeric(0L), p_raw = numeric(0L),
      p_adj = numeric(0L), lfdr = numeric(0L), effect_size = numeric(0L),
      mean_acc_inside = numeric(0L), mean_acc_outside = numeric(0L),
      significant = logical(0L), stringsAsFactors = FALSE
    ))
  }

  has_spurious <- "test_start" %in% names(chunks)

  part_ids <- unique(chunks$id)
  results  <- vector("list", length(part_ids))

  for (pi in seq_along(part_ids)) {
    pid <- part_ids[pi]

    part_rows    <- which(trial$id == pid)
    part_idx     <- trial$trial_idx[part_rows]
    part_resid   <- trial$acc_resid[part_rows]
    part_correct <- trial$correct[part_rows]

    pid_chunks <- chunks[chunks$id == pid, , drop = FALSE]

    # LZ-informed prior for local FDR
    lz_val <- participant$lz[participant$id == pid]
    pi0 <- min(1, max(0, lz_val))
    pi1 <- 1 - pi0

    # For "clean" comparison: exclude trials in ANY chunk for this participant
    if (comparison == "clean") {
      in_any <- rep(FALSE, length(part_idx))
      for (j in seq_len(nrow(pid_chunks))) {
        in_any <- in_any | (part_idx >= pid_chunks$start[j] &
                            part_idx <= pid_chunks$end[j])
      }
      acc_clean <- mean(part_correct[!in_any])
    } else {
      acc_clean <- participant$mean_accuracy[participant$id == pid]
    }

    # Outside accuracy for display: non-chunked trials for this participant.
    # When spurious is active, spurious regions (start:test_start) stay in "outside".
    in_any_display <- rep(FALSE, length(part_idx))
    for (j in seq_len(nrow(pid_chunks))) {
      ch_s <- if (has_spurious) pid_chunks$test_start[j] else pid_chunks$start[j]
      in_any_display <- in_any_display | (part_idx >= ch_s &
                                          part_idx <= pid_chunks$end[j])
    }
    mean_acc_outside <- mean(part_correct[!in_any_display])

    chunk_res <- vector("list", nrow(pid_chunks))
    for (j in seq_len(nrow(pid_chunks))) {
      ch_start <- pid_chunks$start[j]
      ch_end   <- pid_chunks$end[j]

      # Test region: use test_start if available (spurious trim)
      test_start <- if (has_spurious) pid_chunks$test_start[j] else ch_start
      in_chunk <- part_idx >= test_start & part_idx <= ch_end

      inside  <- part_resid[in_chunk]
      outside <- if (comparison == "clean") {
        part_resid[!in_any]
      } else {
        part_resid[!in_chunk]
      }

      n_in  <- length(inside)
      n_out <- length(outside)

      if (n_in < 2L || n_out < 2L) {
        res_row <- data.frame(
          id = pid, start = ch_start, end = ch_end,
          t_stat = NA_real_, df = NA_real_, p_raw = NA_real_,
          effect_size = NA_real_, lfdr = NA_real_,
          mean_acc_inside = mean(part_correct[in_chunk]),
          mean_acc_outside = mean_acc_outside,
          stringsAsFactors = FALSE
        )
        if (has_spurious) {
          res_row$test_start     <- test_start
          res_row$spurious_start <- pid_chunks$spurious_start[j]
        }
        chunk_res[[j]] <- res_row
        next
      }

      tt <- stats::t.test(inside, outside, var.equal = FALSE,
                          alternative = "less")

      t_val  <- tt$statistic[[1L]]
      df_val <- tt$parameter[[1L]]

      # Cohen's d (pooled SD)
      pooled_sd <- sqrt((stats::var(inside) + stats::var(outside)) / 2)
      d <- if (pooled_sd > 0) (mean(inside) - mean(outside)) / pooled_sd else 0

      # Local FDR: LZ-informed Bayesian
      # f0: null density — central t at observed df
      f0 <- stats::dt(t_val, df = df_val)
      # f1: alternative density — non-central t
      #   ncp = expected accuracy drop / theoretical SE
      #   Under H1 chunk accuracy drops to chance (0.5)
      se_theory <- sqrt(0.25 / n_in +
                        acc_clean * (1 - acc_clean) / n_out)
      ncp_val   <- -(acc_clean - 0.5) / se_theory
      ncp_val   <- max(min(ncp_val, 37.62), -37.62)
      f1        <- suppressWarnings(
        stats::dt(t_val, df = df_val, ncp = ncp_val)
      )

      denom    <- pi0 * f0 + pi1 * f1
      lfdr_val <- if (denom > 0) pi0 * f0 / denom else 1

      res_row <- data.frame(
        id = pid, start = ch_start, end = ch_end,
        t_stat = t_val, df = df_val,
        p_raw = tt$p.value, effect_size = d, lfdr = lfdr_val,
        mean_acc_inside = mean(part_correct[in_chunk]),
        mean_acc_outside = mean_acc_outside,
        stringsAsFactors = FALSE
      )
      if (has_spurious) {
        res_row$test_start     <- test_start
        res_row$spurious_start <- pid_chunks$spurious_start[j]
      }
      chunk_res[[j]] <- res_row
    }

    pid_df <- do.call(rbind, chunk_res)
    pid_df$p_adj <- stats::p.adjust(pid_df$p_raw, method = "BH")
    pid_df$significant <- !is.na(pid_df$lfdr) & pid_df$lfdr < fdr_alpha
    results[[pi]] <- pid_df
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
