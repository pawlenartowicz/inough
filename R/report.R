#' Generate interactive HTML report
#'
#' Creates a self-contained HTML file with a participant browser, diagnostic
#' plots, chunk details, and summary statistics.
#'
#' @param x An \code{inough_detected} object.
#' @param ... Arguments passed to methods.
#' @return Invisibly returns the file path.
#' @export
report <- function(x, ...) UseMethod("report")

#' @rdname report
#' @param file Output file path. If \code{NULL} (default), uses a tempfile and
#'   opens in the browser.
#' @param custom_plot Optional per-trial variable to show as an extra panel in
#'   the participant view. A list with three fields:
#'   \itemize{
#'     \item \code{data}: a data.frame with columns \code{id}, \code{trial_idx},
#'       \code{value}. \code{trial_idx} is 1-based within participant, matching
#'       the trial order fed to \code{inough_signals}.
#'     \item \code{title}: string shown as the panel label.
#'     \item \code{description}: string shown as the panel blurb.
#'   }
#' @export
report.inough_detected <- function(x, file = NULL, custom_plot = NULL, ...) {
  open_browser <- is.null(file)
  if (is.null(file)) file <- tempfile(fileext = ".html")

  custom_plot <- validate_custom_plot(custom_plot)

  # --- SUMMARY JSON ---
  participants_list <- build_participant_summary(x)
  summary_data <- list(
    meta = list(
      n_participants  = nrow(x$signals$participant),
      n_trials        = nrow(x$signals$trial),
      package_version = as.character(utils::packageVersion("inough")),
      generated       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      control         = c(as.list(x$control), list(fdr_alpha = x$fdr_alpha)),
      heuristics      = unclass(x$heuristics),
      aggregate       = build_aggregate_summary(participants_list),
      custom_plot     = if (is.null(custom_plot)) NULL else
        list(title = custom_plot$title, description = custom_plot$description)
    ),
    participants = participants_list,
    chunks       = build_chunk_summary(x)
  )
  summary_json <- jsonlite::toJSON(summary_data, auto_unbox = TRUE,
                                   digits = 4, pretty = FALSE, null = "null")

  # --- TRIAL DATA JSON (column-array per participant) ---
  trial_data <- build_trial_data(x, custom_plot = custom_plot)
  trial_json <- jsonlite::toJSON(trial_data, auto_unbox = FALSE,
                                  digits = 4, pretty = FALSE, na = "null")

  # --- Build HTML ---
  template_path <- system.file("report_template.html", package = "inough")
  if (template_path == "") {
    stop("Report template not found. Is the package installed correctly?")
  }
  template <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

  html <- gsub("/*{{SUMMARY_JSON}}*/{}",  summary_json, template, fixed = TRUE)
  html <- gsub("/*{{TRIAL_DATA_JSON}}*/{}", trial_json,  html, fixed = TRUE)

  writeLines(html, file)
  if (open_browser) utils::browseURL(file)

  invisible(file)
}


# --- internal builders ---

build_participant_summary <- function(x) {
  fl    <- flags(x)
  trial <- x$signals$trial
  chunk_fl <- fl[fl$flag_type == "chunk", , drop = FALSE]

  lapply(seq_len(nrow(x$signals$participant)), function(i) {
    p   <- x$signals$participant[i, ]
    pt  <- trial[trial$id == p$id, ]
    n_flagged <- sum(fl$id == p$id)
    sig <- x$results[x$results$id == p$id & x$results$significant, ]
    is_bail     <- p$id %in% x$bailout$id
    bail_reason <- if (is_bail) x$bailout$reason[x$bailout$id == p$id][1L] else NA_character_
    status      <- if (is_bail) "bailout" else if (nrow(sig) > 0L) "flagged" else "clean"

    # Trials removed due to *chunks only* (bailouts handled separately at aggregate level)
    chunk_idx    <- chunk_fl$trial_idx[chunk_fl$id == p$id]
    keep_mask    <- !(pt$trial_idx %in% chunk_idx)
    n_kept       <- sum(keep_mask)
    mean_acc_cor <- if (n_kept > 0L) mean(pt$correct[keep_mask]) else NA_real_

    # Raw (un-normalized) LZ pre/post for diagnostic comparison.
    # We avoid the permutation-normalized value here because the null depends
    # on sequence length, which differs between pre and post.
    lz_raw_pre  <- safe_lz(pt$response)
    lz_raw_post <- if (n_kept >= 2L) safe_lz(pt$response[keep_mask]) else NA_real_

    list(
      id                        = p$id,
      status                    = status,
      n_trials                  = p$n_trials,
      n_trials_corrected        = n_kept,
      mean_accuracy             = round(p$mean_accuracy, 4),
      mean_accuracy_corrected   = round_or_na(mean_acc_cor, 4),
      lz                        = round(p$lz, 4),
      lz_raw                    = round_or_na(lz_raw_pre, 4),
      lz_raw_corrected          = round_or_na(lz_raw_post, 4),
      n_chunks                  = nrow(sig),
      n_flagged_chunk           = length(chunk_idx),
      pct_flagged               = round(100 * n_flagged / p$n_trials, 2),
      bailout_reason            = bail_reason
    )
  })
}

build_aggregate_summary <- function(participants_list) {
  if (!length(participants_list)) {
    return(list(mean_acc_raw = NA_real_,
                mean_acc_no_chunks = NA_real_,
                mean_acc_no_chunks_no_lz_bail = NA_real_,
                n_bail_lz = 0L, n_bail_chance = 0L))
  }
  accs_raw  <- vapply(participants_list, function(q) as.numeric(q$mean_accuracy),            numeric(1L))
  accs_corr <- vapply(participants_list, function(q) as.numeric(q$mean_accuracy_corrected),  numeric(1L))
  reasons   <- vapply(participants_list, function(q) {
    r <- q$bailout_reason; if (is.null(r)) NA_character_ else as.character(r)
  }, character(1L))
  is_lz_bail     <- !is.na(reasons) & reasons %in% c("lz", "both")
  is_chance_bail <- !is.na(reasons) & reasons %in% c("chance", "both")

  list(
    mean_acc_raw                  = round_or_na(mean(accs_raw, na.rm = TRUE), 3),
    mean_acc_no_chunks            = round_or_na(mean(accs_corr, na.rm = TRUE), 3),
    mean_acc_no_chunks_no_lz_bail = round_or_na(mean(accs_corr[!is_lz_bail], na.rm = TRUE), 3),
    n_bail_lz                     = sum(is_lz_bail),
    n_bail_chance                 = sum(is_chance_bail & !is_lz_bail)
  )
}

round_or_na <- function(x, digits) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !is.finite(x)) NA_real_
  else round(x, digits)
}

safe_lz <- function(resp) {
  tryCatch(lz_complexity(as.integer(resp)), error = function(e) NA_real_)
}

build_chunk_summary <- function(x) {
  res <- x$results
  if (nrow(res) == 0L) return(list())

  has_spurious <- "spurious_start" %in% names(res)
  has_lfdr     <- "lfdr" %in% names(res)

  lapply(seq_len(nrow(res)), function(j) {
    r  <- res[j, ]
    out <- list(
      id               = r$id,
      chunk_id         = j,
      start            = r$start,
      end              = r$end,
      significant      = r$significant,
      t_stat           = round(r$t_stat, 3),
      p_raw            = round(r$p_raw, 4),
      p_adj            = round(r$p_adj, 4),
      effect_size      = round(r$effect_size, 3),
      mean_acc_inside  = round(r$mean_acc_inside, 3),
      mean_acc_outside = round(r$mean_acc_outside, 3)
    )
    if (has_lfdr)     out$lfdr           <- round(r$lfdr, 4)
    if (has_spurious) {
      out$spurious_start <- r$spurious_start
      out$test_start     <- r$test_start
    }
    out
  })
}

build_trial_data <- function(x, custom_plot = NULL) {
  w     <- x$control$window_size
  trial <- x$signals$trial
  ids   <- x$signals$participant$id

  out <- stats::setNames(vector("list", length(ids)), ids)
  for (pid in ids) {
    pt <- trial[trial$id == pid, ]

    positions <- seq_len(2L * w + 1L) - (w + 1L)
    weights   <- switch(x$control$window_weight,
      uniform    = rep(1, 2L * w + 1L),
      triangular = (w + 1L) - abs(positions)
    )
    roll_resp <- rolling_wmean(pt$resp_lag1, weights)

    entry <- list(
      trial_idx  = pt$trial_idx,
      correct    = pt$correct,
      acc_resid  = round(pt$acc_resid, 4),
      resp_lag1  = pt$resp_lag1,
      roll_resp  = round(roll_resp, 4),
      abs_resp   = round(abs(roll_resp), 4)
    )

    if (!is.null(custom_plot)) {
      cp_sub <- custom_plot$data[custom_plot$data$id == pid, , drop = FALSE]
      aligned <- rep(NA_real_, nrow(pt))
      if (nrow(cp_sub) > 0L) {
        m <- match(pt$trial_idx, cp_sub$trial_idx)
        aligned <- as.numeric(cp_sub$value)[m]
      }
      entry$custom <- round(aligned, 4)
    }

    out[[pid]] <- entry
  }
  out
}

validate_custom_plot <- function(cp) {
  if (is.null(cp)) return(NULL)
  if (!is.list(cp)) stop("'custom_plot' must be a list with fields data, title, description")
  required <- c("data", "title", "description")
  missing_f <- setdiff(required, names(cp))
  if (length(missing_f) > 0L) {
    stop("'custom_plot' is missing fields: ", paste(missing_f, collapse = ", "))
  }
  if (!is.data.frame(cp$data)) stop("'custom_plot$data' must be a data.frame")
  req_cols <- c("id", "trial_idx", "value")
  missing_c <- setdiff(req_cols, names(cp$data))
  if (length(missing_c) > 0L) {
    stop("'custom_plot$data' is missing columns: ", paste(missing_c, collapse = ", "))
  }
  if (!is.character(cp$title) || length(cp$title) != 1L) {
    stop("'custom_plot$title' must be a single string")
  }
  if (!is.character(cp$description) || length(cp$description) != 1L) {
    stop("'custom_plot$description' must be a single string")
  }
  cp$data <- data.frame(
    id        = as.character(cp$data$id),
    trial_idx = as.integer(cp$data$trial_idx),
    value     = as.numeric(cp$data$value),
    stringsAsFactors = FALSE
  )
  cp
}
