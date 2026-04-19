#' Pipeline tuning parameters
#'
#' @param lz_threshold LZ complexity below this triggers bail-out (default 0.2).
#' @param window_size Half-width of rolling window; total window = 2*w+1
#'   (default 3, giving 7-trial windows).
#' @param sd_threshold Number of SDs above chance that \code{|roll_resp|} must
#'   exceed to flag a candidate region. Under H0 of balanced random responding,
#'   \code{SD = sqrt(sum(w^2)) / sum(w)} where \code{w} are the window weights.
#'   For uniform weights this simplifies to \code{1/sqrt(2*window_size+1)}.
#'   Default 2 (about 5\% false-positive rate per window under H0).
#' @param window_weight Weighting scheme for the rolling window:
#'   \code{"uniform"} (default) or \code{"triangular"} (centre trials weighted
#'   more; reduces edge sensitivity). The SD is computed analytically for both.
#' @param min_chunk Minimum chunk length in trials to retain after merging
#'   (default 6).
#' @param comparison t-test comparison set: \code{"clean"} (chunk vs all
#'   non-flagged trials) or \code{"rest"} (chunk vs everything except that
#'   chunk). Default \code{"clean"}.
#' @return An \code{inough_control} object with a pre-computed
#'   \code{screening_threshold} based on \code{sd_threshold} and the window.
#' @export
inough_control <- function(lz_threshold   = 0.2,
                           window_size    = 3,
                           sd_threshold   = 2,
                           window_weight  = "uniform",
                           min_chunk      = 6,
                           comparison     = "clean") {
  stopifnot(is.numeric(lz_threshold), length(lz_threshold) == 1L,
            lz_threshold > 0, lz_threshold < 1)
  stopifnot(is.numeric(window_size), length(window_size) == 1L,
            window_size >= 1L)
  stopifnot(is.numeric(sd_threshold), length(sd_threshold) == 1L,
            sd_threshold > 0)
  window_weight <- match.arg(window_weight, c("uniform", "triangular"))
  stopifnot(is.numeric(min_chunk), length(min_chunk) == 1L, min_chunk >= 1L)
  stopifnot(comparison %in% c("clean", "rest"))

  w <- as.integer(window_size)
  positions <- seq_len(2L * w + 1L) - (w + 1L)  # -w ... 0 ... w

  weights <- switch(window_weight,
    uniform    = rep(1, 2L * w + 1L),
    triangular = (w + 1L) - abs(positions)
  )

  # Analytical SD of weighted mean of ±1 values under H0 (p = 0.5)
  window_sd          <- sqrt(sum(weights^2)) / sum(weights)
  screening_threshold <- sd_threshold * window_sd

  structure(
    list(
      lz_threshold        = lz_threshold,
      window_size         = w,
      sd_threshold        = sd_threshold,
      window_weight       = window_weight,
      window_sd           = window_sd,
      screening_threshold = screening_threshold,
      min_chunk           = as.integer(min_chunk),
      comparison          = comparison
    ),
    class = "inough_control"
  )
}

#' @export
print.inough_control <- function(x, ...) {
  cat("inough_control:\n")
  for (nm in names(x)) {
    cat("  ", nm, "=", format(x[[nm]]), "\n")
  }
  invisible(x)
}
