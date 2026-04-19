#' Post-detection heuristics
#'
#' Configures optional refinements applied after chunk screening:
#' boundary extension mode and spurious-accuracy trimming.
#'
#' @param boundary_mode How to extend chunk boundaries after detection:
#'   \code{"heuristic"} (default) walks backwards from chunk start to find
#'   where stereotyped responding actually began. \code{"fixed"} extends symmetrically (amount depends on
#'   \code{window_weight} in \code{inough_control}), \code{"heuristic"}
#'   walks backwards from chunk start to find where stereotyped responding
#'   actually began. End boundary always uses fixed extension.
#' @param spurious Logical; enable spurious-start accuracy trimming
#'   (default \code{TRUE}). When \code{TRUE}, the first \code{spurious_n}
#'   trials of each chunk are checked: if accuracy is suspiciously high
#'   (\code{>= k} correct), those trials are excluded from the t-test.
#' @param spurious_n Number of trials at chunk start to inspect (default 6).
#' @param spurious_k Explicit threshold: flag if \code{>= k} correct out of
#'   \code{spurious_n}. When \code{NULL} (default), computed as
#'   \code{ceiling((0.5 + 1.5 * sqrt(0.25 / n)) * n)}.
#' @param min_left Minimum trials remaining after spurious trimming
#'   (default 7). Chunks shorter than this after trimming are dropped.
#' @return An \code{inough_heuristics} object.
#' @export
inough_heuristics <- function(boundary_mode = "heuristic",
                              spurious      = TRUE,
                              spurious_n    = 6L,
                              spurious_k    = NULL,
                              min_left      = 7L) {
  boundary_mode <- match.arg(boundary_mode, c("fixed", "heuristic"))
  stopifnot(is.logical(spurious), length(spurious) == 1L)

  n <- as.integer(spurious_n)
  stopifnot(n >= 2L)
  stopifnot(is.numeric(min_left), length(min_left) == 1L, min_left >= 1L)

  if (spurious) {
    if (!is.null(spurious_k)) {
      k <- as.integer(spurious_k)
      stopifnot(k >= 1L, k <= n)
    } else {
      se <- sqrt(0.25 / n)
      k  <- as.integer(ceiling((0.5 + 1.5 * se) * n))
      k  <- min(k, n)
    }
  } else {
    k <- NA_integer_
  }

  structure(
    list(
      boundary_mode = boundary_mode,
      spurious      = spurious,
      spurious_n    = n,
      spurious_k    = k,
      min_left      = as.integer(min_left)
    ),
    class = "inough_heuristics"
  )
}

#' @export
print.inough_heuristics <- function(x, ...) {
  cat("inough_heuristics:\n")
  cat("  boundary_mode =", x$boundary_mode, "\n")
  if (x$spurious) {
    cat("  spurious      = TRUE\n")
    cat("  spurious_n    =", x$spurious_n, "\n")
    cat("  spurious_k    =", x$spurious_k,
        sprintf("(flag if >= %d/%d correct)\n", x$spurious_k, x$spurious_n))
    cat("  min_left      =", x$min_left, "\n")
  } else {
    cat("  spurious      = FALSE\n")
  }
  invisible(x)
}
