#' Detect inattention episodes
#'
#' Two-stage pipeline: (1) dual-track screening with chunk filtering,
#' (2) formal t-test with FDR correction. Participants with extremely
#' stereotyped responses or chance-level accuracy are bailed out first.
#'
#' @param signals An \code{inough_signals} object from
#'   \code{\link{inough_signals}}.
#' @param fdr_alpha FDR significance level for BH correction (default 0.05).
#' @param control An \code{inough_control} object (see
#'   \code{\link{inough_control}}).
#' @param heuristics An \code{inough_heuristics} object (see
#'   \code{\link{inough_heuristics}}). Controls boundary extension mode
#'   and spurious-accuracy trimming.
#' @return An \code{inough_detected} object.
#' @export
inough_detect <- function(signals, fdr_alpha = 0.2,
                          control    = inough_control(),
                          heuristics = inough_heuristics()) {
  stopifnot(inherits(signals, "inough_signals"))
  stopifnot(is.numeric(fdr_alpha), length(fdr_alpha) == 1L,
            fdr_alpha > 0, fdr_alpha < 1)
  stopifnot(inherits(control, "inough_control"))
  stopifnot(inherits(heuristics, "inough_heuristics"))

  # Pre-pipeline bail-out
  bail     <- bailout(signals$participant, control$lz_threshold)
  bail_ids <- bail$id

  # Stage 1: Screening + merge/filter + boundary extension
  screened   <- screen(signals$trial, signals$participant, bail_ids,
                       control, heuristics)
  candidates <- screened$candidates
  chunks     <- screened$chunks

  # Spurious-accuracy heuristic (between screening and testing)
  if (heuristics$spurious) {
    chunks <- apply_spurious_heuristic(chunks, signals$trial, heuristics)
  }

  # Stage 2: Formal test (Welch t + local FDR)
  results <- test_chunks(chunks, signals$trial, signals$participant,
                         fdr_alpha, control$comparison)

  new_inough_detected(
    signals    = signals,
    bailout    = bail,
    candidates = candidates,
    chunks     = chunks,
    results    = results,
    control    = control,
    heuristics = heuristics,
    fdr_alpha  = fdr_alpha,
    call       = match.call()
  )
}
