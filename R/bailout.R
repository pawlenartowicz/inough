#' Participant-level bail-out
#'
#' Flags entire sessions where trial-level analysis is unreliable: response
#' sequence too stereotyped (low LZ) or accuracy indistinguishable from chance.
#'
#' @param participant Data frame with \code{id}, \code{lz}, \code{n_trials},
#'   \code{mean_accuracy}.
#' @param lz_threshold LZ threshold (from \code{inough_control}).
#' @return Data frame with \code{id}, \code{reason}, \code{lz_value},
#'   \code{accuracy} for bailed-out participants (zero rows if none).
#' @keywords internal
bailout <- function(participant, lz_threshold) {
  low_lz <- participant$lz < lz_threshold

  # 2AFC: chance = 0.5, SE = sqrt(p*(1-p)/n)
  se <- sqrt(0.25 / participant$n_trials)
  at_chance <- participant$mean_accuracy <= (0.5 + se)

  flagged <- low_lz | at_chance

  if (!any(flagged)) {
    return(data.frame(id = character(0L), reason = character(0L),
                      lz_value = numeric(0L), accuracy = numeric(0L),
                      stringsAsFactors = FALSE))
  }

  reason <- ifelse(low_lz & at_chance, "both",
             ifelse(low_lz, "lz", "chance"))

  data.frame(
    id       = participant$id[flagged],
    reason   = reason[flagged],
    lz_value = participant$lz[flagged],
    accuracy = participant$mean_accuracy[flagged],
    stringsAsFactors = FALSE
  )
}
