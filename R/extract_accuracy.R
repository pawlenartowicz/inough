#' Extract accuracy residuals via probit GLMM
#'
#' Fits a probit GLMM (\code{lme4::glmer}) or probit GLM (\code{stats::glm})
#' on accuracy and returns Pearson residuals.
#'
#' @param df Data frame with \code{correct} (0/1), \code{id}, and predictor
#'   columns. Must already contain \code{n_trial} if the formula references it.
#' @param formula Full model formula (built by \code{inough_signals}).
#' @param has_random Logical; whether the formula includes random effects.
#' @return A list with \code{$residuals} (Pearson) and \code{$model}.
#' @keywords internal
extract_accuracy <- function(df, formula, has_random = TRUE) {
  if (!has_random) {
    model <- stats::glm(formula, data = df,
                        family = stats::binomial(link = "probit"))
    return(list(
      residuals = stats::residuals(model, type = "pearson"),
      model     = model
    ))
  }

  # Try default optimizer
  model <- suppressWarnings(
    lme4::glmer(formula, data = df,
                family = stats::binomial(link = "probit"))
  )

  msgs <- model@optinfo$conv$lme4$messages
  if (!is.null(msgs) && length(msgs) > 0L) {
    # Retry with bobyqa
    model <- suppressWarnings(
      lme4::glmer(formula, data = df,
                  family = stats::binomial(link = "probit"),
                  control = lme4::glmerControl(optimizer = "bobyqa"))
    )
    msgs2 <- model@optinfo$conv$lme4$messages
    if (!is.null(msgs2) && length(msgs2) > 0L) {
      warning("GLMM did not converge with default or bobyqa optimizer. ",
              "Residuals may be unreliable. Messages: ",
              paste(msgs2, collapse = "; "), call. = FALSE)
    }
  }

  list(
    residuals = stats::residuals(model, type = "pearson"),
    model     = model
  )
}
