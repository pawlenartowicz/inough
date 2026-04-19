#' Extract inattention signals
#'
#' Fits a probit GLMM on accuracy and computes response bias indicators.
#' This is the first step in the inough pipeline: call \code{inough_signals},
#' then pass the result to \code{\link{inough_detect}}.
#'
#' @param df A data frame with rows ordered by trial within each participant.
#' @param correct Formula. LHS names the accuracy column (0/1 integer). RHS
#'   names design predictors that explain correctness (e.g.,
#'   \code{correct ~ Stim + Weight + Orient + Block}).
#' @param response Formula identifying the response column. Use a two-sided
#'   formula where the RHS names the column (e.g., \code{response ~ answ}).
#'   Must have exactly 2 distinct values; auto-encoded to 0/1 by sorted order.
#' @param id String naming the participant identifier column (default
#'   \code{"ID"}).
#' @param learning_effect Logical. If \code{TRUE} (default), adds \code{n_trial}
#'   as fixed effect and random slope per participant.
#' @param participant_effect Logical. If \code{TRUE} (default), adds random
#'   intercept per participant.
#' @param trial_transform Transformation applied to trial index before rescaling
#'   to [-1, 1]. One of \code{"sqrt"} (default), \code{"log"}, \code{"linear"},
#'   or a function.
#' @return An \code{inough_signals} object.
#' @export
inough_signals <- function(df, correct, response, id = "ID",
                           learning_effect = TRUE,
                           participant_effect = TRUE,
                           trial_transform = "sqrt") {
  stopifnot(is.data.frame(df))

  # --- Parse correct formula ---
  correct_fml <- stats::as.formula(correct)
  if (length(correct_fml) != 3L) {
    stop("'correct' must be a two-sided formula (e.g., correct ~ Stim + Block)")
  }
  correct_col   <- as.character(correct_fml[[2L]])
  predictor_vars <- all.vars(correct_fml[[3L]])
  rhs_str       <- deparse(correct_fml[[3L]], width.cutoff = 500L)


  # --- Parse response formula ---
  response_fml <- stats::as.formula(response)
  if (length(response_fml) == 3L) {
    response_col <- all.vars(response_fml[[3L]])
  } else {
    response_col <- all.vars(response_fml[[2L]])
  }
  if (length(response_col) != 1L) {
    stop("response formula must reference exactly one column, got: ",
         paste(response_col, collapse = ", "))
  }

  # --- Validate columns exist ---
  needed <- unique(c(correct_col, response_col, id, predictor_vars))
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in df: ", paste(missing_cols, collapse = ", "))
  }

  # --- Validate correct column ---
  if (!all(df[[correct_col]] %in% c(0L, 1L))) {
    stop("Column '", correct_col, "' must contain only 0 and 1")
  }

  # --- Validate response column ---
  resp_vals <- sort(unique(df[[response_col]]))
  if (length(resp_vals) != 2L) {
    stop("Column '", response_col, "' must have exactly 2 unique values, got ",
         length(resp_vals))
  }

  # --- Validate id column ---
  if (anyNA(df[[id]]))       stop("NAs found in '", id, "'")
  if (anyNA(df[[correct_col]])) stop("NAs found in '", correct_col, "'")
  if (anyNA(df[[response_col]])) stop("NAs found in '", response_col, "'")

  # --- NAs in predictors: warn and drop ---
  na_pred <- vapply(predictor_vars,
                    function(v) sum(is.na(df[[v]])), integer(1L))
  if (any(na_pred > 0L)) {
    bad <- na_pred[na_pred > 0L]
    warning("Dropping ", sum(bad), " rows with NAs in predictors: ",
            paste(names(bad), "=", bad, collapse = ", "), call. = FALSE)
    keep <- stats::complete.cases(df[, predictor_vars, drop = FALSE])
    df <- df[keep, ]
  }

  # --- Build internal data ---
  id_factor <- factor(df[[id]])
  ids <- levels(id_factor)

  # Auto-encode response to 0/1
  response_encoded <- ifelse(df[[response_col]] == resp_vals[1L], 0L, 1L)

  # Trial transform function
  tfn <- if (is.function(trial_transform)) {
    trial_transform
  } else {
    switch(trial_transform,
           sqrt = sqrt,
           log  = log,
           linear = identity,
           stop("trial_transform must be 'sqrt', 'log', 'linear', or a function"))
  }

  # Generate n_trial: sequential index per participant, transformed, rescaled to [-1,1]
  n_trial <- numeric(nrow(df))
  for (id_val in ids) {
    rows <- which(id_factor == id_val)
    idx  <- seq_along(rows)
    tx   <- tfn(idx)
    rng  <- range(tx)
    n_trial[rows] <- if (rng[1L] == rng[2L]) 0 else 2 * (tx - rng[1L]) / (rng[2L] - rng[1L]) - 1
  }

  # Assemble modelling data frame
  model_df <- df[, predictor_vars, drop = FALSE]
  model_df$correct <- df[[correct_col]]
  model_df$id      <- id_factor
  model_df$n_trial <- n_trial
  model_df$response <- response_encoded

  # --- Build GLMM formula ---
  has_random <- participant_effect || learning_effect

  fml_str <- paste("correct ~", rhs_str)
  if (learning_effect)  fml_str <- paste(fml_str, "+ n_trial")
  if (participant_effect && learning_effect)   fml_str <- paste(fml_str, "+ (1 + n_trial | id)")
  else if (participant_effect)                 fml_str <- paste(fml_str, "+ (1 | id)")
  else if (learning_effect)                    fml_str <- paste(fml_str, "+ (0 + n_trial | id)")
  model_formula <- stats::as.formula(fml_str)

  # --- Fit model ---
  acc  <- extract_accuracy(model_df, model_formula, has_random = has_random)
  bias <- extract_bias(data.frame(id = model_df$id, response = model_df$response))

  # --- Assemble trial data ---
  trial_idx <- integer(nrow(model_df))
  for (id_val in ids) {
    rows <- which(model_df$id == id_val)
    trial_idx[rows] <- seq_along(rows)
  }

  trial <- data.frame(
    id        = model_df$id,
    trial_idx = trial_idx,
    correct   = model_df$correct,
    response  = model_df$response,
    acc_resid = acc$residuals,
    resp_lag1 = bias$trial$resp_lag1,
    stringsAsFactors = FALSE
  )

  # --- Assemble participant data ---
  participant <- data.frame(
    id            = ids,
    lz            = bias$participant$lz[match(ids, bias$participant$id)],
    n_trials      = as.integer(table(model_df$id)[ids]),
    mean_accuracy = vapply(ids, function(i) mean(model_df$correct[model_df$id == i]),
                           numeric(1L)),
    stringsAsFactors = FALSE
  )

  new_inough_signals(
    trial       = trial,
    participant = participant,
    model       = acc$model,
    call        = match.call(),
    formula     = model_formula
  )
}
