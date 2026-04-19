#' Example dual-task data with diverse inattention profiles
#'
#' A minimal anonymized subset of the Dual Task (Gabor orientation under
#' motor interference) dataset, intended for demonstrating the inough
#' pipeline. Twenty participants were sampled to span a range of attention
#' profiles: clean performers, two bail-out cases (one for response
#' stereotypy, one for chance-level accuracy), participants with localized
#' inattention chunks, and participants with extended inattention periods.
#'
#' Participant identifiers have been replaced with arbitrary codes
#' (\code{P01}--\code{P20}) and any session timing information has been
#' removed.
#'
#' @format A data frame with one row per trial and the following columns:
#' \describe{
#'   \item{participant}{Anonymized participant identifier (factor,
#'     \code{P01}--\code{P20}).}
#'   \item{block}{Block index within the session (integer, \code{>= 1};
#'     practice block excluded).}
#'   \item{trial}{Trial index within the participant's session (integer).}
#'   \item{stim}{Stimulus identifier (integer).}
#'   \item{weight}{Stimulus weight / contrast (numeric).}
#'   \item{orient}{Gabor orientation code (integer).}
#'   \item{cue_type}{Cue type code (integer).}
#'   \item{response}{Participant's response (integer, two unique values).}
#'   \item{correct}{Trial accuracy (integer, 0 or 1).}
#' }
#'
#' @source A subset of the Dual Task (\code{s_9}) data collected in the
#'   COST/Kraken consciousness study (Krakow site). Participant IDs have
#'   been re-coded for anonymity.
#'
#' @examples
#' data(task_example)
#' head(task_example)
#'
#' \donttest{
#' signals <- inough_signals(
#'   task_example,
#'   correct  = correct ~ stim + weight + orient + cue_type + block,
#'   response = response ~ response,
#'   id       = "participant"
#' )
#' det <- inough_detect(signals)
#' summary(det)
#' }
"task_example"
