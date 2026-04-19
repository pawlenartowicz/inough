#' Extract response bias signals
#'
#' Computes lag-1 response repetition indicator (per-trial) and normalized
#' Lempel-Ziv complexity of the response sequence (per-participant).
#'
#' LZ is normalized by a permutation null distribution: 400 random 50/50
#' binary sequences at the modal participant length, giving expected value
#' ~1.0 for random sequences regardless of sample size.
#'
#' @param df Data frame with columns \code{id} and \code{response} (0/1).
#' @return A list with:
#'   \describe{
#'     \item{trial}{Data frame with column \code{resp_lag1} (+1 = same as
#'       previous, -1 = different, 0 = first trial).}
#'     \item{participant}{Data frame with columns \code{id} and \code{lz}
#'       (permutation-normalized LZ76 complexity).}
#'   }
#' @keywords internal
extract_bias <- function(df) {
  ids <- unique(df$id)

  # Null distribution: 400 permutations at modal N, 50/50 split (once for all)
  n_per_id <- vapply(ids, function(id) sum(df$id == id), integer(1L))
  modal_n  <- as.integer(names(which.max(table(n_per_id))))
  null_mean <- mean(replicate(400, lz_complexity(sample(0:1, modal_n,
                                                        replace = TRUE))))

  # Per-trial: lag-1 indicator
  resp_lag1 <- integer(nrow(df))
  # Per-participant: LZ complexity
  part_id <- character(length(ids))
  part_lz <- numeric(length(ids))

  for (i in seq_along(ids)) {
    rows <- which(df$id == ids[i])
    resp <- df$response[rows]

    # Lag-1: first trial = 0, same = +1, different = -1
    lag1 <- integer(length(resp))
    if (length(resp) > 1L) {
      same <- resp[-1L] == resp[-length(resp)]
      lag1[-1L] <- ifelse(same, 1L, -1L)
    }
    resp_lag1[rows] <- lag1

    # LZ complexity, permutation-normalized
    part_id[i] <- as.character(ids[i])
    part_lz[i] <- lz_complexity(resp) / null_mean
  }

  list(
    trial = data.frame(resp_lag1 = resp_lag1),
    participant = data.frame(id = part_id, lz = part_lz,
                             stringsAsFactors = FALSE)
  )
}
