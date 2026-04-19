#' Normalized Lempel-Ziv Complexity (LZ76)
#'
#' Computes the normalized LZ76 complexity of a binary sequence. Higher values
#' indicate more random/complex sequences; lower values indicate more
#' predictable/repetitive patterns.
#'
#' @param x Integer vector of 0s and 1s.
#' @return Numeric scalar in \[0, 1\]. Normalized complexity where 1 = maximally
#'   complex (random) and values near 0 = highly predictable.
#' @export
#' @examples
#' lz_complexity(c(0, 0, 0, 0, 0))       # low
#' lz_complexity(c(0, 1, 0, 1, 0, 1))     # low-medium
#' lz_complexity(sample(0:1, 100, TRUE))   # near 1
lz_complexity <- function(x) {
  if (!is.numeric(x) && !is.integer(x)) {
    stop("x must be a numeric/integer vector of 0s and 1s")
  }
  if (length(x) == 0L) {
    stop("x must have length > 0")
  }
  if (!all(x %in% c(0L, 1L))) {
    stop("x must contain only 0s and 1s")
  }

  n <- length(x)
  if (n == 1L) return(1.0)

  # LZ76 sequential decomposition (exhaustive parsing)
  # Parse x into words w1, w2, ... where each word is the shortest string

  # not found as a substring of everything preceding its last character.

  c_n <- 1L  # first symbol is the first word
  p <- 2L    # start position of current word

  while (p <= n) {
    l <- 1L
    found <- TRUE

    while (found && (p + l - 1L) <= n) {
      # Search for x[p:(p+l-1)] in x[1:(p+l-2)] (exhaustive history)
      search_end <- p + l - 2L
      found <- FALSE

      if (search_end >= l) {
        for (j in seq_len(search_end - l + 1L)) {
          if (all(x[j:(j + l - 1L)] == x[p:(p + l - 1L)])) {
            found <- TRUE
            break
          }
        }
      }

      if (found) l <- l + 1L
    }

    c_n <- c_n + 1L
    p <- p + l
  }

  # Entropy-corrected normalization: c(n) / (n / (H * log2(n)))
  # H = binary entropy of empirical response probability.
  # Reduces to c(n)*log2(n)/n for balanced sequences (p=0.5, H=1).
  # Corrects for response bias without permutations.
  p <- mean(x)
  if (p == 0 || p == 1) return(0)
  h <- -p * log2(p) - (1 - p) * log2(1 - p)
  c_n / (n / (h * log2(n)))
}
