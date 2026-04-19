#' Plot per-participant diagnostic panels
#'
#' Produces a stacked 4-panel visualization: accuracy strip, accuracy
#' residuals, lag-1 response, and dual-track z-scores. Flagged chunks
#' are highlighted as red shaded regions.
#'
#' @param x An \code{inough_detected} object.
#' @param id Character scalar — participant ID to plot.
#' @param ... Ignored.
#' @return A \code{patchwork} object (invisibly).
#' @importFrom rlang .data
#' @export
plot.inough_detected <- function(x, id, ...) {
  if (missing(id)) {
    stop("Specify a participant id. Available: ",
         paste(utils::head(x$signals$participant$id, 20), collapse = ", "),
         if (nrow(x$signals$participant) > 20L) ", ..." else "")
  }
  pid <- as.character(id)
  if (!(pid %in% x$signals$participant$id)) {
    stop("id '", pid, "' not found")
  }

  trial <- x$signals$trial[x$signals$trial$id == pid, ]
  part  <- x$signals$participant[x$signals$participant$id == pid, ]
  tidx  <- trial$trial_idx
  is_bailout <- pid %in% x$bailout$id

  sig_chunks <- x$results[x$results$id == pid & x$results$significant, ,
                          drop = FALSE]

  # --- Header tiles ---
  make_tile <- function(label, value, fill = "#e8e8e8") {
    ggplot2::ggplot() +
      ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                        fill = fill, colour = "#cccccc") +
      ggplot2::annotate("text", x = 0.5, y = 0.62, size = 4,
                        colour = "#666666", label = label) +
      ggplot2::annotate("text", x = 0.5, y = 0.3, size = 7,
                        fontface = "bold", label = value) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::theme_void()
  }

  header <- patchwork::wrap_plots(
    make_tile("LZ",       sprintf("%.3f", part$lz),            "#f5e6d8"),
    make_tile("accuracy", sprintf("%.3f", part$mean_accuracy), "#e3edf7"),
    make_tile("chunks",   as.character(nrow(sig_chunks)),       "#e7f0e3"),
    ncol = 3
  ) + patchwork::plot_annotation(
    title = if (is_bailout) paste("id:", pid, "[BAILED OUT]") else paste("id:", pid),
    theme = ggplot2::theme(plot.title = ggplot2::element_text(
      size = 14, face = "bold", hjust = 0.5,
      colour = if (is_bailout) "red" else "black"
    ))
  )

  # Banner for bail-out (shown above panels, not instead of them)
  banner <- if (is_bailout) {
    reason <- x$bailout$reason[x$bailout$id == pid]
    ggplot2::ggplot() +
      ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                        fill = "#ffcccc") +
      ggplot2::annotate("text", x = 0.5, y = 0.5, size = 5,
                        fontface = "bold", colour = "red",
                        label = paste("Bail-out reason:", reason)) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      ggplot2::theme_void()
  } else NULL

  # --- Compute response track ---
  w         <- x$control$window_size
  positions <- seq_len(2L * w + 1L) - (w + 1L)
  weights   <- switch(x$control$window_weight,
    uniform    = rep(1, 2L * w + 1L),
    triangular = (w + 1L) - abs(positions)
  )
  roll_acc  <- rolling_wmean(trial$acc_resid, weights)
  roll_resp <- rolling_wmean(trial$resp_lag1, weights)
  abs_resp  <- abs(roll_resp)

  thresh <- if (identical(x$control$screening_threshold, "lz")) {
    part$lz
  } else {
    x$control$screening_threshold
  }

  # --- Chunk rectangle layers ---
  has_spurious <- "spurious_start" %in% names(sig_chunks)

  # Main tested region
  chunk_rects <- if (nrow(sig_chunks) > 0L) {
    test_s <- if (has_spurious) sig_chunks$test_start else sig_chunks$start
    data.frame(xmin = test_s - 0.5, xmax = sig_chunks$end + 0.5)
  } else NULL

  # Spurious-trimmed region (lighter shading)
  spur_rects <- if (has_spurious && any(sig_chunks$spurious_start)) {
    sp <- sig_chunks[sig_chunks$spurious_start, , drop = FALSE]
    data.frame(xmin = sp$start - 0.5, xmax = sp$test_start - 0.5)
  } else NULL

  add_chunks <- function(p) {
    if (!is.null(spur_rects)) {
      p <- p + ggplot2::annotate("rect",
        xmin = spur_rects$xmin, xmax = spur_rects$xmax,
        ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12)
    }
    if (!is.null(chunk_rects)) {
      p <- p + ggplot2::annotate("rect",
        xmin = chunk_rects$xmin, xmax = chunk_rects$xmax,
        ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.15)
    }
    p
  }

  strip_x <- ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                             axis.text.x  = ggplot2::element_blank())

  # Panel 1: Accuracy strip
  p1 <- add_chunks(
    ggplot2::ggplot(data.frame(x = tidx, y = trial$correct),
                    ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$x, y = 0,
                                          yend = .data$y), linewidth = 0.3) +
      ggplot2::scale_y_continuous(breaks = c(0, 1)) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(y = "correct") +
      ggplot2::theme_minimal() + strip_x
  )

  # Panel 2: Accuracy residuals + rolling mean
  p2 <- add_chunks(
    ggplot2::ggplot(data.frame(x = tidx, y = trial$acc_resid, r = roll_acc),
                    ggplot2::aes(x = .data$x)) +
      ggplot2::geom_point(ggplot2::aes(y = .data$y), size = 0.3, alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$r), colour = "steelblue",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      ggplot2::labs(y = "acc_resid") +
      ggplot2::theme_minimal() + strip_x
  )

  # Panel 3: Lag-1 response + rolling mean
  p3 <- add_chunks(
    ggplot2::ggplot(data.frame(x = tidx, y = trial$resp_lag1, r = roll_resp),
                    ggplot2::aes(x = .data$x)) +
      ggplot2::geom_point(ggplot2::aes(y = .data$y), size = 0.3, alpha = 0.3) +
      ggplot2::geom_line(ggplot2::aes(y = .data$r), colour = "steelblue",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      ggplot2::labs(y = "resp_lag1") +
      ggplot2::theme_minimal() + strip_x
  )

  # Panel 4: |roll_resp| with threshold
  p4 <- add_chunks(
    ggplot2::ggplot(data.frame(x = tidx, y = abs_resp, r = roll_resp),
                    ggplot2::aes(x = .data$x)) +
      ggplot2::geom_line(ggplot2::aes(y = .data$y), colour = "#d95f02",
                          linewidth = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = .data$r), colour = "#7570b3",
                          linewidth = 0.3, alpha = 0.5) +
      ggplot2::geom_hline(yintercept = thresh, linetype = "dashed",
                           colour = "red", alpha = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
      ggplot2::labs(x = "trial", y = "|roll_resp|") +
      ggplot2::theme_minimal()
  )

  panels <- if (is_bailout) {
    list(header, banner, p1, p2, p3, p4)
  } else {
    list(header, p1, p2, p3, p4)
  }
  heights <- if (is_bailout) c(0.6, 0.3, 0.8, 1.2, 1.2, 1.2) else c(0.6, 0.8, 1.2, 1.2, 1.2)

  combined <- patchwork::wrap_plots(panels, ncol = 1, heights = heights)

  print(combined)
  invisible(combined)
}
