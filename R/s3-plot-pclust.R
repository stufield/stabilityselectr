#' S3 Plot Method for Class "pclust"
#'
#' @rdname progenyCluster
#' @param filename Character. A file name and path to save the plot
#' to an external file.
#' @examples
#' # S3 plot method
#' plot(pclust)
#' plot(clustIris)
#'
#' @importFrom graphics par legend lines axis plot segments
#' @importFrom withr local_par defer
#' @importFrom SomaPlotr figure close_figure
#' @export
plot.pclust <- function(x, ..., filename = NULL) {

  figure(filename, width = 10, height = 5)
  local_par(list(mgp = c(2, 0.75, 0), mar = c(3, 4, 3, 1), mfrow = 1:2))
  defer(close_figure(filename))
  labels <- colnames(x$scores)
  ci95   <- data.matrix(x$ci95)
  ylims  <- range(ci95, x$mean.random.scores)
  lw     <- 2
  plot(x$mean.scores, lwd = lw, type = "o",
       col = SomaPlotr::soma_colors[[1L]], pch = 1,
       xaxt = "n", ylim = ylims, ylab = "Stability Score",
       main = "Stability Scores", xlab = "Clusters (k)")
  # bars
  segments(x0 = seq_along(x$mean.scores), y0 = ci95[1L, ],
           x1 = seq_along(x$mean.scores), y1 = ci95[2L, ], lwd = 1.5, col = 1)
  # hats
  segments(x0 = rep(seq_along(x$mean.scores) - 0.05, each = 2),
           y0 = as.numeric(ci95),   # recast
           x1 = rep(seq_along(x$mean.scores) + 0.05, each = 2),
           y1 = as.numeric(ci95),   # recast
           lwd = 1.5, col = 1)
  axis(1, at = seq_len(ncol(x$scores)), labels = labels)
  lines(x$mean.random.scores, type = "o", lwd = lw,
        col = SomaPlotr::soma_colors[[2L]])
  legend("bottomright", legend = c("Progeny Scores", "Null Scores"),
         col = unlist(SomaPlotr::soma_colors[1:2]), pch = 1, lty = 1)
  ylims <- c(min(x$D.max, x$D.gap), max(x$D.max, x$D.gap))
  plot(x$D.max, type = "o", col = SomaPlotr::soma_colors[[4L]], pch = 1,
       main = "Difference Scores", lwd = lw,
       xlab = "Clusters (k)", ylim = ylims, xaxt = "n",
       ylab = "Stability Score")
  axis(1, at = seq_len(ncol(x$scores)), labels = labels)
  lines(x$D.gap, type = "o", lwd = lw, col = SomaPlotr::soma_colors[[5L]], xaxt = "n")
  legend("bottomright", legend = c("Max Dist.", "Gap Dist."),
         col = unlist(SomaPlotr::soma_colors[4:5]), pch = 1, lty = 1)
}
