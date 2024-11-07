#' S3 Plot Method for Class "pclust"
#'
#' @rdname progeny_cluster
#'
#' @param file `character(1)`. A file path to
#'   save the plot to an external file. The extenstion
#'   is used to determine the appropriate graphics device.
#'   Note: be careful of the `...` dimension units differ
#'   among devices!
#'
#' @examples
#' # S3 plot method
#' plot(pclust)
#' plot(clust_iris)
#'
#' # save to file
#' plot(pclust, "cluster.pdf", height = 5, width = 10)
#' plot(pclust, "cluster.png", height = 480, width = 960)
#'
#' @importFrom graphics par legend lines axis plot segments
#' @importFrom withr local_par local_pdf local_png
#' @importFrom withr local_tiff local_bmp local_jpeg local_svg
#'
#' @export
plot.pclust <- function(x, file = NULL, ...) {

  if ( length(ext <- file_ext(file)) > 0L ) {
    local_device <- switch(ext,
                           pdf  = local_pdf,
                           png  = local_png,
                           tiff = local_tiff,
                           bmp  = local_bmp,
                           svg  = local_svg,
                           jpeg = local_jpeg,
                           stop("Invalid file extension: ", value(ext),
                                call. = FALSE))

    local_device(file, ...)
  }

  local_par(
    list(mgp = c(2, 0.75, 0),
         mar = c(3, 4, 3, 1),
         mfrow = 1:2L)
  )

  labels <- colnames(x$scores)
  ci95   <- data.matrix(x$ci95)
  ylims  <- range(ci95, x$mean_random_scores)
  lw     <- 2
  plot(x$mean_scores, lwd = lw, type = "o",
       col = col_palette[[1L]], pch = 1,
       xaxt = "n", ylim = ylims, ylab = "Stability Score",
       main = "Stability Scores", xlab = "Clusters (k)")
  # bars
  segments(x0 = seq_along(x$mean_scores), y0 = ci95[1L, ],
           x1 = seq_along(x$mean_scores), y1 = ci95[2L, ], lwd = 1.5, col = 1)
  # hats
  segments(x0 = rep(seq_along(x$mean_scores) - 0.05, each = 2),
           y0 = as.numeric(ci95),   # recast
           x1 = rep(seq_along(x$mean_scores) + 0.05, each = 2),
           y1 = as.numeric(ci95),   # recast
           lwd = 1.5, col = 1)
  axis(1, at = seq_len(ncol(x$scores)), labels = labels)
  lines(x$mean_random_scores, type = "o", lwd = lw,
        col = col_palette[[2L]])
  legend("bottomright", legend = c("Progeny Scores", "Null Scores"),
         col = unlist(col_palette[1:2L]), pch = 1, lty = 1)
  ylims <- c(min(x$D_max, x$D_gap), max(x$D_max, x$D_gap))
  plot(x$D_max, type = "o", col = col_palette[[4L]], pch = 1,
       main = "Difference Scores", lwd = lw,
       xlab = "Clusters (k)", ylim = ylims, xaxt = "n",
       ylab = "Stability Score")
  axis(1, at = seq_len(ncol(x$scores)), labels = labels)
  lines(x$D_gap, type = "o", lwd = lw, col = col_palette[[5L]], xaxt = "n")
  legend("bottomright", legend = c("Max Dist.", "Gap Dist."),
         col = unlist(col_palette[4:5L]), pch = 1, lty = 1)
}
