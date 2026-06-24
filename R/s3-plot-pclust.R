#' S3 Plot Method for Class "pclust"
#'
#' @rdname progeny_cluster
#'
#' @examples
#' # S3 plot method
#' plot(pclust)
#'
#' plot(clust_iris)
#' @importFrom graphics par legend lines axis plot segments
#' @importFrom ggplot2 aes ggplot geom_point geom_line theme geom_hline
#' @importFrom ggplot2 geom_errorbar scale_color_manual element_text
#' @importFrom withr with_namespace
#'
#' @export
plot.pclust <- function(x, ...) {

  ci95 <- as.data.frame(t(x$ci95_scores))
  pdf <- cbind(progeny_score = x$mean_scores,
               null_score    = x$mean_random_scores,
               ci95) |>
    rn2col("cluster")
  pdf <- pdf |>
    tidyr::pivot_longer(cols = c(null_score, progeny_score),
                        values_to = "score", names_to = "type")
  pdf$type <- factor(pdf$type, levels = c("progeny_score", "null_score"))
  p1 <- pdf |>
    ggplot(aes(x = cluster, y = score, color = type, group = type)) +
    geom_point(size = 3, alpha = 0.75) +
    scale_color_manual(values = unname(col_palette[1:2L]),
                       labels = function(x) gsub("_", " ", x)) +
    geom_line() +
    labs(y = "Stability Score", x = "Clusters (k)",
         color = NULL, title = "Stability Scores") +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.1) +
    theme(plot.title = element_text(hjust = 0.5))

  pdf2 <- cbind(max_dist = x$D_max, gap_dist = x$D_gap) |>
    as.data.frame() |> rn2col("cluster") |>
    tidyr::pivot_longer(cols = c(max_dist, gap_dist),
                        names_to = "type", values_to = "score")

  p2 <- pdf2 |>
    ggplot(aes(x = cluster, y = score, color = type, group = type)) +
    geom_point(size = 3, alpha = 0.75) +
    scale_color_manual(values = unname(col_palette[4:5L]),
                       labels = function(x) gsub("_", " ", x)) +
    geom_line() +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5,
               linetype = "dashed") +
    labs(y = "Stability Score", x = "Clusters (k)",
         color = NULL, title = "Difference Scores") +
    theme(plot.title = element_text(hjust = 0.5))

  p2
  withr::with_namespace("patchwork", p1 + p2)
}
