#' Perform Stability Clustering
#'
#' Partitioning Around Medoids (PAM) is used both because
#'   is uses a more robust measurement of the cluster centers (medoids) and
#'   because this implementation keeps the cluster labels consistent across
#'   runs, a key feature in calculating the across run stability. This does
#'   not occur using [kmeans()] where the initial cluster
#'   labels are arbitrarily assigned.
#'
#' @family cluster
#'
#' @inheritParams progeny_cluster
#' @param k `integer(1)`. The number of clusters.
#' @param iter `integer(1)`. The number of random subset iterations to perform.
#' @return A n x (k + 1) dimensional `tibble` of clustering probabilities for
#'   each `k`, plus a final column named `"ProbK"`, which indicates
#'   the "most probable" cluster membership for that sample.
#'
#' @note How do we make sure clusters are indexed the same as what
#'   comes out of k-means? Worried about index errors (but seems ok for now).
#' @author Stu Field
#' @seealso [pam()]
#'
#' @references Hastie, et al. 2009.
#'
#' @examples
#' stab_clust <- withr::with_seed(999, stability_cluster(progeny_data, k = 3, iter = 500))
#' table(actual = rep(1:3L, each = 50L), predicted = stab_clust$ProbK)
#'
#' stab_clust$true_cluster <- rep(1:3L, each = 50L)
#'
#' # View the stable clusters
#' stab_clust
#'
#' # View the incorrectly clustered samples (n = 6)
#' filter(stab_clust, ProbK != true_cluster)
#'
#' # Plot Stability Clusters
#' cols <- c("#24135F", "#00A499", "#840B55")
#' withr::with_par(list(mgp = c(2.00, 0.75, 0.0), mar = c(3, 4, 3, 1), mfrow = 1:2L), {
#'   plot(progeny_data,
#'        col = cols[stab_clust$true_cluster],
#'        bg  = cols[stab_clust$true_cluster],
#'        pch = stab_clust$true_cluster + 20,
#'        lwd = 1, cex = 1.75, main = "Simulated 3 Cluster Data")
#'   plot(progeny_data,
#'        col = cols[stab_clust$ProbK],
#'        bg  = cols[stab_clust$true_cluster],
#'        pch = stab_clust$true_cluster + 20,
#'        lwd = 2.5, cex = 1.5, main = "Stability Clustering")
#' })
#' @importFrom cluster pam
#' @importFrom tibble as_tibble
#' @export
stability_cluster <- function(data, k, iter = 100) {

  data <- data.matrix(data, rownames.force = FALSE)
  n    <- nrow(data)

  cluster_mat <- replicate(iter,  {
      s <- sample(n, floor(n / 2), replace = FALSE)
      # don't use; inconsistent labels
      # clust_idx1 <- stats::kmeans(data[ s,], k)$cluster # nolint: commented_code_linter.
      clust_idx1 <- pam(data[ s, ], k = k)$cluster
      clust_idx2 <- pam(data[-s, ], k = k)$cluster
      out <- numeric(n)
      out[s]  <- clust_idx1
      out[-s] <- clust_idx2
      out
    })

  # map over the rows of custer counts
  res <- apply(cluster_mat, 1, function(.x) {
      prop <- prop.table(table(.x, useNA = "no"))
      c(prop, which.max(prop))
    }) |> t() |>
    data.frame() |> as_tibble()
  names(res) <- gsub("^X", "k=", names(res))
  dplyr::rename(res, ProbK = ncol(res))
}
