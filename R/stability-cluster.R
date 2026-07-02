#' Perform Stability Clustering
#'
#' Partitioning Around Medoids (PAM) is used both because it uses
#'   a more robust measurement of the cluster centers (medoids) and
#'   because this implementation keeps the cluster labels consistent
#'   across runs, a key feature in calculating the across run
#'   stability. This does not occur using [kmeans()] where the
#'   initial cluster labels are arbitrarily assigned.
#'
#' @family cluster
#'
#' @inheritParams progeny_cluster
#' @inheritParams stability_selection
#'
#' @param k `integer(1)`. The number of clusters.
#'
#' @return A \eqn{n \times (k + 1)} dimensional `tibble` of clustering
#'   probabilities for each `k`, plus a final column
#'   named `prob_k`, which indicates the "most probable"
#'   cluster membership for that sample.
#'
#' @note How do we make sure clusters are indexed the same
#'   as what comes out of k-means? Could be susceptible to index
#'   errors (but seems ok for now).
#'
#' @author Stu Field
#' @seealso [cluster::pam()]
#'
#' @references Hastie, et al. 2009.
#'
#' @examples
#' stab_clust <- stability_cluster(progeny_data, k = 3L, n_iter = 750L)
#'
#' stab_clust$true_cluster <- rep(1:3L, each = 50L)
#'
#' # View stable clusters
#' stab_clust
#'
#' # confusion matrix
#' stab_clust |>
#'   with(table(truth = true_cluster, predicted = prob_k))
#'
#' # View the incorrectly clustered samples (n = 5)
#' filter(stab_clust, prob_k != true_cluster)
#'
#' # Plot Stability Clusters
#' cols <- c("#24135F", "#00A499", "#840B55")
#' withr::with_par(list(mgp = c(2.00, 0.75, 0.0),
#'                      mar = c(3, 4, 3, 1), mfrow = 1:2L), {
#'   plot(progeny_data,
#'        col = cols[stab_clust$true_cluster],
#'        bg  = cols[stab_clust$true_cluster],
#'        pch = stab_clust$true_cluster + 20,
#'        lwd = 1, cex = 1.75, main = "Simulated 3 Cluster Data")
#'   plot(progeny_data,
#'        col = cols[stab_clust$prob_k],
#'        bg  = cols[stab_clust$true_cluster],
#'        pch = stab_clust$true_cluster + 20,
#'        lwd = 2.5, cex = 1.5, main = "Stability Clustering")
#' })
#' @importFrom cluster pam
#' @importFrom tibble as_tibble
#' @importFrom parallel mclapply
#' @importFrom dplyr pick everything mutate
#' @export
stability_cluster <- function(data, k, n_iter = 100L, r_seed = 1234) {

  data <- data.matrix(data, rownames.force = FALSE)
  n    <- nrow(data)

  # set seed for reproducibility
  # Claude says "L'Ecuyer-CMRG" is better cross
  # platform & in parallel processing
  rng <- RNGkind()
  withr::local_seed(r_seed, .rng_kind = "L'Ecuyer-CMRG")
  withr::defer(restore_rng_kind(rng))

  cluster_mat <- parallel::mclapply(1:n_iter, function(i) {
      s <- sample(n, floor(n / 2), replace = FALSE)
      # do not use -> inconsistent labels:
      # clust_idx1 <- stats::kmeans(data[ s,], k)$cluster # nolint: commented_code_linter.
      clust_idx1 <- pam(data[ s, ], k = k)$cluster
      clust_idx2 <- pam(data[-s, ], k = k)$cluster
      out     <- numeric(n)
      out[s]  <- clust_idx1
      out[-s] <- clust_idx2
      out
    }, mc.set.seed = TRUE, mc.cores = get_cores()) |>
    do.call(what = rbind) |> t()

  fix_nm <- function(x, y) paste0(y, x)

  # map over the rows of cluster counts
  apply(cluster_mat, 1, function(.x) {
      c(prop.table(table(.x, useNA = "no")))
    }) |> t() |> as_tibble() |>
    set_Names(fix_nm, "k=") |>
    mutate(prob_k = apply(pick(everything()), 1, which.max))
}
