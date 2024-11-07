#' Calculate Empirical FDR Break Points
#'
#' Calculates the stability selection threshold, the mean number of
#'   false positive selected features (empirical), and the number of
#'   selected features for specified FDR break points.
#'   Relies on [calc_emp_fdr()] to calculate the mean false discovery
#'   based on permutations during the stability selection algorithm.
#'
#' @family empirical FDR
#' @inheritParams calc_emp_fdr
#'
#' @param fdr_breaks `numeric(n)`. A vector specifying the desired mean number
#'   of empirical false positives at which to calculate various thresholds.
#' @return A list consisting of:
#'   \item{n_selected}{A vector of the number of features selected at each
#'     empirical stability selection threshold}
#'   \item{meanFPs}{A vector of the mean number of false positive selected
#'     features at each empirical stability selection threshold}
#'   \item{breaks}{A `tibble` of containing empirical false positive
#'     summary statistics at each FDR specified break point}
#'
#' @author Stu Field
#' @seealso [get_stable_features()], [stability_selection()]
#'
#' @examples
#' # l1-logistic
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y        <- sample(1:2, n_samples, replace = TRUE)
#' })
#' stab_sel <- stability_selection(x, y, "l1-logistic", num_iter = 25,
#'                                 num_perms = 25,
#'                                 r_seed = 101, parallel = TRUE)
#' calc_emp_fdr_breaks(stab_sel)
#' @importFrom dplyr filter bind_rows
#' @importFrom tibble tibble
#' @export
calc_emp_fdr_breaks <- function(x, thresh_seq = seq(1, 0.1, by = -0.01),
                                fdr_breaks = c(0.5, 1, 2, 3, 5)) {
  fdr_data <- tibble(
    MeanFPs    = unname(calc_emp_fdr(x, thresh_seq, warn = FALSE)),
    n_selected = get_threshold_features(x$stabpath_matrix,  # pass matrix here! S3 method
                                        thresh_vec = thresh_seq,
                                        warn = FALSE) |>
      vapply(nrow, 1L, USE.NAMES = FALSE),
    piThresh   = thresh_seq
  )
  break_data <- setNames(fdr_breaks, fdr_breaks) |>
    lapply(function(.x) {
      # all rows where MeanFPs > each break pt; pull first row at break
      dplyr::filter(fdr_data, MeanFPs > .x)[1L, ]
      }) |> bind_rows(.id = "FDR_breaks")
  break_data$FDR_breaks <- as.numeric(break_data$FDR_breaks)
  list(fdr_data = fdr_data, breaks = break_data)
}
