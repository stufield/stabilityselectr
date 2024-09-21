#' Calculate Empirical FDR Break Points
#'
#' Calculates the stability selection threshold, the mean number of
#'   false positive selected features (empirical), and the number of
#'   selected features for specified FDR break points.
#'   Relies on [calcEmpFDR()] to calculate the mean false discovery
#'   based on permutations during the stability selection algorithm.
#'
#' @family empirical FDR
#' @inheritParams calcEmpFDR
#' @param fdr.breaks Numeric. A vector specifying the desired mean number
#'   of empirical false positives at which to calculate various thresholds.
#' @return A list consisting of:
#'   \item{n_selected}{A vector of the number of features selected at each
#'     empirical stability selection threshold}
#'   \item{meanFPs}{A vector of the mean number of false positive selected
#'     features at each empirical stability selection threshold}
#'   \item{breaks}{A `tibble` of containing empirical false positive
#'     summary statistics at each FDR specified break point}
#' @author Stu Field
#' @seealso [getStableFeatures()], [stabilitySelection()]
#' @examples
#' # l1-logistic
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y        <- sample(1:2, n_samples, replace = TRUE)
#' })
#' stab_sel <- stabilitySelection(x, y, "l1-logistic", num.iter = 25,
#'                                num.perms = 25,
#'                                r.seed = 101, parallel = TRUE)
#' calcEmpFDRbreaks(stab_sel)
#' @importFrom dplyr filter bind_rows
#' @importFrom tibble tibble
#' @export
calcEmpFDRbreaks <- function(x, thresh.seq = seq(1, 0.1, by = -0.01),
                             fdr.breaks = c(0.5, 1, 2, 3, 5)) {
  fdr_data <- tibble(
    MeanFPs    = unname(calcEmpFDR(x, thresh.seq, warn = FALSE)),
    n_selected = getThresholdFeatures(x$stabpath.matrix,  # pass matrix here! S3 method
                                      thresh.vec = thresh.seq,
                                      warn = FALSE) |>
      vapply(nrow, 1L, USE.NAMES = FALSE),
    piThresh   = thresh.seq
  )
  break_data <- setNames(fdr.breaks, fdr.breaks) |>
    lapply(function(.x) {
      # all rows where MeanFPs > each break pt; pull first row at break
      dplyr::filter(fdr_data, MeanFPs > .x)[1L, ]
      }) |> bind_rows(.id = "FdrBreaks")
  break_data$FdrBreaks <- as.numeric(break_data$FdrBreaks)
  list(fdr_data = fdr_data, breaks = break_data)
}
