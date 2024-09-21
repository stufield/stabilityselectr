#' Stable Features at Numerous Thresholds
#'
#' [getThresholdFeatures()] calculates the features at or above a given
#'   selection probability threshold. It is a thin wrapper around
#'   [getStableFeatures()] returning a list of thresholds corresponding
#'   to `thresh.vec`.
#'
#' @rdname getStableFeatures
#' @param thresh.vec Numeric vector of threshold values.
#' @param ... Additional arguments passed to [getStableFeatures()],
#'   typically the `add.features =` and `warn =` arguments.
#' @return For [getThresholdFeatures()], a *list* of data frames for various
#'   thresholds corresponding to the entries of `thresh.vec`.
#' @examples
#' # l1-logistic (default)
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y           <- sample(1:2, n_samples, replace = TRUE)
#' })
#' stab_sel <- stabilitySelection(x, y, "l1-logistic", r.seed = 101)
#' getThresholdFeatures(stab_sel, seq(0.6, 0.9, 0.1))
#' @importFrom stats setNames
#' @export
getThresholdFeatures <- function(x, thresh.vec, ...) {
  # use this construction for the '...' to pass thru properly ...
  lapply(
    setNames(thresh.vec, paste0("thresh_", thresh.vec)),
    function(.x) getStableFeatures(x, thresh = .x, ...)
  )
}
