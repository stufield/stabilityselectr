#' Stable Features at Given Thresholds
#'
#' [get_threshold_features()] calculates the features at or above a given
#'   selection probability threshold. It is a thin wrapper around
#'   [get_stable_features()] returning a list of thresholds corresponding
#'   to `thresh_vec`.
#'
#' @rdname get_stable_features
#'
#' @param thresh_vec `numeric(1)`. A vector of threshold values.
#' @param ... Additional arguments passed to [get_stable_features()],
#'   typically the `add_features =` and `warn =` arguments.
#'
#' @return For [get_threshold_features()], a *list* of
#'   data frames for various thresholds corresponding
#'   to the entries of `thresh_vec`.
#'
#' @examples
#' # l1-logistic
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y           <- sample(1:2, n_samples, replace = TRUE)
#' })
#' stab_sel <- stability_selection(x, y, "l1-logistic", r_seed = 101)
#' get_threshold_features(stab_sel, seq(0.6, 0.9, 0.1))
#' @importFrom stats setNames
#' @export
get_threshold_features <- function(x, thresh_vec, ...) {
  # use this construction for the '...' to pass thru properly ...
  lapply(
    setNames(thresh_vec, paste0("thresh_", thresh_vec)),
    function(.x) get_stable_features(x, thresh = .x, ...)
  )
}
