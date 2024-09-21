#' Calculate Empirical Number of False Positives
#'
#' Calculate the mean number of false positive features from
#'   a permutation analysis performed during a stability selection
#'   run. This assumes a permutation set was generated during
#'   stability selection, (i.e. `num.perms > 0`).
#'
#' @family empirical FDR
#' @param x An object of class `stab_sel` generated via [stabilitySelection()].
#' @param thresh.seq  Numeric sequence in \verb{[0, 1]} specifying
#'   the thresholds to evaluate.
#' @param warn Logical. Should warnings be triggered if
#'   mean of `< 5` permutations is being returned?
#' @return A named vector indicating the average number (counts) of false
#'   positive features selected at the various thresholds specified
#'   by `thresh.seq`.
#' @author Stu Field, Michael R. Mehan
#' @seealso [getThresholdFeatures()]
#' @examples
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y  <- sample(1:2, n_samples, replace = TRUE)
#' })
#' ss <- stabilitySelection(x, y, "l1-logistic", num.iter = 25,
#'                          num.perms = 25, r.seed = 101, parallel = TRUE)
#' calcEmpFDR(ss, seq(0.5, 0.9, 0.1))
#' @export
calcEmpFDR <- function(x, thresh.seq, warn = TRUE) {

  if ( !is.stab_sel(x) ) {
    stop("`x` argument *must* be class `stab_sel`.", call. = FALSE)
  }

  if ( missing(thresh.seq) ) {
    stop("Must pass `thresh.seq =` argument.", call. = FALSE)
  }

  if ( !x$perm.data || is.null(x$permpath.list) ) {
    stop("No permuted data ...\n",
         "Please set `num.perms > 0` in `stabilitySelection()`",
         call. = FALSE)
  }

  # count FP counts at various thresholds
  perm_path_n <- lapply(x$permpath.list, function(.x) {
      getThresholdFeatures(.x, thresh.vec = thresh.seq, warn = FALSE) |>
        vapply(nrow, 1L)
    }) |> data.frame() # keeps rownames; `map_df()` does not

  if ( ncol(perm_path_n) < 5 && warn ) {
    warning("* Returning mean of < 5 permutations ...", call. = FALSE)
  }
  apply(perm_path_n, 1, mean)
}
