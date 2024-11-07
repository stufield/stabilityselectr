#' Calculate Empirical Number of False Positives
#'
#' Calculate the mean number of false positive features from
#'   a permutation analysis performed during a stability selection
#'   run. This assumes a permutation set was generated during
#'   stability selection, (i.e. `num_perms > 0`).
#'
#' @family empirical FDR
#'
#' @param x A `stab_sel` class object generated
#'   via [stability_selection()].
#' @param thresh_seq `numeric(n)`. A sequence in
#'   \verb{[0, 1]} specifying the thresholds to evaluate.
#' @param warn `logical(1)`. Should warnings be triggered if
#'   mean of `< 5` permutations is being returned?
#'
#' @return A named vector indicating the average number
#'   (counts) of false positive features selected at the
#'   various thresholds specified by `thresh_seq`.
#'
#' @author Stu Field, Michael R. Mehan
#' @seealso [get_threshold_features()]
#'
#' @examples
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y  <- sample(1:2, n_samples, replace = TRUE)
#' })
#' ss <- stability_selection(x, y, "l1-logistic", num_iter = 25,
#'                           num_perms = 25, r_seed = 101, parallel = TRUE)
#' calc_emp_fdr(ss, seq(0.5, 0.9, 0.1))
#' @export
calc_emp_fdr <- function(x, thresh_seq, warn = TRUE) {

  if ( !is_stab_sel(x) ) {
    stop("`x` argument *must* be class `stab_sel`.", call. = FALSE)
  }

  if ( missing(thresh_seq) ) {
    stop("Must pass `thresh_seq =` argument.", call. = FALSE)
  }

  if ( !x$perm_data || is.null(x$permpath_list) ) {
    stop("No permuted data ...\n",
         "Please set `num_perms > 0` in `stability_selection()`",
         call. = FALSE)
  }

  # count FP counts at various thresholds
  perm_path_n <- lapply(x$permpath_list, function(.x) {
      get_threshold_features(.x, thresh_vec = thresh_seq, warn = FALSE) |>
        vapply(nrow, 1L)
    }) |> data.frame() # keeps rownames; `map_df()` does not

  if ( ncol(perm_path_n) < 5L && warn ) {
    signal_info("* Returning mean of < 5 permutations ...",
                call. = FALSE)
  }
  apply(perm_path_n, 1, mean)
}
