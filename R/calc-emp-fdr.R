#' Calculate Empirical Number of False Positives
#'
#' Calculate the mean number of false positive features from
#'   a permutation analysis performed during a stability selection
#'   run. This assumes a permutation set was generated during
#'   stability selection, (i.e. `n_perm > 0`).
#'
#' @param x A `stab_sel` class object generated
#'   via [stability_selection()].
#'
#' @param thresh_seq `numeric(n)`. A sequence in
#'   \verb{[0, 1]} specifying the thresholds.
#'   See `thresh` in [get_stable_features()]
#'
#' @param warn `logical(1)`. Should warnings be triggered if
#'   mean of `< 5` permutations is being returned?
#'
#' @return A named vector indicating the average number
#'   (counts) of false positive features selected at the
#'   various thresholds specified by `thresh_seq`.
#'
#' @author Stu Field, Michael R. Mehan
#'
#' @order 1
#' @seealso [stability_selection()], [get_stable_features()]
#'
#'
#' @examples
#' n_feat      <- 20
#' n_samples   <- 100
#' x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#' colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#' y           <- sample(1:2, n_samples, replace = TRUE)
#'
#' ss <- stability_selection(x, y, n_iter = 25, n_perm = 50,
#'                           r_seed = 101, parallel = TRUE)
#'
#' calc_emp_fdr(ss, seq(0.5, 0.9, 0.1))
#'
#' @export
calc_emp_fdr <- function(x, thresh_seq, warn = TRUE) {

  .check_perm(x)

  if ( missing(thresh_seq) ) {
    stop("Must pass `thresh_seq =` param", call. = FALSE)
  }

  # count FP counts at various thresholds
  perm_path_n <- lapply(x$permpath_list, function(.x) {
      get_stable_features(.x, thresh = thresh_seq, warn = FALSE) |>
        vapply(nrow, 1L, USE.NAMES = TRUE)
    }) |> data.frame() # keeps rownames; `map_df()` does not

  if ( ncol(perm_path_n) < 5L && warn ) {
    signal_info("Returning mean of < 5 permutations ...",
                call. = FALSE)
  }

  apply(perm_path_n, 1L, mean)
}
