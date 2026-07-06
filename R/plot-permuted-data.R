
#' @describeIn calc_emp_fdr
#'   Plot the permutation paths for a `stab_sel` object.
#'   These paths are the stability selection paths
#'   of the `n` class scrambled permutations, i.e. the null.
#'
#' @param which `integer(1)`. Which of the null hypothesis
#'   permuted stability paths to plot.
#'
#' @examples
#' # Plot the permuted data individually
#' plot_permuted_data(ss, 3L)   # choose 3rd permutation
#' @importFrom graphics plot
#' @export
plot_permuted_data <- function(x, which = NULL, ...) {
  stopifnot(
    "`x` must have permuted data. Please set `n_perm`." = x$perm_data
  )
  if ( is.null(which) ) {
    L   <- length(x$permpath_list)
    idx <- readline(
             paste(
               "Which of the", L, "permuted null stability",
               "paths do you want to plot (integer)? "
             )
           )
    which <- as.numeric(idx)
    if ( which > L ) {
      stop( "Must choose a value in [1, ", L, "].", call. = FALSE)
    }
  }
  x$lambda          <- x$perm_lambda
  x$stabpath_matrix <- x$permpath_list[[which]]
  plot(x, plot_perm = FALSE, ...)
}
