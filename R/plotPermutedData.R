
#' @describeIn plotEmpFDR
#' Plot the permutation paths for an object of class `stab_sel`.
#' These paths are the stability selection paths of the `n` class scrambled
#' permutations, i.e. the null.
#'
#' @param which Integer. Indicating which of the null hypothesis
#'   permuted stability paths to plot.
#' @param ... Additional arguments passed to [plot.stab_sel()].
#' @examples
#' # Plot the permuted data individually
#' plotPermutedData(stab_sel, 3L)   # choose 3rd permutation
#' @importFrom graphics plot
#' @export
plotPermutedData <- function(x, which = NULL, ...) {
  stopifnot(x$perm.data)
  if ( is.null(which) ) {
    L   <- length(x$permpath.list)
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
  x$lambda          <- x$perm.lambda
  x$stabpath.matrix <- x$permpath.list[[which]]
  plot(x, plot.perm = FALSE, ...)
}
