
#' Check for `stabpath_matrix` "type"
#' Internal for trapping bad objects.
#' @noRd
is_stabpath_matrix <- function(x) {
  inherits(x, "matrix") &     # it's a matrix class
    !is.null(rownames(x)) &   # has rownames (features)
    ncol(x) > 10 &            # has more than 10 iterations
    sum(is.na(x)) == 0 &      # no NAs
    sum(x < 0) == 0 &         # all in [0,1]
    sum(x > 1) == 0
}

#' Internal for calculating the AUC for all stability paths
#'
#' @param x A `stabpath_matrix` entry of a `stab_path` object.
#' @param values A vector of values to iterate over, typically lambda_norm.
#' @importFrom utils head tail
#'
#' @noRd
calc_path_auc <- function(x, values) {
  stopifnot(is_stabpath_matrix(x))
  sum_diffs <- head(values, -1L) + tail(values, -1L)
  # flip matrix to perform row-wise operations column-wise
  vapply(
    data.frame(t(x)),
    function(.x) as.double(diff(.x) %*% sum_diffs) / 2,
    0.1
  )
}

get_analytes <- getFromNamespace("get_analytes", "helpr")

log_rfu <- function(x) {
  cls <- class(x)
  cols <- get_analytes(x)
  for ( i in cols ) x[[i]] <- log10(x[[i]])
  structure(x, class = cls)
}

col_palette <- list(
  purple     = "#24135F",
  lightgreen = "#00A499",
  lightgrey  = "#707372",
  magenta    = "#840B55",
  lightblue  = "#006BA6",
  yellow     = "#D69A2D",
  darkgreen  = "#007A53",
  darkblue   = "#1B365D",
  darkgrey   = "#54585A",
  blue       = "#004C97"
)
