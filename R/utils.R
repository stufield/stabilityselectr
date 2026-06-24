
#' Check for `stabpath_matrix` "type"
#'   internal for trapping bad objects.
#' @noRd
is_stabpath_matrix <- function(x) {
  inherits(x, "matrix") &     # it's a matrix class
    !is.null(rownames(x)) &   # has rownames (features)
    ncol(x) > 10L &           # has more than 10 iterations
    sum(is.na(x)) == 0 &      # no NAs
    sum(x < 0) == 0 &         # all in [0,1]
    sum(x > 1) == 0
}

get_cores <- function() {
  sys_ok <- grepl("linux|darwin", R.version$os)   # linux or MacOS?
  chk    <- Sys.getenv("_R_CHECK_LIMIT_CORES_")   # CRAN CHECK variable
  if ( !sys_ok ) {
    signal_info(
      "Windows detected. Continuing with `cores = 1L`."
    )
    1L
  } else if ( !requireNamespace("parallel", quietly = TRUE) ) {
    signal_info(
      "The `parallel` package is not installed ...\n",
      "Continuing with `cores = 1L`."
    )
    1L
  } else if ( nzchar(chk) && chk == "TRUE" ) {
    2L   # use 2 cores in CRAN/Travis/AppVeyor
  } else {
    max(1L, parallel::detectCores() - 2L) # leave 2 cores for system
  }
}


#' Internal for calculating the AUC for all stability paths
#'
#' @param x A `stabpath_matrix` entry of a `stab_path` object.
#'
#' @param values A vector of values to iterate over,
#'   typically `lambda_norm`.
#'
#' @importFrom utils head tail
#' @noRd
calc_path_auc <- function(x, values) {
  stopifnot(
    "Must pass `x` as a `stability path matrix`." = is_stabpath_matrix(x)
  )
  sum_diffs <- head(values, -1L) + tail(values, -1L)
  # flip matrix to perform row-wise operations column-wise
  vapply(
    data.frame(t(x)),
    function(.x) as.double(diff(.x) %*% sum_diffs) / 2,
    0.1
  )
}

log_rfu <- function(x) {
  cls <- class(x)
  cols <- utils::getFromNamespace("get_analytes", "helpr")(x)
  for ( i in cols ) x[[i]] <- log10(x[[i]])
  structure(x, class = cls)
}

col_palette <- c(
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
