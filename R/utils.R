
.check_perm <- function(x) {
  stopifnot(
    "`x` must be a `stab_sel` object."         = is_stab_sel(x),
    "`x` must be creted with `n_perm = TRUE`." = x$perm_data,
    "`x` is missing permuted (null) data."     = !is.null(x$permpath_list)
  )
}

#' @importFrom stats runif
#' @noRd
.calcW <- function(p, alpha, Pw, kernel) {
  if ( is.na(Pw) && kernel != "cox" ) {
    W <- runif(p, alpha, 1.0)
  } else {
    W    <- rep(1.0, length = p)
    draw <- runif(p)
    # as per the paper (RKD: 2013-11-08)
    alpha_star <- ifelse(kernel %in% c("ridge", "lasso"), 1.0 / alpha, alpha)
    W[draw < Pw] <- alpha_star
  }
  W
}

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
  } else if ( nzchar(chk) && chk == "TRUE" ) {
    2L   # use 2 cores in CRAN/Travis/AppVeyor
  } else if ( identical(Sys.getenv("TESTTHAT"), "true") ) {
    2L   # use 2 cores if testing
  } else {
    max(1L, parallel::detectCores() - 2L) # leave 2 cores for system
  }
}

restore_rng_kind <- getFromNamespace("restore_rng_kind", "withr")

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
