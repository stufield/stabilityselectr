#' Calculate Stable Features
#'
#' Returns a data frame object of all features with a maximum
#'   selection probability greater than a minimum threshold.
#'
#' A "stable feature" is defined as a feature with a maximum selection
#'   probability greater than a supplied threshold. This function returns a
#'   `data.frame` of all features that satisfy this criterion along with the
#'   maximum selection probability and the upper bound on the false discover
#'   rate. This false discovery rate bound is only defined for `thresh > 0.5`,
#'   it is otherwise undefined.
#'
#' *IMPORTANT!* If you pass to the `matrix` method, there
#'   is no permutation analysis performed, i.e. no `$EmpFDR` column
#'   in the returned data frame. This calculation takes a long
#'   time and is not always desired, so this method offers the
#'   user a control mechanism for the output behavior.
#'
#' @param x An `stab_sel` class object OR a matrix containing
#'   selection probabilities, i.e. the `stabpath_matrix` entry
#'   of a `stab_sel` class object.
#' @param thresh `numeric(1)` in \verb{[0, 1]}. Minimum selection
#'   probability threshold.
#' @param add_features `character(n)`. A string of additional features
#'   to *force* into the resulting table, irrespective of their
#'   threshold. Used mostly in the S3 plot method to see a given
#'   stability path of a feature not meeting a threshold cutoff.
#'   Must be exact string match.
#' @param warn `logical(1)`. Should warnings be triggered if no stable
#'   features were found at the specified threshold OR if the
#'   FDR upper bound is undefined at thresholds `<= 0.5`?
#'
#' @return A two column data frame containing maximum selection probabilities
#'   and FDR upper bounds when appropriate, see `Details`.
#'
#' @seealso [stability_selection()]
#' @author Stu Field
#' @examples
#' # l1-logistic
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y           <- sample(1:2, n_samples, replace = TRUE)
#' })
#'
#' stab_sel <- stability_selection(x, y)
#'
#' # Stable features at thresh = 0.55
#' get_stable_features(stab_sel, 0.55)
#' @export
get_stable_features <- function(x, thresh, add_features, warn) {
  UseMethod("get_stable_features")
}

#' @noRd
#' @export
get_stable_features.default <- function(x, thresh, add_features, warn) {
  stop(
    "Couldn't find a S3 method for this class object: ",
    value(class(x)),
    call. = FALSE
  )
}

#' @noRd
#' @export
get_stable_features.stab_sel <- function(x, thresh = 0.75, add_features = NULL,
                                         warn = interactive()) {

  # pass thru to S3 matrix method
  df <- get_stable_features(x$stabpath_matrix,
                            thresh       = thresh,
                            add_features = add_features,
                            warn         = warn)

  # the next part is key:
  # it only happens in the `get_stable_features.stab_sel()` method.
  # in the `get_stable_features.matrix()` method we do not calculate empirical FDRs
  # this saves time when not desired and provides a method to control the behavior
  if ( x$perm_data && nrow(df) > 0 ) {
    # unsure which is correct
    denom <- seq_len(nrow(df))
    denom <- nrow(df)    # use this one for now; sgf
    df$EmpFDR <- calc_emp_fdr(x, thresh_seq = df$MaxSelectProb, warn = FALSE) / denom
  }
  df
}

#' If passing the `stabpath_matrix` element of a `stab_sel` object.
#' **This is the main workhorse**
#'
#' @importFrom dplyr filter arrange desc mutate row_number
#' @noRd
#' @export
get_stable_features.matrix <- function(x, thresh = 0.75, add_features = NULL,
                                       warn = interactive()) {

  if ( !is_stabpath_matrix(x) ) {
    stop(
      "The object `x` passed to `get_stable_features()` does not appear to be ",
      "a 'stabpath_matrix'-like object (", value(class(x)), ") ...\n",
      "Are you sure this is a matrix of stability paths?", call. = FALSE
    )
  }

  df <- data.frame(MaxSelectProb = apply(x, 1, max)) |> rn2col("feature")
  df <- dplyr::arrange(df, desc(MaxSelectProb)) |>
    dplyr::filter(MaxSelectProb >= thresh | feature %in% add_features)

  if ( nrow(df) == 0L ) {
    if ( warn ) {
      warning("No stable features at `thresh = ",
              value(round(thresh, 3L)), "`",
              call. = FALSE)
    }
    df <- dplyr::mutate(df, FDRbound = numeric(0))
  } else if ( thresh <= 0.5 ) {
    if ( warn ) {
      warning("FDR upper bound not defined for `thresh <= 0.5`",
              call. = FALSE)
    }
    df <- dplyr::mutate(df, FDRbound = NA_real_)
  } else {
    n_feat <- nrow(x)
    df <- dplyr::mutate(
      df, FDRbound = row_number() / n_feat^2 / (2 * thresh - 1)
    )
  }
  col2rn(df, "feature")
}
