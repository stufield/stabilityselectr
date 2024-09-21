#' Calculate Stable Features
#'
#' Returns a data frame object of all features with a maximum
#' selection probability greater than a minimum threshold.
#'
#' A "stable feature" is defined as a feature with a maximum selection
#' probability greater than a supplied threshold. This function returns a
#' data.frame of all features that satisfy this criterion along with the
#' maximum selection probability and the upper bound on the false discover
#' rate. This false discovery rate bound is only defined for `thresh > 0.5`,
#' it is otherwise undefined.
#'
#' *IMPORTANT!* If you pass to the `matrix` method, there
#'   is no permutation analysis performed, i.e. no `$EmpFDR` column
#'   in the returned data frame. This calculation takes a long
#'   time and is not always desired, so this method offers the
#'   user a control mechanism for the output behavior.
#'
#' @param x An object of class `stab_sel` OR a matrix containing
#'   selection probabilities, i.e. the `stabpath.matrix` entry
#'   of a `stab_sel` class object.
#' @param thresh Numeric in \verb{[0, 1]}. Minimum selection
#'   probability threshold.
#' @param add.features Character. A string of additional features
#'   to *force* into the resulting table, irrespective of their
#'   threshold. Used mostly in the S3 plot method to see a given
#'   stability path of a feature not meeting a threshold cutoff.
#'   Must be exact string match.
#' @param warn Logical. Should warnings be triggered if no stable
#'   features were found at the specified threshold OR if the
#'   FDR upper bound is undefined at thresholds `<= 0.5`.
#' @return A two column data frame containing maximum selection probabilities
#'   and FDR upper bounds when appropriate, see `Details`.
#' @seealso [stabilitySelection()]
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
#' stab_sel <- stabilitySelection(x, y)
#'
#' # Stable features at thresh = 0.55
#' getStableFeatures(stab_sel, 0.55)
#' @export
getStableFeatures <- function(x, thresh, add.features, warn) {
  UseMethod("getStableFeatures")
}

#' S3 method default
#' @noRd
#' @export
getStableFeatures.default <- function(x, thresh, add.features, warn) {
  stop(
    "Couldn't find a S3 method for this class object: ", value(class(x)),
    call. = FALSE
  )
}

#' S3 `getStableFeatures()` method for `stab_sel` class.
#' @noRd
#' @export
getStableFeatures.stab_sel <- function(x, thresh = 0.75, add.features = NULL,
                                       warn = interactive()) {

  # pass thru to S3 matrix method
  df <- getStableFeatures(x$stabpath.matrix,
                          thresh       = thresh,
                          add.features = add.features,
                          warn         = warn)

  # the next part is key:
  # it only happens in the `getStableFeatures.stab_sel()` method.
  # in the `getStableFetures.matrix()` method we do not calculate empirical FDRs
  # this saves time when not desired and provides a method to control the behavior
  if ( x$perm.data && nrow(df) > 0 ) {
    # unsure which is correct
    denom <- seq_len(nrow(df))
    denom <- nrow(df)    # use this one for now; sgf
    df$EmpFDR <- calcEmpFDR(x, thresh.seq = df$MaxSelectProb, warn = FALSE) / denom
  }
  df
}

#' S3 `getStableFeatures()` method for `matrix` class.
#'
#' If passing the `stabpath.matrix` element of a `stab_sel` object.
#' **This is the main workhorse**
#'
#' @importFrom dplyr filter arrange desc mutate row_number
#' @noRd
#' @export
getStableFeatures.matrix <- function(x, thresh = 0.75, add.features = NULL,
                                     warn = interactive()) {

  if ( !is_stabpath.matrix(x) ) {
    stop(
      "The object `x` passed to `getStableFeatures()` does not appear to be ",
      "a 'stabpath.matrix'-like object (", value(class(x)), ") ...\n",
      "Are you sure this is a matrix of stability paths?", call. = FALSE
    )
  }

  df <- data.frame(MaxSelectProb = apply(x, 1, max)) |> rn2col("feature")
  df <- dplyr::arrange(df, desc(MaxSelectProb)) |>
    dplyr::filter(MaxSelectProb >= thresh | feature %in% add.features)

  if ( nrow(df) == 0 ) {
    if ( warn ) {
      warning("No stable features at `thresh = ", value(round(thresh, 3)), "`",
              call. = FALSE)
    }
    df <- dplyr::mutate(df, FDRbound = numeric(0))
  } else if ( thresh <= 0.5 ) {
    if ( warn ) {
      warning("FDR upper bound not defined for `thresh <= 0.5`", call. = FALSE)
    }
    df <- dplyr::mutate(df, FDRbound = NA_real_)
  } else {
    n_feat <- nrow(x)
    df <- dplyr::mutate(df, FDRbound = row_number() / n_feat^2 / (2 * thresh - 1))
  }
  col2rn(df, "feature")
}
