#' Stability Selection
#'
#' Performs stability selection on a set of predictive features and a response
#' variable. Stability selection is performed using a user-specified kernel.
#' For classification problems the `l1-logistic` kernel should be used and the
#' response should be a vector or factor of with two class labels. For lasso
#' the response column should be a numeric vector. For "Cox" the response is a
#' *two column matrix* containing the event time in the first column and the
#' censoring indicator in the second column. Lastly, "ridge" regression was
#' implemented as of Jan 21, 2014.
#
#' The randomized lasso is used if the `alpha` parameter is set to a value
#' less than 1. In a randomized Lasso the model coefficients are randomly
#' re-weighted when calculating the regularization term. This weighting can be
#' performed in two different ways. If `Pw = NA` then these random
#' weights are sampled uniformly between `alpha` and 1. If `Pw` is
#' supplied, then the random weights are chosen to be `alpha` with
#' probability `Pw` and 1 otherwise. The latter choice is used in Theorem
#' 2 in Meinshausen and Buhlmann. Recommended values of `alpha` and
#' `Pw` are \verb{[0.5, 0.2]}.
#'
#' Stability selection can be performed on multiple cores by setting
#' `parallel = TRUE`. This functionality requires
#' [parallel::mclapply()] from the \pkg{parallel} package.
#' This is *not* available for Windows based OS.
#'
#' @param x A numeric \eqn{n x p} matrix of predictive features containing `n`
#' observation rows and `p` feature columns. Alternatively, an object of class
#' `stab_sel` if passing to one of the S3 generic methods.
#' @param y The response variable. If kernel is "l1-logistic" then a
#' vector of binary class labels. If kernel is "Cox" then it is a two column
#' matrix with the event time in the first column and the censoring
#' indicator (1 = event, 0 = censored) in the second column.
#' @param kernel Character. A string describing the underlying model
#' used for selection. Options are:
#' \itemize{
#'   \item "l1-logistic" (default)
#'   \item "lasso"
#'   \item "Cox"
#'   \item "ridge"
#'   \item "multinomial"
#'   \item "pca.sd"
#'   \item "pca.thresh"
#' }
#' @param num.iter Integer. Defining the number of sub-sampling iterations for
#'   the stability selection.
#' @param parallel Logical. Should parallel processing via multiple cores be
#'   implemented? Must be on a Linux platform and have the \pkg{parallel}
#'   package installed. Otherwise defaults to 1 core.
#' @param alpha Numeric. Value defining the weakness parameter for the
#'   randomized regularization. This is the minimum random weight applied to each
#'   beta coefficient in the regularization.
#' @param Pw Numeric. Value defining probability of a weak weight, see
#'   `alpha`. If `Pw = NA` then the coefficient weights are sampled
#'   uniformly from `alpha` to 1.
#' @param standardize Logical. Whether the data should be centered and scaled.
#' @param lambda.min.ratio The minimum value of lambda/max(lambda) to use during
#'   the selection procedure. See [glmnet()].
#' @param num.perms Integer. The number of permutations to use
#'   in calculating the empirical false positive rate.
#' @param beta.threshold Numeric. Floating point value defining selection
#'   levels for `ridge regression`. Since ridge regression will not zero out
#'   coefficients, selection of coefficient curves by selection probability
#'   is not effective. Any variable having a coefficient with absolute value
#'   greater than or equal to `beta.threshold` will be selected.
#' @param elastic.alpha Numeric. Floating point value between 0 and 1. When 0,
#'   the results of [glmnet()] are equivalent to Ridge regression.
#'   When 1, the results are equivalent to Lasso. Any value between 0 and 1
#'   creates a compromise between L1 and L2 penalty.
#' @param lambda.pad The lambda path is padded with high values of lambda in
#'   order to produce a more appealing plot. Occasionally, the degree of
#'   padding needs to be adjusted in order to produce better resolution at low
#'   values of lambda. Typical values for this parameter are 20 (default), 15,
#'   10, or 5.
#' @param impute.outliers Logical. Should statistical outliers
#'   (\eqn{3 * \sigma}) be imputed to approximate a Gaussian distribution
#'   during stability selection?
#' @param impute.n.sigma Numeric. Standard deviation outlier threshold for
#'   imputing outliers if `impute.outliers = TRUE`, ignored otherwise.
#' @param r.seed Seed for the random number generator, allowing for
#'   reproducibility of results.
#' @param ... Additional arguments passed to one of the S3 methods for
#'   objects of class `stab_sel`, generics include:
#'   * [plot.stab_sel()]
#'   * [print.stab_sel()]
#'   * [summary.stab_sel()]
#' @return An object of class `stab_sel`:
#'   \item{stabpath.matrix}{A matrix of \eqn{features x lambda.seq} containing
#'     stability selection probabilities. A row in this matrix corresponds to a
#'     stability selection path for a single feature.}
#'   \item{lambda}{the sequence of lambdas used for regularization. They
#'     correspond to the columns of `stabpath.matrix`.}
#'   \item{alpha}{the weakness parameter provided in the call.}
#'   \item{Pw}{the weak weight probability provided in the call.}
#'   \item{kernel}{the kernel used (e.g. l1-logistic).}
#'   \item{num.iter}{Integer. The number of iterations used in computing
#'     the stability paths.}
#'   \item{standardize}{should the data be standardized prior to analysis?}
#'   \item{lambda.min.ratio}{?}
#'   \item{perm.data}{Logical. Is there permuted data to perform empirical FDR?}
#'   \item{permpath.list}{list containing information to calculated the
#'     permutation paths of the empirical false positive rate.}
#'   \item{perm.lambda}{The lambda used in the permuted lists.}
#'   \item{permpath.max}{max lambda for the permuted lists (I think).}
#'   \item{beta}{A matrix of the betas calculated during the selection process.}
#'   \item{r.seed}{The random seed used.}
#' @author Michael R. Mehan, Stu Field, and Robert Kirk DeLisle
#' @seealso [glmnet()], [getStableFeatures()]
#' @references Meinshausen, N. and Buhlmann, P. (2010), Stability selection.
#'   Journal of the Royal Statistical Society: Series B (Statistical
#'   Methodology), 72: 417-473. doi: 10.1111/j.1467-9868.2010.00740.x
#' @examples
#' # l1-logistic
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samp      <- 100
#'   x           <- matrix(rnorm(n_samp * n_feat), n_samp, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y           <- sample(1:2, n_samp, replace = TRUE)
#'   stab_sel    <- stabilitySelection(x, y, kernel = "l1-logistic", r.seed = 101)
#' })
#'
#' # Cox
#' xcox <- stripMeta(log10(sim_test_data))
#' # Note this works because colnames are already "time" and "status". In real
#' # datasets, need to rename the final matrix as "time" and "status".
#' ycox <- select(sim_test_data, time, status) |> as.matrix()
#' stab_sel_cox <- stabilitySelection(xcox, ycox, kernel = "Cox", r.seed = 101)
#' @importFrom glmnet glmnet
#' @importFrom stats runif setNames
#' @importFrom tibble tibble as_tibble
#' @export
stabilitySelection <- function(x, y = NULL,
                               kernel = c("l1-logistic", "lasso", "ridge",
                                          "Cox", "pca.sd", "pca.thresh",
                                          "multinomial"),
                               num.iter = 100,
                               parallel = FALSE,
                               alpha = 0.8, Pw = 0.5, num.perms = 0,
                               standardize = TRUE,
                               lambda.min.ratio = 0.1,
                               beta.threshold = 0L,
                               elastic.alpha = 1.0,
                               lambda.pad = 20,
                               impute.outliers = FALSE,
                               impute.n.sigma = 3,
                               r.seed = sample(1000, 1), ...) {

  x <- data.matrix(x)      # turn into data matrix

  if ( any(is.na(x)) ) {
    stop(
      "There are NAs detected in the data matrix ...\n",
      "This will cause `glmnet()` to fail ...\n",
      "Please remove the feature or use `splyr::imputeNAs()`",
      call. = FALSE
    )
  }

  kernel <- match.arg(kernel)

  if ( impute.outliers ) {
    x <- apply(x, 2, imputeOutliers, n.sigma = impute.n.sigma)
  }

  # Checks for parallel processing
  if ( parallel ) {
    sys_ok <- grepl("linux|darwin", R.version$os)   # linux or MacOS?
    chk    <- Sys.getenv("_R_CHECK_LIMIT_CORES_")   # CRAN CHECK variable
    if ( !sys_ok ) {   # If Windows, don't do parallel processing
      signal_info(
        "Windows OS detected. Setting 'cores = 1L' for serial processing."
      )
      n_cores <- 1L
    } else if ( !requireNamespace("parallel", quietly = TRUE) ) {
      signal_info(
        "Did not find `parallel` package installed on this Linux system ...\n",
        "Setting 'cores = 1L' for serial processing."
      )
      n_cores <- 1L
    } else if ( nzchar(chk) && chk == "TRUE" ) {
      # use 2 cores in CRAN/Travis/AppVeyor
      n_cores <- 2L
    } else {
      # use almost all cores in devtools::test()
      n_cores <- max(1L, parallel::detectCores() - 2L)
    }
  } else {
    n_cores <- 1L
  }

  signal_done(
    "Using kernel:", value(kernel), "and", value(n_cores),
    ifelse(n_cores > 1, "cores (parallel)", "core (serial)")
  )

  # set random seed here to ensure reproducibility of stability paths
  set.seed(r.seed,                    # parallel reproducibility; RNGkind()
           kind = ifelse(n_cores > 1, "L'Ecuyer", "default"))


  if ( num.perms > 0 ) {
    perm_seq <- seq_len(num.perms)
    perm.y   <- # assigned by if statement
      if ( inherits(y, "matrix") ) {
        lapply(perm_seq,
               function(i) y[sample(1:nrow(y), nrow(y)), ]) # nolint: seq_linter.
      } else if ( is.vector(y) || is.factor(y) ) {
        lapply(perm_seq,
               function(i) y[sample(1:length(y), length(y))]) # nolint: seq_linter.
      } else {
        list()
      }
  }


  # assigned by if/else
  lambda_seq <- if ( kernel == "l1-logistic" ) {
    y <- factor(y)
    p <- ncol(x)
    if ( is.na(Pw) ) {
      W <- stats::runif(p, alpha, 1)
    } else {
      W    <- rep(1, length = p)
      draw <- stats::runif(p, 0, 1)
      W[draw < Pw] <- alpha
    }
    glmnet::glmnet(x, y, family = "binomial",
                   standardize = standardize,
                   lambda.min.ratio = lambda.min.ratio,
                   penalty.factor = W)$lambda

  } else if ( kernel == "multinomial") {
    y <- factor(y)
    p <- ncol(x)
    if ( is.na(Pw) ) {
      W <- stats::runif(p, alpha, 1)
    } else {
      W    <- rep(1, length = p)
      draw <- stats::runif(p, 0, 1)
      W[draw < Pw] <- alpha
    }
    glmnet::glmnet(x, y, family = "multinomial", standardize = standardize,
                   lambda.min.ratio = lambda.min.ratio, penalty.factor = W)$lambda

  } else if ( kernel == "pca.sd" ) {
    if ( floor(alpha) != alpha ) {
      stop(
        "Bad alpha: ", alpha,
        "\nPlease set to number of principal components to consider.",
        call. = FALSE
      )
    }
    max.lambda <- 10
    exp(
      seq(log(max.lambda),                      # log-space
          log(max.lambda * lambda.min.ratio),   # hi resolution at low values
          length.out = 100)
      )                                         # back to linear space

  } else if ( kernel == "pca.thresh" ) {
    if ( floor(alpha) != alpha ) {
      stop(
        "Bad alpha: ", alpha,
        "\nPlease set to number of principal components to consider.",
        call. = FALSE
      )
    }
    max.lambda <- 0.5
    exp(
      seq(log(max.lambda),
        log(max.lambda * lambda.min.ratio),
        length.out = 100)
      )

  } else if ( kernel == "lasso" ) {
    p <- ncol(x)
    if ( is.na(Pw) ) {
      W <- stats::runif(p, alpha, 1)
    } else {
      W    <- rep(1, length = p)
      draw <- stats::runif(p, 0, 1)
      W[draw < Pw] <- 1.0 / alpha    # as per the paper version - RKD 2013 Nov 8
    }
    glmnet::glmnet(x, y, family = "gaussian",
                   standardize = standardize,
                   lambda.min.ratio = lambda.min.ratio,
                   penalty.factor = W,
                   alpha = elastic.alpha)$lambda

  } else if ( kernel == "ridge" ) {
    if ( elastic.alpha != 0 ) {
      stop(
        "Invalid `elastic.alpha =` argument for 'ridge' kernels (",
        elastic.alpha, "). Please set `elastic.alpha = 0`.", call. = FALSE
      )
    }
    p <- ncol(x)
    if ( is.na(Pw) ) {
      W <- stats::runif(p, alpha, 1)
    } else {
      W    <- rep(1, length = p)
      draw <- stats::runif(p, 0, 1)
      W[draw < Pw] <- 1.0 / alpha     # as per the paper version - RKD 2013 Nov 8
    }
    glmnet::glmnet(x, y, family = "gaussian",
                   standardize = standardize,
                   penalty.factor = W, alpha = 0)$lambda

  } else if ( kernel == "Cox" ) {
    nsteps <- 100
    S <- survival::Surv(y[, 1L], y[, 2L])
    y <- cbind(time = S[, 1L], status = S[, 2L])
    glmnet::glmnet(x, y, nlambda = nsteps, family = "cox",
                   standardize = standardize,
                   lambda.min.ratio = lambda.min.ratio)$lambda
  }


  if ( !kernel %in% c("pca.thresh", "pca.sd") ) {
    # pad lambda at the upper end
    seq_shift <- lambda_seq[1L] / lambda_seq[2L]
    pad <- rev(seq_len(lambda.pad)) |>  # default 1:20 & invert 20:1
      vapply(function(.x) lambda_seq[1L] * seq_shift^.x, 0.1)
    lambda_seq <- c(pad, lambda_seq)
  }

  # assigned by if/else
  perm_lambda <-
    if ( length(lambda_seq) > 20L ) {
      lambda_seq[round(seq(1, length(lambda_seq),
                           length.out = round(length(lambda_seq) / 6)))]
    } else {
      lambda_seq
    }

  permpath_list <- NULL

  # THE ACTION!
  if ( n_cores > 1 ) {        # if parallel processing
    stab_time <- system.time(
      stabpath_mat <- parallel::mclapply(1:num.iter, function(i) {
             calcStabilityPaths(x, y,
                                kernel         = kernel,
                                lambda.seq     = lambda_seq,
                                alpha          = alpha,
                                Pw             = Pw,
                                beta.threshold = beta.threshold,
                                elastic.alpha  = elastic.alpha,
                                standardize    = standardize)
        }, mc.cores = n_cores) |> Reduce(f = "+") # add them all up
    )
    signal_done(
      "Stablity path run time:", value(round(stab_time["elapsed"], 4))
    )

    if ( num.perms > 0 ) {
      perm_time <- system.time(
        permpath_list <- lapply(perm_seq, function(.j) {
            parallel::mclapply(1:num.iter, function(i) {
                calcStabilityPaths(x, perm.y[[.j]],
                                   kernel         = kernel,
                                   lambda.seq     = perm_lambda,
                                   alpha          = alpha,
                                   Pw             = Pw,
                                   beta.threshold = beta.threshold,
                                   standardize    = standardize,
                                   elastic.alpha  = elastic.alpha)
              }, mc.cores = n_cores) |> Reduce(f = "+")
          }) |> setNames(sprintf("Perm_%03i", perm_seq))
      )
      signal_done(
        "Perm path run time:", value(round(perm_time["elapsed"], 4))
      )
    }

  } else {
    stabpath_mat <- replicate(num.iter, simplify = FALSE,
                      calcStabilityPaths(x, y,
                        kernel         = kernel,
                        lambda.seq     = lambda_seq,
                        alpha          = alpha,
                        Pw             = Pw,
                        standardize    = standardize,
                        beta.threshold = beta.threshold,
                        elastic.alpha  = elastic.alpha)) |> Reduce(f = "+")
    if ( num.perms > 0 ) {
      permpath_list <- lapply(perm_seq, function(.x) {
          replicate(num.iter, simplify = FALSE,
            calcStabilityPaths(x, perm.y[[.x]],
                               kernel         = kernel,
                               lambda.seq     = perm_lambda,
                               alpha          = alpha,
                               Pw             = Pw,
                               standardize    = standardize,
                               beta.threshold = beta.threshold,
                               elastic.alpha  = elastic.alpha)) |> Reduce(f = "+")
        }) |>
        setNames(sprintf("Perm_%03i", perm_seq))
    }
  }

  stabpath_mat <- stabpath_mat / num.iter / 2
  rownames(stabpath_mat) <- colnames(x)

  if ( num.perms > 0 ) {
    permpath_list <- lapply(permpath_list, function(.mat) {
      structure(.mat / num.iter / 2, dimnames = list(colnames(x), NULL))
    })
  }

  if ( is.null(permpath_list) ) {
    perm_max_mat <- tibble(Feature = character(0))
  } else {
    perm_max_mat <- lapply(permpath_list, function(.x) apply(.x, 1, max)) |>
      data.frame() |>        # preserve feature rownames
      rn2col("Feature") |>   # move features to column
      as_tibble()            # df => tibble
  }

  # have to get the betas for the appropriate lambda sequence
  beta <- switch(kernel,
         "l1-logistic" = glmnet::glmnet(x, y, family = "binomial",
                                        standardize = standardize,
                                        lambda = lambda_seq,
                                        penalty.factor = W)$beta,
         "multinomial" = glmnet::glmnet(x, y, family = "multinomial",
                                        standardize = standardize,
                                        lambda = lambda_seq,
                                        penalty.factor = W)$beta,
         "lasso"       = glmnet::glmnet(x, y, family = "gaussian",
                                        standardize = standardize,
                                        lambda = lambda_seq,
                                        penalty.factor = W,
                                        alpha = elastic.alpha)$beta,
         "ridge"       = glmnet::glmnet(x, y, family = "gaussian",
                                        standardize = standardize,
                                        penalty.factor = W, alpha = 0)$beta,
         "Cox"         = glmnet::glmnet(x, y, nlambda = nsteps,
                                        family = "cox",
                                        standardize = standardize,
                                        lambda = lambda_seq)$beta,
         "pca.sd"      = NULL,
         "pca.thresh"  = NULL)

  list(stabpath.matrix  = stabpath_mat,
       lambda           = lambda_seq,
       alpha            = alpha,
       Pw               = Pw,
       kernel           = kernel,
       num.iter         = num.iter,
       standardize      = standardize,
       impute.outliers  = impute.outliers,
       lambda.min.ratio = lambda.min.ratio,
       perm.data        = length(permpath_list) > 0L & num.perms > 0,
       permpath.list    = permpath_list,
       perm.lambda      = perm_lambda,
       permpath.max     = perm_max_mat,
       beta             = beta,
       r.seed           = r.seed) |> addClass("stab_sel")
}


#' Test for object type "stab_sel"
#'
#' The [is.stab_sel()] function checks whether
#' an object is of class `stab_sel`. See [inherits()].
#' @rdname stabilitySelection
#' @return The `is.stab_sel` function returns a logical boolean.
#' @examples
#' # Test for class `stab_sel`
#' is.stab_sel(stab_sel)
#'
#' @export
is.stab_sel <- function(x) inherits(x, "stab_sel")


#' @describeIn stabilitySelection
#' The S3 `print` method for class `stab_sel`.
#'
#' @return The S3 print method returns:
#' \item{Stability Selection Kernel}{The kernel used in the stability selection
#'       algorithm.}
#' \item{Weakness}{The weakness used (`alpha` argument).}
#' \item{Weakness Probability}{The probability of the weakness being applied
#'       (`Pw =` argument).}
#' \item{Number of Iterations}{Number of iterations in the selection
#'       (`num.iter =` argument).}
#' \item{Standardized}{Was the data standardized prior to stability selection?}
#' \item{Imputed Outliers}{Were statistical outliers imputed to the Gaussian
#'       approximation prior to stability selection?}
#' \item{Lambda Max}{The maximum lambda (tuning parameter) used.}
#' \item{Lambda Max}{The maximum lambda (tuning parameter) used.}
#' \item{Random Seed}{The seed passed to the random number generator
#'       for subset selection.}
#' @examples
#' # S3 print method
#' stab_sel
#'
#' @export
print.stab_sel <- function(x, ...) {
  writeLines(
    signal_rule(
      sprintf("Stability Selection (%s: %s)",
              add_color("Kernel", "cyan"),
              add_color(x$kernel, "green")),
      line_col = "blue", lty = "double")
  )
  key <- c(
    "Weakness (alpha)",
    "Weakness Probability (Pw)",
    "Number of Iterations",
    "Standardized ",
    "Imputed Outliers",
    "Lambda Max",
    "Lambda Min Ratio",
    "Permuted Data",
    "Random Seed"
    ) |> pad(27)
  value <- list(
    round(x$alpha, 2L),
    round(x$Pw, 2L),
    round(x$num.iter),
    ifelse(x$standardize, "Yes", "No"),
    ifelse(x$impute.outliers, "Yes", "No"),
    round(x$lambda[1L], 4L),
    round(x$lambda.min.ratio, 1L),
    ifelse(x$perm.data, "Yes", "No"),
    x$r.seed
  )
  liter(key, value, function(.x, .y) {
    writeLines(paste(add_color(symbl$bullet, "red"), .x, value(.y)))
  })
  writeLines(signal_rule(line_col = "green", lty = "double"))
  invisible(x)
}


#' @describeIn stabilitySelection
#' The S3 `summary` method for class `stab_sel`.
#'
#' @param object An object of class `stab_sel`.
#' @param thresh A numeric minimum selection probability threshold.
#' This value can also be a vector of values in \verb{[0, 1]},
#' but ideally greater than 0.50.
#' @note Additional features can be passed as strings to the summary method
#' via the `add.features` argument.
#' @examples
#' # S3 summary method
#' summary(stab_sel, thresh = 0.6)
#' summary(stab_sel, thresh = 0.8, add.features = "feat_c")   # force feat_c into table
#' @importFrom dplyr select everything matches
#' @export
summary.stab_sel <- function(object, ..., thresh) {
  if ( missing(thresh) ) {
    stop(
      "Must pass a stability threshold to S3 summary method. ",
      "Please pass `thresh =`.", call. = FALSE
    )
  }

  lambda_norm <- object$lambda / max(object$lambda)
  df          <- getStableFeatures(object, thresh = thresh, ...)

  if ( nrow(df) > 0 ) {
    df$AUC <- object$stabpath.matrix[rownames(df), ] |>  # reordering features
      calc_path_auc(values = lambda_norm)
  }

  dplyr::select(df, -matches("FDR"), everything()) |>
    rn2col("feature") |>
    as_tibble()
}


#' @describeIn stabilitySelection
#' The S3 `plot` method plots the selection paths for the features. This
#' plot closely resembles a lasso coefficient plot with the regularization
#' parameter (lambda) plotted on x-axis and the feature selection probability
#' (rather than the model coefficient) is plotted on the y-axis.
#'
#' Plots the regularization parameter (lambda) on the x-axis and the selection
#' probability on the y-axis. The regularization parameter is plotted as
#' lambda/max(lambda) so that it is in the range from 1 to 0. The selection
#' probability corresponds to the number of times a particular marker was
#' chosen at a given value of lambda. Each line in the plot is a marker and
#' represents the stability selection path over the range of regularization
#' parameter. All features that have a maximum selection probability greater
#' than `thresh` (shown as a dotted horizontal line) are colored and labeled
#' and the remaining features are colored gray and unlabeled. Additionally, you
#' can provide a set of custom labels that will be colored and labeled
#' regardless of their max selection probability. Each feature is labeled with
#' a capital letter and the full name of the feature is indicated in the legend
#' along with the AUC for its curve in parentheses.
#'
#' @param custom.labels a character vector of additional features to label in
#' the plot, see Details.
#' @param main optional title for the plot (defaults depend on the kernel
#' used)
#' @param sort.by.AUC Logical. If `TRUE`, entries in the legend will be sorted by their
#' curve AUC values which are in parentheses following the variable name in the
#' legend.
#' @param ln.cols A vector of colors to be used as line colors in plotting.
#' Colors are recycled as necessary. See [SomaPlotr::soma_colors].
#' @param add.perm Logical. Should empirical false discovery lines
#' from the null permutation be added to the plot (if permutation was
#' performed)? This can be time consuming depending on the number of
#' permutations performed, so the default is `FALSE`.
#' @param emp.thresh a vector describing the empirical threshold values
#' to be used (`default = seq(1, 0.1, by = 0.01)`).
#' @return The S3 `plot` method returns a data frame of all features
#' labeled in the plot, either by having a maximum selection probability
#' greater than `thresh` or because they were provided as
#' `custom.labels`.
#' @examples
#' # S3 plot method
#' plot(stab_sel, thresh = 0.8)
#'
#' @importFrom utils head
#' @importFrom dplyr case_when arrange select mutate desc left_join
#' @importFrom SomaPlotr theme_soma
#' @export
plot.stab_sel <- function(x, thresh = 0.60,
                          custom.labels = NULL,
                          main = NULL,
                          sort.by.AUC = TRUE,
                          ln.cols = unlist(SomaPlotr::soma_colors),
                          add.perm = FALSE,
                          emp.thresh = seq(1, 0.1, by = -0.01), ...) {

  stopifnot(is_stabpath.matrix(x$stabpath.matrix))   # sanity catch
  lambda_norm <- x$lambda / max(x$lambda)

  if ( is.null(main) ) {
    main <- sprintf("Stability Paths (Weakness = %0.2f, Pw = %0.2f)",
                    x$alpha, x$Pw)
  }

  custom_features <- intersect(custom.labels, rownames(x$stabpath.matrix))
  summary_df      <- summary(x, thresh    = thresh,
                             add.features = custom_features,
                             warn         = FALSE)

  if ( sort.by.AUC ) {
    summary_df <- dplyr::arrange(summary_df, dplyr::desc(AUC))
  }

  L           <- nrow(summary_df)
  legend_cols <- dplyr::case_when(L > 12 & L < 24 ~ 2,
                                  L > 24 & L < 36 ~ 3,
                                  L > 36 ~ 4,
                                  TRUE ~ 1)

  # pre-process the x-axis values for dplyr::left_join to the prob. paths
  # The names will be X1 ... Xn; following conversion to tbl_df
  x_df <- data.frame(column      = paste0("X", seq_len(ncol(x$stabpath.matrix))),
                     lambda_norm = lambda_norm)

  auc_df <- dplyr::select(summary_df, feature, AUC)

  # need to have special getSeqId version here so that non-seqId
  # containing features (e.g. clinical variables) do not return NA
  .get_seq <- function(.x) sub("\\.", "-", sub("^seq\\.", "", .x))

  # get order of feature labels
  label_order <- sprintf("%s (%0.2f)", .get_seq(auc_df$feature), auc_df$AUC)

  line_cols <- rep(unname(ln.cols), times = 5)[1:L]
  names(line_cols) <- label_order

  if ( L < nrow(x$stabpath.matrix) ) {  # if un-selected features, append color map
    line_cols <- c(line_cols, "Unselected (NA)" = "grey80")
  }

  data <- data.frame(x$stabpath.matrix) |>
    rn2col("feature") |>
    dplyr::left_join(auc_df, by = "feature") |>
    dplyr::mutate(
      feat_sel = dplyr::case_when(
        feature %in% summary_df$feature ~ .get_seq(feature),
        TRUE ~ "Unselected"),
      label = factor(sprintf("%s (%0.2f)", feat_sel, AUC), levels = names(line_cols))
    ) |>
    tidyr::pivot_longer(c(-feature, -feat_sel, -AUC, -label),
                        names_to = "column", values_to = "prob") |>
    dplyr::left_join(x_df, by = "column")

  # sanity check for safety
  stopifnot(all(names(line_cols) %in% unique(data$label)))

  p <- data |>
    ggplot2::ggplot(ggplot2::aes(x = lambda_norm, y = prob)) +
    ggplot2::geom_line(ggplot2::aes(group = feature, colour = label)) +
    ggplot2::scale_colour_manual(values = line_cols) +
    ggplot2::scale_x_reverse(expand = c(0, 0)) +
    ggplot2::geom_hline(yintercept = thresh, colour = "black",
                        linetype = "solid") +
    ggplot2::guides(colour = ggplot2::guide_legend(title = "Feature",
                                                   ncol = legend_cols)) +
    ggplot2::labs(x = "lambda / max(lambda)",
                  y = "Selection Probability",
                  title = main) +
    theme_soma(legend.position = "right") +
    NULL

  if ( add.perm ) {
    if ( !x$perm.data ) {
      warning(
        "No permutation data present in `stab_sel` object. ",
        "Skipping permutation cutoffs.", call. = FALSE
      )
    } else {
      emp_breaks_list <- calcEmpFDRbreaks(x, emp.thresh)
      emp_breaks_list$breaks$line <- "dotted"
      break_cols <- unlist(SomaPlotr::soma_colors, use.names = FALSE) |> head(5L)
      p <- p +
        ggplot2::geom_hline(
          data = emp_breaks_list$breaks,
          mapping = ggplot2::aes(yintercept = piThresh,
                                 linetype = factor(MeanFPs)),
          colour = break_cols
          ) +
        ggplot2::scale_linetype_manual(
          name = "MeanFPs",
          values = rep("dashed", length(break_cols))
          ) +
        ggplot2::guides(linetype = ggplot2::guide_legend(
          override.aes = list(colour = break_cols)
            )
          )
    }
  }
  return(p)
}
