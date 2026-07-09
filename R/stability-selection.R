#' Stability Selection
#'
#' Performs stability selection on a set of predictive features
#'   and a response variable. Stability selection is performed using
#'   a user-specified kernel. For classification problems the
#'   `binomial` kernel should be used and the response should be
#'   a vector or factor of with two class labels. For lasso the
#'   response column should be a numeric vector. For "cox" the
#'   response is a *two column matrix* containing the event time in
#'   the first column and the censoring indicator in the second column.
#
#' The randomized lasso is used if the `alpha` parameter is set to
#'   a value less than 1. In a randomized Lasso the model
#'   coefficients are randomly re-weighted when calculating the
#'   regularization term. This weighting can be performed in two
#'   different ways. If `Pw = NA` then these random weights are
#'   sampled uniformly between `alpha` and `1.0`. If `Pw` is supplied,
#'   then the random weights are chosen to be `alpha` with
#'   probability `Pw` and `1.0` otherwise. The latter choice is used
#'   in Theorem 2 in Meinshausen and Buhlmann (2010). Recommended
#'   values of `alpha` and `Pw` are \verb{[0.5, 0.2]} respectively.
#'
#' Stability selection can be run in parallel using multiple forked
#'   processes by setting `parallel = TRUE`. This requires
#'   [parallel::mclapply()] and is *not* available for Windows.
#'
#' @order 1
#' @family stability
#'
#' @param x A numeric \eqn{n \times p} matrix of predictive features
#'   containing `n` observation rows and `p` feature columns.
#'   Alternatively, a `stab_sel` class object if passing to one of
#'   its generic S3 methods.
#'
#' @param y The response variable. If kernel is "binomial", see
#'   [glmnet::glmnet()] for options. If kernel is "cox", a two column
#'   matrix with the event time in the first column and the censoring
#'   indicator (1 = event, 0 = censored) in the second column.
#'
#' @param kernel `character(1)`. A string describing the underlying
#'   model used for selection. Current options are:
#'   \itemize{
#'     \item "binomial" (logistic)
#'     \item "lasso"
#'     \item "ridge"
#'     \item "cox"
#'     \item "multinomial"
#'     \item "pca.sd"
#'     \item "pca.thresh"
#'   }
#'
#' @param n_iter `integer(1)`. Defining the number of
#'   sub-sampling iterations during each selection.
#'
#' @param parallel `logical(1)`. Should parallel processing
#'   via multiple cores be implemented? Must be on Linux or
#'   MacOS platform and have the \pkg{parallel}
#'   package installed. Otherwise defaults to 1 core.
#'
#' @param alpha `numeric(1)`. Value defining the weakness
#'   parameter for the randomized regularization.
#'   This is the minimum random weight applied to each
#'   beta coefficient in the regularization.
#'
#' @param Pw `numeric(1)`. Value defining probability of
#'   a weak weight, see `alpha`. If `Pw = NA` then the
#'   coefficient weights are sampled uniformly from `alpha` to `1.0`.
#'
#' @param standardize `logical(1)`. Whether the data should
#'   be centered and scaled. Passed to [glmnet()].
#'
#' @param lambda_min_ratio The minimum value of
#'   \eqn{\lambda / max(\lambda)} to use during
#'   the selection procedure. See [glmnet()].
#'
#' @param n_perm `integer(1)`. The number of permutations
#'   to use in calculating the empirical false positive rate.
#'
#' @param beta_threshold `numeric(1)`. Floating point value
#'   defining selection levels for `ridge regression`.
#'   Since ridge regression will not zero out coefficients,
#'   selection of coefficient curves by selection probability
#'   is not effective. Any variable having a coefficient with
#'   absolute value greater than or equal to `beta_threshold`
#'   will be selected.
#'
#' @param elastic_alpha `numeric(1)`. A value in `[0, 1]`. When 0,
#'   the results of [glmnet()] are equivalent to Ridge regression.
#'   When 1, the results are equivalent to Lasso. Any value between
#'   0 and 1 creates a compromise between the L1 and L2 penalty.
#'
#' @param lambda_pad `numeric(1)`. The lambda path is padded
#'   with high lambda values to produce more appealing stability
#'   paths for plotting. Occasionally, the degree of padding
#'   needs adjustment to produce better resolution at lower lambda.
#'   Typical values are: 20 (default), 15, 10, or 5.
#'
#' @param impute_outliers `logical(1)`. Should statistical
#'   outliers (\eqn{3 \times \sigma}) be imputed to approximate
#'   a Gaussian distribution during stability selection?
#'   See [wranglr::impute_outliers()].
#'
#' @param impute_n_sigma `numeric(1)`. Threshold standard deviation
#'   for outlier detection when imputing outliers. Ignored if
#'   `impute_outliers = FALSE`.
#'
#' @param r_seed `integer(1)`. Seed for the random number
#'   generator, for reproducibility.
#'
#' @param ... Additional arguments passed to one of the S3 methods
#'   for `stab_sel` class objects, generics include:
#'   * [plot.stab_sel()]
#'   * [print.stab_sel()]
#'   * [summary.stab_sel()]
#'
#' @return A `stab_sel` class object:
#'   \item{stabpath_matrix}{A matrix of \eqn{p \times lambda\_seq}
#'     containing stability selection probabilities. A row in this
#'     matrix corresponds to a stability selection path for a single feature.}
#'   \item{lambda}{the sequence of lambdas used for regularization. They
#'     correspond to the columns of `stabpath_matrix`.}
#'   \item{alpha}{the weakness parameter provided in the call.}
#'   \item{Pw}{the weak weight probability provided in the call.}
#'   \item{kernel}{the kernel used (e.g. `binomial`).}
#'   \item{n_iter}{The number of iterations used in computing
#'     the stability paths.}
#'   \item{standardize}{should the data be standardized prior to analysis?}
#'   \item{lambda_min_ratio}{?}
#'   \item{perm_data}{Logical. Is there permuted data to perform empirical FDR?}
#'   \item{permpath_list}{list containing information to calculated the
#'     permutation paths of the empirical false positive rate.}
#'   \item{perm_lambda}{The lambda used in the permuted lists.}
#'   \item{permpath_max}{max lambda for the permuted lists (I think).}
#'   \item{beta}{A matrix of the estimated betas calculated during
#'     the selection process.}
#'   \item{r_seed}{The random seed used.}
#'
#' @author Stu Field, Michael R. Mehan, and Robert Kirk DeLisle
#' @seealso [glmnet()], [get_stable_features()]
#'
#' @references Meinshausen, N and P Buhlmann. (2010). Stability selection.
#'   Journal of the Royal Statistical Society: Series B (Statistical
#'   Methodology), **72**: 417-473. doi: 10.1111/j.1467-9868.2010.00740.x
#'
#' @examples
#' # logistic regression
#' n_feat      <- 20L
#' n_samp      <- 2500L
#' x           <- matrix(rnorm(n_samp * n_feat), n_samp, n_feat)
#' fn <- function() paste0(sample(letters, 4L), collapse = "")
#' colnames(x) <- replicate(n_feat, paste0("f_", fn()))
#' y           <- sample(1:2, n_samp, replace = TRUE)
#' stab_sel    <- stability_selection(x, y)
#'
#' # Cox
#' xcox <- feature_matrix(stabilityselectr:::log_rfu(simdata))
#'
#' # Note: this works because colnames are already "time" and "status".
#' #   In 'real' datasets, you may need to rename the final matrix to
#' #   `time` and `status`.
#'
#' ycox <- data.matrix(select(simdata, time, status))
#' stab_sel_cox <- stability_selection(xcox, ycox, kernel = "cox", r_seed = 3)
#' @importFrom glmnet glmnet
#' @importFrom stats runif setNames
#' @importFrom parallel mclapply
#' @importFrom tibble tibble as_tibble
#' @export
stability_selection <- function(x, y = NULL,
                                kernel = c("binomial", "lasso", "ridge",
                                           "cox", "multinomial",
                                           "pca.sd", "pca.thresh"),
                                n_iter = 100, parallel = FALSE,
                                alpha = 0.5, Pw = 0.2,
                                n_perm = 0, standardize = TRUE,
                                lambda_min_ratio = 0.1,
                                beta_threshold = 0L,
                                elastic_alpha = 1.0,
                                lambda_pad = 20,
                                impute_outliers = FALSE,
                                impute_n_sigma = 3,
                                r_seed = 1234) {

  x <- data.matrix(x)  # convert to data matrix if df

  if ( any(is.na(x)) ) {
    ft_na <- names(which(apply(x, 2, function(.x) any(is.na(.x)))))
    stop(
      "NAs detected in `x`, this will cause `glmnet()` to fail\n",
      "Please check: ", value(ft_na), "\n",
      "Please remove the feature or use `wranglr::imputeNAs()`",
      call. = FALSE
    )
  }

  kernel <- match.arg(kernel)

  if ( impute_outliers ) {
    x <- apply(x, 2, impute_outliers, n_sigma = impute_n_sigma)
  }

  # Checks for parallel processing
  n_cores <- ifelse(parallel, get_cores(), 1L)

  signal_done(
    "Using kernel:", value(kernel), "and", value(n_cores),
    ifelse(n_cores > 1L, "cores (parallel)", "core (serial)")
  )

  # set seed for reproducibility
  # Claude says "L'Ecuyer-CMRG" is better cross
  # platform & in parallel processing
  rng <- RNGkind("L'Ecuyer-CMRG")
  withr::local_seed(r_seed)
  withr::defer(restore_rng_kind(rng))

  if ( n_perm > 0L ) {
    perm_seq <- seq_len(n_perm)
    if ( inherits(y, "matrix") ) {
      perm_y <- lapply(perm_seq, function(i) {
        y[sample(1:nrow(y), nrow(y)), ]}) # nolint: seq_linter.
    } else if ( is.vector(y) || is.factor(y) ) {
      perm_y <- lapply(perm_seq, function(i) {
        y[sample(1:length(y), length(y))]}) # nolint: seq_linter.
    } else {
      perm_y <- list()
    }
  }

  W <- .calcW(ncol(x), alpha, Pw, kernel)
  lambda_seq <- .get_lambda_seq(x = x, y = y, kernel = kernel, Pw = Pw,
                                standardize = standardize, alpha = alpha,
                                W = W, elastic_alpha = elastic_alpha,
                                lambda_min_ratio = lambda_min_ratio,
                                lambda_pad = lambda_pad)

  if ( length(lambda_seq) > 20L ) {
    perm_lambda <- lambda_seq[round(seq(1, length(lambda_seq),
                              length.out = round(length(lambda_seq) / 6)))]
  } else {
    perm_lambda <- lambda_seq
  }

  permpath_list <- NULL

  # THE ACTION ----
  stab_time <- system.time(
    stabpath_mat <- mclapply(seq(n_iter), function(i) {
      withr::local_seed(r_seed + i)
      calc_stability_paths(x, y,
                           kernel         = kernel,
                           lambda_seq     = lambda_seq,
                           alpha          = alpha,
                           Pw             = Pw,
                           beta_threshold = beta_threshold,
                           elastic_alpha  = elastic_alpha,
                           standardize    = standardize)
    }, mc.set.seed = FALSE, mc.cores = n_cores) |> Reduce(f = "+")
  )
  signal_done(
    "Stablity path run time:",
    paste0(value(round(stab_time["elapsed"], 4L)), "s")
  )

  if ( n_perm > 0L ) {
    perm_time <- system.time(
      permpath_list <- lapply(perm_seq, function(.j) {
        mclapply(seq(n_iter), function(i) {
          withr::local_seed(r_seed + i)
          calc_stability_paths(x, perm_y[[.j]],
                               kernel         = kernel,
                               lambda_seq     = perm_lambda,
                               alpha          = alpha,
                               Pw             = Pw,
                               standardize    = standardize,
                               beta_threshold = beta_threshold,
                               elastic_alpha  = elastic_alpha)
          }, mc.set.seed = FALSE, mc.cores = n_cores) |> Reduce(f = "+")
      }) |> setNames(sprintf("Perm_%03i", perm_seq))
    )
    signal_done(
      "Perm path run time:",
      paste0(value(round(perm_time["elapsed"], 4L)), "s")
    )
  }

  stabpath_mat <- stabpath_mat / n_iter / 2
  rownames(stabpath_mat) <- colnames(x)

  if ( n_perm > 0L ) {
    permpath_list <- lapply(permpath_list, function(.mat) {
      structure(.mat / n_iter / 2, dimnames = list(colnames(x), NULL))
    })
    perm_max_mat <- permpath_list |>
      lapply(function(.x) apply(.x, 1, max)) |>
      data.frame() |>        # preserve feature rownames
      rn2col("Feature") |>   # move features to column
      as_tibble()            # df => tibble
  } else {
    perm_max_mat <- tibble(Feature = character(0))
  }

  # get the betas for the appropriate lambda sequence
  glmnet_args <- list(
    x = x,
    y = y,
    standardize = standardize,
    lambda = lambda_seq,
    penalty.factor = W
  )

  glmnet_args <- switch(kernel,
    "binomial" = {
      glmnet_args$family  = kernel
      glmnet_args},
    "multinomial" = {
      glmnet_args$family  = kernel
      glmnet_args},
    "lasso"       = {
      glmnet_args$family = "gaussian"
      glmnet_args$alpha  = elastic_alpha
      glmnet_args},
    "ridge"       = {
      glmnet_args$family = "gaussian"
      glmnet_args$alpha  = 0
      glmnet_args},
    "cox"         = {
      glmnet_args$family   = kernel
      glmnet_args$nlambda  = 100   # nsteps
      glmnet_args$cox.ties = "efron"
      glmnet_args},
    NULL)  # pca kernels

  beta <- tryCatch(do.call(glmnet::glmnet, glmnet_args)$beta,
                   error = function(e) NULL) # NULL for pca kernels

  list(stabpath_matrix  = stabpath_mat,
       lambda           = lambda_seq,
       alpha            = alpha,
       Pw               = Pw,
       kernel           = kernel,
       n_iter           = n_iter,
       standardize      = standardize,
       impute_outliers  = impute_outliers,
       lambda_min_ratio = lambda_min_ratio,
       perm_data        = length(permpath_list) > 0L & n_perm > 0L,
       permpath_list    = permpath_list,
       perm_lambda      = perm_lambda,
       permpath_max     = perm_max_mat,
       beta             = beta,
       r_seed           = r_seed) |> add_class("stab_sel")
}
