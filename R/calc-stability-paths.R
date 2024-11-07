#' Calculate stability paths
#'
#' Internal function called by stability selection that returns the stability
#'   path matrix for *1* iteration (+ half/- half).
#'   See [stability_selection()] for more details.
#'
#' @param x A numeric matrix of predictive features.
#'   Contains a row for each observation.
#' @param y Response variable. If kernel is `l1-logistic`, it is a
#'   binary vector of class labels. If kernel is `Cox`, then it is a two column
#'   matrix of (**NEED HELP HERE -- not sure if correct; untested**).
#' @param kernel Character. String describing the underlying model
#'   used for selection. Options are:
#'   * l1-logistic
#'   * lasso
#'   * ridge
#'   * Cox
#'   * multinomial
#'   * pca.sd
#'   * pca.thresh
#' @param lambda_seq `numeric(n)`. A vector of lambdas used for regularization.
#' @param num_iter `integer(1)`. The number of subsampling iterations for
#'   the stability selection.
#' @param alpha `numeric(n)`. Value defining the weakness parameter for the
#'   randomized regularization. This is the minimum random weight applied to each
#'   beta coefficient in the regularization.
#' @param Pw `numeric(n)`. Value defining probability of a weak weight, see
#'   `alpha`. If `Pw = NA` then the coefficient weights are sampled
#'   uniformly from `alpha` to 1.
#' @param standardize `logical(1)`. Value passed to [glmnet()]. If
#'   `TRUE`, values will be standardized by [glmnet()].
#' @param beta_threshold `numeric(1)`. A value defining selection threshold for
#'   `Ridge regression` only. Since Ridge regression will not zero
#'   out coefficients, selection of coefficient curves by selection
#'   probability is not effective. Any variable having a coefficient with
#'   absolute value greater than or equal to `beta_threshold` will be selected.
#' @param elastic_alpha `numeric(1)`. A value between 0 and 1. When 0, the
#'   results of `glmnet` are equivalent to Ridge Regression. When 1, the results
#'   are equivalent to Lasso. Any value between 0 and 1 creates a compromise
#'   between L1 and L2 penalty.
#'
#' @return A matrix `features x lambda_seq` containing binary stability
#'   selections (`0 = no; 1 = yes`) for each feature at each value of lambda.
#'   Each row corresponds to a stability selection path for
#'   a single feature. **Note:** the values are for *both* halves of the
#'   stability selection, so can range from 0 (not selected) to
#'   2 (selected both times).
#' @author Michael R. Mehan, Stu Field, and Robert Kirk DeLisle
#'
#' @examples
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y           <- sample(1:2, n_samples, replace = TRUE)
#' })
#'
#' lambda_vec <- glmnet::glmnet(x, y, family = "binomial",
#'                              standardize = TRUE, lambda.min.ratio = 0.1,
#'                              penalty.factor = stats::runif(ncol(x)))$lambda
#'
#' path_matrix <- withr::with_seed(101,
#'   calc_stability_paths(x, y, "l1-logistic",
#'                        standardize = TRUE,
#'                        lambda_seq = lambda_vec,
#'                        alpha = 0.8, Pw = 0.5))
#'
#' # Cox kernel example
#' xcox <- strip_meta(log_rfu(sim_adat))
#'
#' # Note this works because colnames are already "time" and "status". In real
#' # datasets, need to rename the final matrix as "time" and "status".
#' ycox <- select(sim_adat, time, status) |> as.matrix()
#' cox_pm <- calc_stability_paths(xcox, ycox, "Cox", standardize = TRUE,
#'                                lambda_seq = c(0, 1, 100), alpha = 0.8, Pw = 0.5)
#' @importFrom stats prcomp runif
#' @importFrom glmnet glmnet
#' @noRd
calc_stability_paths <- function(x, y = NULL, kernel, lambda_seq, alpha, Pw,
                                 standardize, elastic_alpha, beta_threshold = 0L) {

  if ( beta_threshold == 0L && kernel == "ridge" ) {
    stop(
      "A `beta_threshold = 0` performs no feature selecting. Please ",
      "set a value of `beta_threshold > 0`.", call. = FALSE
    )
  }
  # nolint start (vars used in child scope)
  nobs         <- nrow(x)
  p            <- ncol(x)
  n_half1      <- floor(nobs / 2)
  lambda_seq_n <- length(lambda_seq)
  # nrow = no. features
  # ncol = lambda penalty sequence
  half1 <- sample(seq_len(nobs), n_half1, replace = FALSE)
  x1    <- x[half1, ]
  x2    <- x[-half1, ]

  if ( kernel == "Cox" ) {
    y1 <- y[half1, ]
    y2 <- y[-half1, ]
  } else if ( kernel %in% c("pca.sd", "pca.thresh") ) {
    # else is for unsupervised methods (PCA) which have no y
    NULL
  } else {
    y1 <- y[half1]
    y2 <- y[-half1]
  }
  # nolint end

  .stabPathFun <- switch(kernel,
    "l1-logistic" = .calc_l1,
    "lasso"       = .calc_lasso,
    "multinomial" = .calc_multinomial,
    "pca.sd"      = .calc_pca_sd,
    "pca.thresh"  = .calc_pca_thresh,
    "ridge"       = .calc_ridge,
    "Cox"         = .calc_Cox
  )

  .stabPathFun()
}


#' L1-Logistic Regression
#' @noRd
.calc_l1 <- function(e = parent.frame()) {
   # e: get all objects from parent environment

  W <- .calcW(e$p, e$alpha, e$Pw, e$kernel)

  # get selected features (0, 1)
  first_half <- glmnet::glmnet(e$x1, e$y1, family = "binomial",
                               lambda = e$lambda_seq,
                               standardize = e$standardize,
                               penalty.factor = W)$beta |>
    as.logical() |> as.numeric()
  second_half <- glmnet::glmnet(e$x2, e$y2, family = "binomial",
                                lambda = e$lambda_seq,
                                standardize = e$standardize,
                                penalty.factor = W)$beta |>
    as.logical() |> as.numeric()
  # convert to same dim matrix as stabpath_matrix
  first_half_mat  <- matrix(first_half, nrow = e$p)
  second_half_mat <- matrix(second_half, nrow = e$p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  # range: [0, 2]
  return(first_half_mat + second_half_mat)
}


#' Lasso Regression
#' @noRd
.calc_lasso <- function(e = parent.frame()) {
  # e: get all objects from parent environment

  W <- .calcW(e$p, e$alpha, e$Pw, e$kernel)

  first_half <- glmnet::glmnet(e$x1, e$y1, family = "gaussian",
                               lambda = e$lambda_seq,
                               standardize = e$standardize,
                               penalty.factor = W,
                               alpha = e$elastic_alpha)$beta |>
    as.logical() |> as.numeric()
  second_half <- glmnet::glmnet(e$x2, e$y2, family = "gaussian",
                                lambda = e$lambda_seq,
                                standardize = e$standardize,
                                penalty.factor = W,
                                alpha = e$elastic_alpha)$beta |>
    as.logical() |> as.numeric()
  # convert to same dim matrix as stabpath_matrix
  first_half_mat  <- matrix(first_half, nrow = e$p)   # should error out
  second_half_mat <- matrix(second_half, nrow = e$p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  return(first_half_mat + second_half_mat)
}


#' Multinomial Regression
#' @noRd
.calc_multinomial <- function(e = parent.frame()) {
  # e: get all objects from parent environment

  W <- .calcW(e$p, e$alpha, e$Pw, e$kernel)

  # get selected features (0, 1)
  # nolint next: undesirable_linter.
  first_half <- sapply(glmnet::glmnet(e$x1, e$y1, family = "multinomial",
                                      lambda = e$lambda_seq,
                                      standardize = e$standardize,
                                      penalty.factor = W)$beta, as.logical)

  # nolint next: undesirable_linter.
  second_half <- sapply(glmnet::glmnet(e$x2, e$y2, family = "multinomial",
                                       lambda = e$lambda_seq,
                                       standardize = e$standardize,
                                       penalty.factor = W)$beta, as.logical)

  # selected (0,1) in ANY of the classes & convert to matrix same dim as stabpath_matrix
  first_half_mat  <- matrix(as.numeric(apply(first_half, 1, any)), nrow = e$p)
  second_half_mat <- matrix(as.numeric(apply(second_half, 1, any)), nrow = e$p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  return(first_half_mat + second_half_mat)
}


#' PCA threshold
#' @noRd
.calc_pca_thresh <- function(e = parent.frame()) {
  # e: get all objects from parent environment
  pr1   <- stats::prcomp(e$x1, center = TRUE, scale. = e$standardize)
  stab1 <- lapply(e$lambda_seq, function(lambda) {
    apply(pr1$rotation[, 1:e$alpha], 2, function(.x) {
          abs(.x) > lambda }) |>    # lambda is loadings threshold
          rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  pr2   <- stats::prcomp(e$x2, center = TRUE, scale. = e$standardize)
  stab2 <- lapply(e$lambda_seq, function(lambda) {
    apply(pr2$rotation[, 1:e$alpha], 2, function(.x) {
          abs(.x) > lambda }) |>    # lambda is loadings threshold
          rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  return(stab1 + stab2)
}


#' PCA standard deviation
#' @noRd
.calc_pca_sd <- function(e = parent.frame()) {
  # e: get all objects from parent environment
  pr1   <- stats::prcomp(e$x1, center = TRUE, scale. = e$standardize)
  stab1 <- lapply(e$lambda_seq, function(lambda) {
    apply(pr1$rotation[, 1:e$alpha], 2, function(.x) {
          # MAD approx. for normal dist; Robust alternative to ML
          coefs <- helpr::fit_gauss(.x, mad = TRUE)
          as.numeric(abs(.x - coefs[1L]) > coefs[2L] * lambda) # lambda is SDs from mean.
          }) |> rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  pr2   <- stats::prcomp(e$x2, center = TRUE, scale. = e$standardize)
  stab2 <- lapply(e$lambda_seq, function(lambda) {
    apply(pr2$rotation[, 1:e$alpha], 2, function(.x) {
          coefs <- helpr::fit_gauss(.x, mad = TRUE)
          as.numeric(abs(.x - coefs[1L]) > coefs[2L] * lambda) # lambda is SDs from mean.
          }) |> rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  return(stab1 + stab2)
}


#' Ridge Regression
#' @noRd
.calc_ridge <- function(e = parent.frame()) {
  # e: get all objects from parent environment

  W <- .calcW(e$p, e$alpha, e$Pw, e$kernel)

  # all variables will (almost) always have non-zero coefficients
  # in ridge regression; implement a threshold for selection
  first_half <- apply(glmnet::glmnet(e$x1, e$y1, family = "gaussian",
                                     alpha = 0, lambda = e$lambda_seq,
                                     standardize = e$standardize,
                                     penalty.factor = W)$beta,
                      1:2,
                      function(.x) ifelse(abs(.x) < e$beta_threshold, 0L, .x)) |>
                  as.logical() |> as.numeric()

  second_half <- apply(glmnet::glmnet(e$x2, e$y2, family = "gaussian",
                                      alpha = 0, lambda = e$lambda_seq,
                                      standardize = e$standardize,
                                      penalty.factor = W)$beta,
                       1:2,
                       function(.x) ifelse(abs(.x) < e$beta_threshold, 0L, .x)) |>
                   as.logical() |> as.numeric()

  # convert to same dim matrix as stabpath_matrix
  first_half_mat  <- matrix(first_half, nrow = e$p)
  second_half_mat <- matrix(second_half, nrow = e$p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  return(first_half_mat + second_half_mat)
}


#' Cox model
#' @noRd
.calc_Cox <- function(e = parent.frame()) {
  # e: get all objects from parent environment

  err_cnt <- 0
  W       <- .calcW(e$p, e$alpha, e$Pw, e$kernel)

  # Try catch here to handle glmnet failures ...
  # first half
  first_half <- tryCatch({
    glmnet::glmnet(e$x1, e$y1, family = "cox", lambda = e$lambda_seq,
                   standardize = FALSE, penalty.factor = W)$beta
    }, error = function(err) {
      # error handler picks up where error was generated
      print(sprintf("`calc_stability_paths()` glmnet error: %s", err))
      rep(0, length = e$p)
    }
  )

  if ( sum(first_half == 0) == e$p ) {
    err_cnt <- err_cnt + 1
  }

  W <- .calcW(e$p, e$alpha, e$Pw, e$kernel)

  # second half
  second_half <- tryCatch({
    glmnet::glmnet(e$x2, e$y2, family = "cox", lambda = e$lambda_seq,
                   standardize = FALSE, penalty.factor = W)$beta
    }, error = function(err) {
      # error handler picks up where error was generated
      print(sprintf("`calc_stability_paths()` glmnet error: %s", err))
      rep(0, length = e$p)
    }
  )

  if ( sum(second_half == 0) == e$p ) {
    err_cnt <- err_cnt + 1
  }

  if ( err_cnt > 0 ) {
    warning("Detected ", value(err_cnt), " `glmnet()` *Cox* errors",
            call. = FALSE)
  }
  first_half_mat  <- matrix(as.numeric(as.logical(first_half)), nrow = e$p)
  second_half_mat <- matrix(as.numeric(as.logical(second_half)), nrow = e$p)
  return(first_half_mat + second_half_mat)
}


#' @importFrom stats runif
#' @noRd
.calcW <- function(p, alpha, Pw, kernel) {

  if ( is.na(Pw) && kernel != "Cox" ) {
    W <- stats::runif(p, alpha, 1)
  } else {
    W    <- rep(1, length = p)
    draw <- stats::runif(p, 0, 1)
    # as per the paper; RKD 2013 Nov 8
    W[draw < Pw] <- ifelse(kernel %in% c("ridge", "lasso"), 1L / alpha, alpha)
  }
  W
}
