#' Calculate stability paths
#'
#' Internal function called by stability selection that
#'   returns the stability path matrix for *1* iteration (+ half/- half).
#'   See [stability_selection()] for more details.
#'
#' @inheritParams stability_selection
#'
#' @return A matrix `features x lambda_seq` containing binary
#'   stability selections (`0 = no; 1 = yes`) for each feature
#'   at each value of lambda. Each row corresponds to a
#'   stability selection path for a single feature.
#'   **Note:** the values are for *both* halves of the
#'   stability selection, so values can be 0 (not selected), 1,
#'   or 2 (selected in both halves).
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
#'   calc_stability_paths(x, y, kernel = "l1-logistic", standardize = TRUE,
#'                        alpha = 0.8, Pw = 0.8, lambda_seq = lambda_vec))
#'
#' # Cox kernel example
#' xcox <- feature_matrix(stabilityselectr:::log_rfu(simdata))
#'
#' # Note this works because colnames are already "time" and "status". In real
#' # datasets, need to rename the final matrix as "time" and "status".
#' ycox <- select(simdata, time, status) |> as.matrix()
#' cox_pm <- calc_stability_paths(xcox, ycox, "cox", standardize = TRUE,
#'                                lambda_seq = c(0, 1, 100), alpha = 0.8, Pw = 0.5)
#' @importFrom stats prcomp runif
#' @importFrom glmnet glmnet
#' @noRd
calc_stability_paths <- function(x, y = NULL, kernel, lambda_seq, alpha, Pw,
                                 standardize, elastic_alpha, beta_threshold = 0L) {

  if ( kernel == "ridge" && beta_threshold == 0L ) {
    stop(
      "A `beta_threshold = 0` performs no feature selecting. Please ",
      "set a value of `beta_threshold > 0`.", call. = FALSE
    )
  }

  # nolint start
  nobs         <- nrow(x)
  p            <- ncol(x)
  n_half1      <- floor(nobs / 2)
  lambda_seq_n <- length(lambda_seq)
  # nrow = no. features
  # ncol = lambda penalty sequence
  half1 <- sample(seq_len(nobs), n_half1, replace = FALSE)
  x1    <- x[half1, ]
  x2    <- x[-half1, ]

  if ( kernel %in% c("pca.sd", "pca.thresh") ) {
    # for unsupervised methods (PCA) which have no y
    NULL
  } else {
    if (is.null(dim(y))) {  # response a vector
      y1 <- y[half1]
      y2 <- y[-half1]
    } else {                # response not a vector
      y1 <- y[half1, ]
      y2 <- y[-half1, ]
    }
  }
  # nolint end

  .stab_path_type <- switch(kernel,
    "l1-logistic" = .calc_l1,
    "lasso"       = .calc_lasso,
    "multinomial" = .calc_multinomial,
    "pca.sd"      = .calc_pca_sd,
    "pca.thresh"  = .calc_pca_thresh,
    "ridge"       = .calc_ridge,
    "cox"         = .calc_cox,
    stop("Could not find appropriate ", value("kernel"),
         ". Please check call.", call. = FALSE)
  )

  env <- function() parent.frame(n = 1)
  get("attach")(env(), name = "ss_tmp", pos = 2L,
                warn.conflicts = FALSE)
  on.exit(get("detach")("ss_tmp", character.only = TRUE))
  .stab_path_type()
}


#' L1-Logistic Regression
#' @noRd
.calc_l1 <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  # get selected features (0, 1)
  first_half <- glmnet::glmnet(x1, y1, family = "binomial",
                               lambda = lambda_seq,
                               standardize = standardize,
                               penalty.factor = W)$beta |>
    as.logical() |> as.numeric()
  second_half <- glmnet::glmnet(x2, y2, family = "binomial",
                                lambda = lambda_seq,
                                standardize = standardize,
                                penalty.factor = W)$beta |>
    as.logical() |> as.numeric()
  # convert to same dim matrix as stabpath_matrix
  first_half_mat  <- matrix(first_half, nrow = p)
  second_half_mat <- matrix(second_half, nrow = p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  # range: [0, 2]
  first_half_mat + second_half_mat
}


#' Lasso Regression
#' @noRd
.calc_lasso <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  first_half <- glmnet::glmnet(x1, y1, family = "gaussian",
                               lambda = lambda_seq,
                               standardize = standardize,
                               penalty.factor = W,
                               alpha = elastic_alpha)$beta |>
    as.logical() |> as.numeric()
  second_half <- glmnet::glmnet(x2, y2, family = "gaussian",
                                lambda = lambda_seq,
                                standardize = standardize,
                                penalty.factor = W,
                                alpha = elastic_alpha)$beta |>
    as.logical() |> as.numeric()
  # convert to same dim matrix as stabpath_matrix
  first_half_mat  <- matrix(first_half, nrow = p)   # should error out
  second_half_mat <- matrix(second_half, nrow = p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  first_half_mat + second_half_mat
}


#' Multinomial Regression
#' @noRd
.calc_multinomial <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  # get selected features (0, 1)
  # nolint next: undesirable_linter.
  first_half <- sapply(glmnet::glmnet(x1, y1, family = "multinomial",
                                      lambda = lambda_seq,
                                      standardize = standardize,
                                      penalty.factor = W)$beta, as.logical)

  # nolint next: undesirable_linter.
  second_half <- sapply(glmnet::glmnet(x2, y2, family = "multinomial",
                                       lambda = lambda_seq,
                                       standardize = standardize,
                                       penalty.factor = W)$beta, as.logical)

  # selected (0,1) in ANY of the classes & convert to matrix same dim as stabpath_matrix
  first_half_mat  <- matrix(as.numeric(apply(first_half, 1, any)), nrow = p)
  second_half_mat <- matrix(as.numeric(apply(second_half, 1, any)), nrow = p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  first_half_mat + second_half_mat
}


#' PCA threshold
#' @noRd
.calc_pca_thresh <- function() {
  pr1   <- stats::prcomp(x1, center = TRUE, scale. = standardize)
  stab1 <- lapply(lambda_seq, function(lambda) {
    apply(pr1$rotation[, 1:alpha], 2, function(.x) {
          abs(.x) > lambda }) |>    # lambda is loadings threshold
          rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  pr2   <- stats::prcomp(x2, center = TRUE, scale. = standardize)
  stab2 <- lapply(lambda_seq, function(lambda) {
    apply(pr2$rotation[, 1:alpha], 2, function(.x) {
          abs(.x) > lambda }) |>    # lambda is loadings threshold
          rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  stab1 + stab2
}


#' PCA standard deviation
#' @noRd
.calc_pca_sd <- function() {
  pr1   <- stats::prcomp(x1, center = TRUE, scale. = standardize)
  stab1 <- lapply(lambda_seq, function(lambda) {
    apply(pr1$rotation[, 1:alpha], 2, function(.x) {
          # MAD approx. for normal dist; Robust alternative to ML
          coefs <- helpr::fit_gauss(.x, mad = TRUE)
          as.numeric(abs(.x - coefs[1L]) > coefs[2L] * lambda) # lambda is SDs from mean.
          }) |> rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  pr2   <- stats::prcomp(x2, center = TRUE, scale. = standardize)
  stab2 <- lapply(lambda_seq, function(lambda) {
    apply(pr2$rotation[, 1:alpha], 2, function(.x) {
          coefs <- helpr::fit_gauss(.x, mad = TRUE)
          as.numeric(abs(.x - coefs[1L]) > coefs[2L] * lambda) # lambda is SDs from mean.
          }) |> rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  stab1 + stab2
}


#' Ridge Regression
#' @noRd
.calc_ridge <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  # all variables will (almost) always have non-zero coefficients
  # in ridge regression; implement a threshold for selection
  first_half <- apply(glmnet::glmnet(x1, y1, family = "gaussian",
                                     alpha = 0, lambda = lambda_seq,
                                     standardize = standardize,
                                     penalty.factor = W)$beta,
                      1:2,
                      function(.x) ifelse(abs(.x) < beta_threshold, 0L, .x)) |>
                  as.logical() |> as.numeric()

  second_half <- apply(glmnet::glmnet(x2, y2, family = "gaussian",
                                      alpha = 0, lambda = lambda_seq,
                                      standardize = standardize,
                                      penalty.factor = W)$beta,
                       1:2,
                       function(.x) ifelse(abs(.x) < beta_threshold, 0L, .x)) |>
                   as.logical() |> as.numeric()

  # convert to same dim matrix as stabpath_matrix
  first_half_mat  <- matrix(first_half, nrow = p)
  second_half_mat <- matrix(second_half, nrow = p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  first_half_mat + second_half_mat
}


#' Cox model
#' @noRd
.calc_cox <- function() {
  err_cnt <- 0L

  # tryCatch to handle glmnet failures ...
  # first half
  first_half <- tryCatch({
    glmnet::glmnet(
      x1, y1, family = "cox", standardize = FALSE,
      lambda = lambda_seq, cox.ties = "breslow",
      penalty.factor = .calcW(p, alpha, Pw, kernel)
    )$beta
    }, error = function(err) {
      # error handler picks up where error was generated
      print(sprintf("`calc_stability_paths()` glmnet error: %s", err))
      rep(0, length = p)
    }
  )

  if ( sum(first_half == 0) == p ) {
    err_cnt <- err_cnt + 1L
  }


  # second half
  second_half <- tryCatch({
    glmnet::glmnet(
      x2, y2, family = "cox", standardize = FALSE,
      lambda = lambda_seq, cox.ties = "breslow",
      penalty.factor = .calcW(p, alpha, Pw, kernel)
    )$beta
    }, error = function(err) {
      print(sprintf("`calc_stability_paths()` glmnet error: %s", err))
      rep(0, length = p)
    }
  )

  if ( sum(second_half == 0) == p ) {
    err_cnt <- err_cnt + 1L
  }

  if ( err_cnt > 0L ) {
    warning("Detected ", value(err_cnt), " `glmnet()` *cox* errors",
            call. = FALSE)
  }
  first_half_mat  <- matrix(as.numeric(as.logical(first_half)), nrow = p)
  second_half_mat <- matrix(as.numeric(as.logical(second_half)), nrow = p)
  first_half_mat + second_half_mat
}


#' @importFrom stats runif
#' @noRd
.calcW <- function(p, alpha, Pw, kernel) {

  if ( is.na(Pw) && kernel != "cox" ) {
    W <- stats::runif(p, alpha, 1)
  } else {
    W    <- rep(1, length = p)
    draw <- stats::runif(p, 0, 1)
    # as per the paper; RKD 2013 Nov 8
    W[draw < Pw] <- ifelse(kernel %in% c("ridge", "lasso"), 1L / alpha, alpha)
  }
  W
}
