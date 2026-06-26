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
calc_stability_paths <- function(x, y = NULL, kernel,
                                 lambda_seq, alpha, Pw,
                                 standardize = NULL,
                                 elastic_alpha = NULL,
                                 beta_threshold = 0L) {

  if ( kernel == "ridge" && beta_threshold == 0L ) {
    stop(
      "In ridge regression, a `beta_threshold = 0` performs no ",
      "feature selection.\nPlease set a value of `beta_threshold > 0`.",
      call. = FALSE
    )
  }

  # nolint start
  nobs         <- nrow(x)
  p            <- ncol(x)
  n_half1      <- floor(nobs / 2)
  lambda_seq_n <- length(lambda_seq)
  half1 <- sample(seq_len(nobs), n_half1, replace = FALSE)
  x1    <- x[half1, ]
  x2    <- x[-half1, ]

  if ( kernel %in% c("pca.sd", "pca.thresh") ) {
    NULL   # unsupervised methods (PCA) have no y
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
  h1 <- glmnet::glmnet(x1, y1, family = "binomial",
                       lambda = lambda_seq,
                       standardize = standardize,
                       penalty.factor = W)$beta |>
    as.logical() |> as.numeric()
  h2 <- glmnet::glmnet(x2, y2, family = "binomial",
                       lambda = lambda_seq,
                       standardize = standardize,
                       penalty.factor = W)$beta |>
    as.logical() |> as.numeric()
  # convert to same dim matrix as stabpath_matrix
  first_half  <- matrix(h1, nrow = p)
  second_half <- matrix(h2, nrow = p)

  # add no. of selections for each feature by each lambda to stabpath_matrix
  # range: [0, 2]
  first_half + second_half
}


#' Lasso Regression
#' @noRd
.calc_lasso <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  h1 <- glmnet::glmnet(x1, y1, family = "gaussian",
                       lambda = lambda_seq,
                       standardize = standardize,
                       penalty.factor = W,
                       alpha = elastic_alpha)$beta |>
    as.logical() |> as.numeric()
  h2 <- glmnet::glmnet(x2, y2, family = "gaussian",
                       lambda = lambda_seq,
                       standardize = standardize,
                       penalty.factor = W,
                       alpha = elastic_alpha)$beta |>
    as.logical() |> as.numeric()
  first_half  <- matrix(h1, nrow = p)
  second_half <- matrix(h2, nrow = p)
  first_half + second_half
}


#' Ridge Regression
#' @noRd
.calc_ridge <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  # all variables will (almost) always have non-zero coefficients
  # in ridge regression; implement a threshold for selection
  h1 <- apply(glmnet::glmnet(x1, y1, family = "gaussian",
                             alpha = 0, lambda = lambda_seq,
                             standardize = standardize,
                             penalty.factor = W)$beta,
                      1:2,
                      function(.x) ifelse(abs(.x) < beta_threshold, 0L, .x)) |>
                  as.logical() |> as.numeric()

  h2 <- apply(glmnet::glmnet(x2, y2, family = "gaussian",
                             alpha = 0, lambda = lambda_seq,
                             standardize = standardize,
                             penalty.factor = W)$beta,
                       1:2,
                       function(.x) ifelse(abs(.x) < beta_threshold, 0L, .x)) |>
                   as.logical() |> as.numeric()

  first_half  <- matrix(h1, nrow = p)
  second_half <- matrix(h2, nrow = p)
  first_half + second_half
}


#' Cox Regression
#' @noRd
.calc_cox <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  h1 <- tryCatch({   # tryCatch to handle glmnet failures
    glmnet::glmnet(
      x1, y1, family = "cox", standardize = FALSE,
      lambda = lambda_seq, cox.ties = "efron",
      penalty.factor = W)$beta
    }, error = function(e) {
      cat("\n`glmnet()` error in `calc_stability_paths()`:\n",
          e$message, "\n")
      rep.int(0, p)
    }
  )

  err_cnt <- 0L

  if ( sum(h1 == 0) == p ) {
    err_cnt <- err_cnt + 1L
  }

  h2 <- tryCatch({
    glmnet::glmnet(
      x2, y2, family = "cox", standardize = FALSE,
      lambda = lambda_seq, cox.ties = "efron",
      penalty.factor = W)$beta
    }, error = function(e) {
      cat("\n`glmnet()` error in `calc_stability_paths()`:\n",
          e$message, "\n")
      rep.int(0, p)
    }
  )

  if ( sum(h2 == 0) == p ) {
    err_cnt <- err_cnt + 1L
  }

  if ( err_cnt > 0L ) {
    warning("Detected ", value(err_cnt), " `glmnet()` *cox* errors",
            call. = FALSE)
  }

  first_half  <- matrix(as.numeric(as.logical(h1)), nrow = p)
  second_half <- matrix(as.numeric(as.logical(h2)), nrow = p)
  first_half + second_half
}


#' Multinomial Regression
#' @noRd
.calc_multinomial <- function() {

  W <- .calcW(p, alpha, Pw, kernel)

  # get selected features (0, 1)
  # nolint next: undesirable_linter.
  h1 <- sapply(glmnet::glmnet(x1, y1, family = "multinomial",
                              lambda = lambda_seq,
                              standardize = standardize,
                              penalty.factor = W)$beta, as.logical)

  # nolint next: undesirable_linter.
  h2 <- sapply(glmnet::glmnet(x2, y2, family = "multinomial",
                              lambda = lambda_seq,
                              standardize = standardize,
                              penalty.factor = W)$beta, as.logical)

  # selected (0,1) in ANY of the classes & convert to matrix same dim as stabpath_matrix
  first_half  <- matrix(as.numeric(apply(h1, 1, any)), nrow = p)
  second_half <- matrix(as.numeric(apply(h2, 1, any)), nrow = p)
  first_half + second_half
}


#' PCA threshold
#' @noRd
.calc_pca_thresh <- function() {
  pr1 <- stats::prcomp(x1, center = TRUE, scale. = standardize)
  pr2 <- stats::prcomp(x2, center = TRUE, scale. = standardize)
  first_half <- lapply(lambda_seq, function(lambda) {
    apply(pr1$rotation[, 1:alpha], 2, function(.x) {
          abs(.x) > lambda }) |>    # lambda is loadings threshold
          rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  second_half <- lapply(lambda_seq, function(lambda) {
    apply(pr2$rotation[, 1:alpha], 2, function(.x) {
          abs(.x) > lambda }) |>    # lambda is loadings threshold
          rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  first_half + second_half
}


#' PCA standard deviation
#' @noRd
.calc_pca_sd <- function() {
  pr1 <- stats::prcomp(x1, center = TRUE, scale. = standardize)
  pr2 <- stats::prcomp(x2, center = TRUE, scale. = standardize)
  first_half <- lapply(lambda_seq, function(lambda) {
    apply(pr1$rotation[, 1:alpha], 2, function(.x) {
          # MAD approx. for normal dist; Robust alternative to ML
          coefs <- helpr::fit_gauss(.x, mad = TRUE)
          as.numeric(abs(.x - coefs[1L]) > coefs[2L] * lambda) # lambda is SDs from mean.
         }) |> rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  second_half <- lapply(lambda_seq, function(lambda) {
    apply(pr2$rotation[, 1:alpha], 2, function(.x) {
          coefs <- helpr::fit_gauss(.x, mad = TRUE)
          as.numeric(abs(.x - coefs[1L]) > coefs[2L] * lambda) # lambda is SDs from mean.
         }) |> rowSums() |> as.logical() |> as.numeric()
  }) |> data.frame() |> as.matrix() |> unname()
  first_half + second_half
}
