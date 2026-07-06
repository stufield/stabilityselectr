
#' Internal to determine the lambda seq object
#'
#' @importFrom glmnet glmnet
#' @noRd
.get_lambda_seq <- function(x, y, kernel, Pw, standardize, alpha, W,
                            elastic_alpha, lambda_min_ratio, lambda_pad) {

  if ( kernel == "binomial" ) {

    lambda_seq <- glmnet::glmnet(
      x, y, family = kernel, standardize = standardize,
      lambda.min.ratio = lambda_min_ratio,
      penalty.factor = W)$lambda

  } else if ( kernel == "multinomial") {

    lambda_seq <- glmnet::glmnet(
      x, y, family = kernel, standardize = standardize,
      lambda.min.ratio = lambda_min_ratio,
      penalty.factor = W)$lambda

  } else if ( kernel == "pca.sd" ) {
    if ( floor(alpha) != alpha ) {
      stop(
        "Bad alpha: ", alpha,
        "\nPlease set to number of principal components to consider.",
        call. = FALSE
      )
    }
    max_lambda <- 10
    lambda_seq <- exp(
      seq(log(max_lambda), log(max_lambda * lambda_min_ratio),
          length.out = 100L) # hi resolution at low values
      ) # back to linear space

  } else if ( kernel == "pca.thresh" ) {

    if ( floor(alpha) != alpha ) {
      stop(
        "Bad alpha: ", alpha,
        "\nPlease set to number of principal components to consider.",
        call. = FALSE
      )
    }
    max_lambda <- 0.5
    lambda_seq <- exp(
      seq(log(max_lambda), log(max_lambda * lambda_min_ratio),
          length.out = 100)
      )

  } else if ( kernel == "lasso" ) {

    lambda_seq <- glmnet::glmnet(
      x, y, family = "gaussian", standardize = standardize,
      lambda.min.ratio = lambda_min_ratio, penalty.factor = W,
      alpha = elastic_alpha)$lambda

  } else if ( kernel == "ridge" ) {
    if ( elastic_alpha != 0 ) {
      stop(
        "Invalid `elastic_alpha =` argument for 'ridge' kernels (",
        value(elastic_alpha), "). Please set `elastic_alpha = 0`.",
        call. = FALSE
      )
    }
    lambda_seq <- glmnet::glmnet(
      x, y, family = "gaussian", standardize = standardize, alpha = 0,
      penalty.factor = W)$lambda

  } else if ( kernel == "cox" ) {

    stopifnot(
      "`y` response must be named `time` and `status` respectively." =
        identical(colnames(y), c("time", "status"))
    )

    lambda_seq <- glmnet::glmnet(
      x, y, nlambda = 100, family = kernel, standardize = standardize,
      cox.ties = "efron", # penalty.factor = W, (this differs from others)
      lambda.min.ratio = lambda_min_ratio)$lambda

  }

  if ( !kernel %in% c("pca.thresh", "pca.sd") ) {
    # need to pad lambda at the upper end
    seq_shift <- lambda_seq[1L] / lambda_seq[2L]
    pad <- rev(seq_len(lambda_pad)) |>  # default 1:20 & invert 20:1
      vapply(function(.x) lambda_seq[1L] * seq_shift^.x, 0.1)
    lambda_seq <- c(pad, lambda_seq)
  }

  lambda_seq
}
