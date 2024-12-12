
# Testing ----
# L1-logistic ----
test_that("`calc_stability_paths()` generates expected output for kernel = l1-logistic", {
  withr::with_seed(101, {
    lambda_vec <- glmnet::glmnet(x, y, family = "binomial",
                                 standardize = TRUE, lambda.min.ratio = 0.1,
                                 penalty.factor = stats::runif(ncol(x)))$lambda
    path_matrix <- calc_stability_paths(x, y, kernel = "l1-logistic",
                                        lambda_seq = lambda_vec,
                                        alpha = 0.8, Pw = 0.5,
                                        elastic_alpha = 1,     # not used
                                        standardize = TRUE)
  })
  expect_true(is.matrix(path_matrix))
  expect_equal(dim(path_matrix), c(n_feat, length(lambda_vec)))
  expect_equal(sum(path_matrix), 81)
  expect_true(all(path_matrix >= 0))          # range in [0,2]
  expect_true(all(path_matrix <= 2))          # range in [0,2]
  expect_equal(c(table(path_matrix)), c("0" = 1921, "1" = 77, "2" = 2))
  expect_equal(diag(path_matrix), rep(0, n_feat))
  expect_equal(colSums(path_matrix),
               rep(c(0, 1, 2, 3, 5, 6), c(74, 4, 11, 1, 8, 2)))
  # which features are never selected
  expect_equal(which(rowSums(path_matrix) == 0),
               c(1:8, 11:12, 14:16, 19, 20))  # zero seln
  expect_equal(rowSums(path_matrix),
               c(rep(0, 8), 26, 10, 0, 0, 22, 0, 0, 0, 13, 10, 0, 0))
})

# Lasso ----
test_that("`calc_stability_paths()` generates expected output for kernel = lasso", {
  L <- 100 + 20
  withr::with_seed(1, {
    y <- rnorm(n_samples)
    path_matrix <- calc_stability_paths(x, y, kernel = "lasso",
                                        lambda_seq = seq(0.25, 0.016, length.out = L),
                                        alpha = 0.8, Pw = 0.5,
                                        elastic_alpha = 1,    # 1 is lasso
                                        beta_threshold = 0,
                                        standardize = TRUE)
  })
  expect_true(is.matrix(path_matrix))
  expect_equal(dim(path_matrix), c(n_feat, L))
  expect_equal(sum(path_matrix), 1607)
  expect_true(all(path_matrix >= 0))          # range in [0,2]
  expect_true(all(path_matrix <= 2))          # range in [0,2]
  expect_equal(c(table(path_matrix)), c("0" = 1199, "1" = 795, "2" = 406))
  expect_equal(diag(path_matrix), c(0, 1, 0, 0, 1, 0, 0, 0, 1, 1, rep(0, 10)))
  expect_equal(colSums(path_matrix),
    c(3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6,
      6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8,
      8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 13,
      13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 15, 16, 16, 16, 17, 19, 19,
      19, 19, 19, 19, 19, 19, 19, 19, 20, 21, 21, 21, 21, 21, 21, 21, 21, 22,
      22, 22, 23, 24, 25, 25, 25, 26, 26, 26, 29, 29, 31, 31, 31, 31, 33, 34, 34)
  )
  # which features are never selected
  expect_equal(which(rowSums(path_matrix) == 0), integer(0))  # zero seln
  expect_equal(rowSums(path_matrix),
               c(23, 192, 28, 68, 199, 165, 9, 18, 134, 120, 90, 55, 112,
                 60, 55, 100, 70, 99, 3, 7))
})

# Cox ----
test_that("`calc_stability_paths()` generates expected output with kernel = 'Cox'", {
  xcox <- feature_matrix(log_rfu(simdata))
  ycox <- select(simdata, time, status) |> as.matrix()
  cox_pm <- withr::with_seed(
    1234,
    calc_stability_paths(xcox, ycox, kernel = "Cox", standardize = TRUE,
                       lambda_seq = c(0.01, 0.1), elastic_alpha = 1,     # not used
                       alpha = 0.8, Pw = 0.5)
  )
  expect_true(is.matrix(cox_pm))
  expect_equal(dim(cox_pm), c(40, 2))
  expect_equal(sum(cox_pm), 19)
  expect_true(all(cox_pm >= 0))          # range in [0,2]
  expect_true(all(cox_pm <= 2))          # range in [0,2]
  expect_equal(c(table(cox_pm)), c("0" = 63, "1" = 15, "2" = 2))
  expect_equal(diag(cox_pm), c(0, 0))
  expect_equal(colSums(cox_pm), c(0, 19))

  # which features are never selected
  expect_equal(which(rowSums(cox_pm) == 0),
               c(1L, 2L, 3L, 4L, 5L, 8L, 10L, 16L, 19L, 21L, 22L, 23L, 24L,
                 26L, 27L, 29L, 31L, 33L, 34L, 35L, 38L, 39L, 40L))
  # zero seln
  expect_equal(rowSums(cox_pm),
               c(0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 2, 1, 1, 1, 0, 1, 1, 0, 1, 0,
                 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0))
})

# Ridge ----
test_that("`calc_stability_paths()` generates expected output for kernel = ridge", {
  L <- 100 + 20
  withr::with_seed(10, {
    y <- rnorm(n_samples)
    path_matrix <- calc_stability_paths(x, y, kernel = "ridge",
                                        lambda_seq = seq(0.25, 0.016, length.out = L),
                                        alpha = 0.8, Pw = 0.5,
                                        elastic_alpha = 0,   # 0 is ridge
                                        beta_threshold = 0.1,
                                        standardize = TRUE)
  })
  expect_true(is.matrix(path_matrix))
  expect_equal(dim(path_matrix), c(n_feat, L))
  expect_equal(sum(path_matrix), 2346)
  expect_true(all(path_matrix >= 0))          # range in [0,2]
  expect_true(all(path_matrix <= 2))          # range in [0,2]
  expect_equal(c(table(path_matrix)), c("0" = 621, "1" = 1212, "2" = 567))
  expect_equal(diag(path_matrix), c(0, 1, 2, 0, 0, 0, 2, 2, 1, 1, 1, 1,
                                    2, 1, 0, 1, 1, 0, 1, 1))
  expect_equal(colSums(path_matrix),
               c(18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
                 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
                 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
                 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
                 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 20,
                 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21,
                 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
                 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24))
  # which features are never selected
  expect_equal(which(rowSums(path_matrix) == 0), c(5, 6, 18))  # zero seln
  expect_equal(rowSums(path_matrix),
               c(30, 124, 240, 31, 0, 0, 240, 240, 120, 120, 120, 120, 240,
                 146, 38, 166, 120, 0, 120, 131))

  expect_error(
    calc_stability_paths(x, y, kernel = "ridge",
                       lambda_seq = seq(0.25, 0.016, length.out = L),
                       alpha = 0.8, Pw = 0.5,
                       elastic_alpha = 0,   # 0 is ridge
                       beta_threshold = 0,  # This must be > 0
                       standardize = TRUE),
    paste("A `beta_threshold = 0` performs no feature selecting.",
          "Please set a value of `beta_threshold > 0`.")
    )
})

# Multinomial ----
test_that("`calc_stability_paths()` generates expected output for kernel = multinomial", {
  L <- 100 + 20
  withr::with_seed(10, {
    y <- sample(1:2, n_samples, replace = TRUE)
    path_matrix <- calc_stability_paths(x, y, "multinomial",
                                      lambda_seq = seq(0.25, 0.016, length.out = L),
                                      alpha = 0.8, Pw = 0.5,
                                      elastic_alpha = 1,     # not used
                                      standardize = TRUE)
  })
  expect_true(is.matrix(path_matrix))
  expect_equal(dim(path_matrix), c(n_feat, L))
  expect_equal(sum(path_matrix), 754)
  expect_true(all(path_matrix >= 0))          # range in [0,2]
  expect_true(all(path_matrix <= 2))          # range in [0,2]
  expect_equal(c(table(path_matrix)), c("0" = 1744, "1" = 558, "2" = 98))
  expect_equal(diag(path_matrix), rep(0, 20))
  expect_equal(colSums(path_matrix),
               c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2,
                 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5,
                 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 7, 10, 10, 10, 11,
                 12, 12, 12, 12, 12, 12, 13, 13, 15, 15, 16, 16, 17, 18, 18,
                 19, 22, 22, 22, 23, 23, 23, 23, 25, 26, 28))
  # which features are never selected
  expect_equal(which(rowSums(path_matrix) == 0), 19)  # zero seln
  expect_equal(rowSums(path_matrix),
               c(11, 30, 75, 33, 18, 105, 1, 57, 83, 28, 31, 15,
                 58, 10, 94, 16, 13, 73, 0, 3))
})

# PCA thresh ----
test_that("`calc_stability_paths` generates expected output for kernel = pca.thresh", {
  max_lambda <- 0.5
  lambda_vec <- seq(log(max_lambda),
                    log(max_lambda * 0.1),
                    length.out = 100) |> exp()

  path_matrix <- withr::with_seed(
    1,
    calc_stability_paths(x, kernel = "pca.thresh",
                       lambda_seq = lambda_vec,
                       alpha = 2, Pw = 0.5,
                       beta_threshold = 0,
                       elastic_alpha = 1,     # not used
                       standardize = TRUE)
  )
  expect_true(is.matrix(path_matrix))
  expect_equal(dim(path_matrix), c(n_feat, length(lambda_vec)))
  expect_equal(sum(path_matrix), 2736)
  expect_true(all(path_matrix >= 0))          # range in [0,2]
  expect_true(all(path_matrix <= 2))          # range in [0,2]
  expect_equal(c(table(path_matrix)), c("0" = 474, "1" = 316, "2" = 1210))
  expect_equal(diag(path_matrix), c(rep(0, 8), 2, rep(0, 9), 2, 0))
  expect_equal(colSums(path_matrix),
               c(0, 0, 0, 1, 1, 1, 1, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7,
                 7, 8, 10, 14, 17, 19, 19, 21, 21, 21, 22, 25, 25, 26, 27,
                 27, 27, 27, 28, 28, 28, 29, 30, 31, 33, 33, 33, 33, 33, 33,
                 34, 34, 34, 34, 34, 35, 36, 36, 38, 38, 38, 38, 38, 38, 38,
                 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39,
                 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
                 39, 39, 39, 39, 39, 39))
  # which features are never selected
  expect_equal(which(rowSums(path_matrix) == 0), integer(0))  # zero seln
  expect_equal(rowSums(path_matrix),
               c(135, 178, 142, 153, 133, 116, 115, 155, 185, 77, 156, 137,
                 157, 127, 89, 105, 114, 147, 184, 131))
})

# PCA SD ----
test_that("`calc_stability_paths` generates expected output for kernel = pca.sd", {
  max_lambda <- 10
  lambda_vec <- seq(log(max_lambda),
                    log(max_lambda * 0.1),
                    length.out = 100) |> exp()

  path_matrix <- withr::with_seed(
    1,
    calc_stability_paths(x, y, "pca.sd",
                       lambda_seq = lambda_vec,
                       alpha = 3, Pw = 0.5,
                       beta_threshold = 0,
                       elastic_alpha = 1,     # not used
                       standardize = TRUE)
  )
  expect_true(is.matrix(path_matrix))
  expect_equal(dim(path_matrix), c(n_feat, length(lambda_vec)))
  expect_equal(sum(path_matrix), 614)
  expect_true(all(path_matrix >= 0))          # range in [0,2]
  expect_true(all(path_matrix <= 2))          # range in [0,2]
  expect_equal(c(table(path_matrix)), c("0" = 1542, "1" = 302, "2" = 156))
  expect_equal(diag(path_matrix), rep(0, n_feat))
  expect_equal(colSums(path_matrix),
               c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 1, 2, 2, 2, 3, 4, 6, 7, 8, 9, 9, 11, 11,
                 11, 12, 12, 13, 13, 15, 16, 18, 19, 20, 22, 22, 22, 23,
                 24, 24, 25, 27, 28, 28, 28, 28, 28, 30, 31))
  # which features are never selected
  expect_equal(which(rowSums(path_matrix) == 0), integer(0))  # zero seln
  expect_equal(rowSums(path_matrix),
               c(34, 51, 9, 33, 40, 2, 2, 20, 46, 32, 24, 13,
                 27, 33, 57, 57, 48, 33, 45, 8))
})
