
# Setup ----
# silence signalling from `signal_done()`
withr::local_options(list(signal.quiet = TRUE))

# Testing ----
# stability_selection ----
test_that("`stability_selection()` generates correct object", {
  expect_s3_class(ss, "stab_sel")
  expect_type(ss, "list")
  expect_true(is.matrix(ss$stabpath_matrix))
  expect_true(is_stabpath_matrix(ss$stabpath_matrix))
  # snapshot the classes
  expect_snapshot( lapply(ss, class) )

  expect_length(ss, 15L)
  expect_equal(ss$alpha, 0.8)
  expect_equal(ss$Pw, 0.5)
  expect_equal(ss$n_iter, 100L)
  expect_true(ss$standardize)
  expect_equal(ss$lambda_min_ratio, 0.1)
  expect_length(ss$lambda, 120L)
  expect_equal(ss$r_seed, 101)
  expect_equal(ss$alpha, 0.8)
  expect_equal(ss$kernel, "l1-logistic")
  expect_false(ss$perm_data)
  expect_null(ss$permpath_list)
  expect_s3_class(ss$permpath_max, "tbl_df")
  expect_equal(dim(ss$permpath_max), 0:1L) # empty tibble; no perm
  expect_named(ss$permpath_max, "Feature")
  expect_false(ss$impute_outliers)
})

# Testing values ----
test_that("`stability_selection()` generates the correct values", {
  expect_equal(sum(ss$lambda), 5.8775715729)
  expect_equal(sum(ss$perm_lambda), 1.0036573800)
  expect_length(ss$perm_lambda, 20L)
  expect_equal(dim(ss$stabpath_matrix), c(ncol(x), length(ss$lambda)))
  expect_equal(rownames(ss$stabpath_matrix),
               paste0("feat_", letters[1:20L]))
  expect_equal(diag(ss$stabpath_matrix),
               c(0.035, 0.015, 0.03, 0.08, 0.02, 0.04, 0.005,
                 0.005, 0.03, 0.17, 0.055, 0.105, 0.395, 0.015,
                 0.045, 0.035, 0.05, 0.25, 0.14, 0.375))
  expect_equal(sum(ss$stabpath_matrix), 1216.42)

  withr::with_options(list(pillar.sigfig = 6L),
    expect_snapshot(
      print(
        tibble::enframe(rowSums(ss$stabpath_matrix), "feat", "rowsum"),
        n = Inf)
    )
  )
})

# Testing outlier imputation ----
test_that("the `stab_select` object is created correctly when imputing outliers", {
  ss1 <- stability_selection(x, y, kernel = "l1-logistic",
                            impute_outliers = TRUE, r_seed = 101)
  expect_s3_class(ss1, "stab_sel")
  expect_equal(dim(ss1$stabpath_matrix), dim(ss$stabpath_matrix))
  expect_true(ss1$impute_outliers)
  expect_equal(sum(ss1$stabpath_matrix), 1216.755)
  expect_equal(
    rownames(ss1$stabpath_matrix), rownames(ss$stabpath_matrix)
  )
})

# S3 summary method ----
test_that("the S3 `summary` method generates correct output", {
  withr::local_options(list(pillar.sigfig = 6L))
  expect_snapshot(
    summary(ss, thresh = 0.49, warn = TRUE)   # thresh < 0.5 warning
  )
  expect_snapshot( summary(ss, thresh = 0.9) )
})

# S3 print method ----
test_that("`stab_sel` S3 print method returns expected known output", {
  # reinstate signalling for the print method
  withr::with_options(
    list(signal.quiet = FALSE),
    expect_snapshot_output(ss)
  )
})

# NAs throw errors ----
test_that("`stability_selection()` with NAs throws an error", {
  x[sample(length(x), 50)] <- NA # insert random NAs in data matrix
  expect_error(
    stability_selection(x, y),
    "NAs detected in `x`, this will cause `glmnet()` to fail\nPlease check: ",
    fixed = TRUE
  )
})

# Trips correct errors ----
test_that("`stability_selection()` trips the correct errors if kernel = ridge", {
  expect_error(
    stability_selection(x, y, kernel = "ridge",
                       lambda_seq = lambda_vec,
                       elastic_alpha = 1L,     # ridge is alpha = 0!
                       standardize = TRUE),
    paste("Invalid `elastic_alpha =` argument for 'ridge' kernels (1).",
          "Please set `elastic_alpha = 0`."), fixed = TRUE
  )
})

test_that("`stability_selection()` trips the correct errors if kernel = pca.thresh", {
  expect_error(
    stability_selection(x, y, kernel = "pca.thresh",
                       lambda_seq = lambda_vec,
                       alpha = 0.8, Pw = 0.5,   # alpha must be PC dim
                       standardize = TRUE),
    "Bad alpha: 0.8\nPlease set to number of principal components to consider\\."
  )
})

test_that("`stability_selection()` trips the correct error for `kernel = pca.sd`", {
  expect_error(
    stability_selection(x, y, kernel = "pca.sd",
                       lambda_seq = lambda_vec,
                       alpha = 0.8, Pw = 0.5, # alpha must be PC dim
                       standardize = TRUE),
    "Bad alpha: 0.8\nPlease set to number of principal components to consider\\.",
  )
})

# Cox kernel ----
test_that("`stability_selection()` generates expected values for the Cox kernel", {
  xcox   <- feature_matrix(log_rfu(simdata))
  ycox   <- survival::Surv(simdata$time, simdata$status)  # a matrix
  ss_cox <- stability_selection(xcox, ycox, kernel = "cox", r_seed = 101)

  expect_equal(sum(ss_cox$lambda), 13.5338163)
  expect_equal(sum(ss_cox$perm_lambda), 2.31104197)
  expect_length(ss_cox$perm_lambda, 20L)
  expect_equal(dim(ss_cox$stabpath_matrix),
               c(ncol(xcox), length(ss_cox$lambda)))
  expect_equal(sum(ss_cox$stabpath_matrix), 26.595)

  withr::with_options(list(pillar.sigfig = 6L),
    expect_snapshot(
      print(
        tibble::enframe(rowSums(ss_cox$stabpath_matrix), "feat", "rowsum"),
        n = Inf)
    )
  )
})
