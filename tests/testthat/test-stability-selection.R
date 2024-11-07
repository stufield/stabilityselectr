
# Setup ----
# silence signalling from `signal_done()`
withr::local_options(list(signal.quiet = TRUE))

# Testing ----
# stability_selection ----
test_that("`stability_selection()` generates correct object", {
  expect_s3_class(ss, "stab_sel")
  expect_type(ss, "list")
  expect_true(is_stabpath_matrix(ss$stabpath_matrix))
  expect_length(ss, 15L)
  expect_true(is.matrix(ss$stabpath_matrix))
  expect_named(ss, c("stabpath_matrix", "lambda",
                     "alpha", "Pw",
                     "kernel", "num_iter",
                     "standardize", "impute_outliers",
                     "lambda_min_ratio", "perm_data",
                     "permpath_list", "perm_lambda",
                     "permpath_max", "beta", "r_seed")
    )
  expect_equal(ss$alpha, 0.8)
  expect_equal(ss$Pw, 0.5)
  expect_equal(ss$num_iter, 100L)
  expect_true(ss$standardize)
  expect_equal(ss$lambda_min_ratio, 0.1)
  expect_length(ss$lambda, 120L)
  expect_equal(ss$r_seed, 101)
  expect_equal(ss$alpha, 0.8)
  expect_equal(ss$kernel, "l1-logistic")
  expect_false(ss$perm_data)
  expect_null(ss$permpath_list)
  expect_s3_class(ss$permpath_max, "tbl_df")
  expect_equal(dim(ss$permpath_max), 0:1L)   # empty tibble; no perm
  expect_named(ss$permpath_max, "Feature")
  expect_false(ss$impute_outliers)
})

# Testing values ----
test_that("`stability_selection()` generates the correct values", {
  expect_equal(sum(ss$lambda), 5.87757157290224)
  expect_equal(sum(ss$perm_lambda), 1.00365738006801)
  expect_length(ss$perm_lambda, 20L)
  expect_equal(dim(ss$stabpath_matrix), c(20L, 120L))
  expect_equal(rownames(ss$stabpath_matrix),
               paste0("feat_", letters[1:20L]))
  expect_equal(diag(ss$stabpath_matrix),
               c(0.035, 0.015, 0.03, 0.08, 0.02, 0.04, 0.005, 0.005, 0.03,
                 0.17, 0.055, 0.105, 0.395, 0.015, 0.045, 0.035, 0.05, 0.25,
                 0.14, 0.375))
  expect_equal(sum(ss$stabpath_matrix), 1216.42)
  expect_equal(rowSums(ss$stabpath_matrix),
               c(64.96, 52.73, 60.875, 77.055, 54.44, 59.56, 53.285,
                 53.49, 49.225, 75.985, 55.3, 63.265, 86.24, 46.3,
                 47.75, 47.83, 56.51, 65.8, 66.645, 79.175),
               ignore_attr = TRUE
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
  expect_equal(rownames(ss1$stabpath_matrix), rownames(ss$stabpath_matrix))
})

# S3 summary method ----
test_that("the S3 `summary` method generates correct output", {
  expect_warning(
    summ <- summary(ss, thresh = 0.49, warn = TRUE),
    "FDR upper bound not defined for `thresh <= 0.5`"
  )   # thresh < 0.5 warning
  expect_true(is.na(sum(summ$FDRbound)))
  expect_equal(dim(summ), c(20, 4))
  summ <- summary(ss, thresh = 0.9)
  expect_s3_class(summ, "tbl_df")
  expect_equal(dim(summ), c(6, 4))
  out <- data.frame(stringsAsFactors = FALSE,
          feature     = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s", "feat_m"),
    MaxSelectProb     = c(0.945, 0.935, 0.915, 0.915, 0.905, 0.9),
                  AUC = c(0.422318758134717, 0.415582158399216,
                          0.342093188198169, 0.430051773514685,
                          0.343445346775382, 0.410396059422385),
             FDRbound = c(0.003125, 0.00625, 0.009375, 0.0125, 0.015625, 0.01875)
   )
  # make df because tibbles have weird floating point errors in testing
  expect_equal(out, data.frame(summ))
})

# S3 print method ----
test_that("`stab_sel` S3 print method returns expected known output", {
  # reinstate signalling for the print method
  withr::with_options(list(signal.quiet = FALSE), expect_snapshot_output(ss))
})

# NAs throw errors ----
test_that("`stability_selection()` with NAs throws an error", {
  x[sample(length(x), 50)] <- NA       # random NAs spread in data matrix
  expect_error(
    stability_selection(x, y),
    paste0("There are NAs detected in the data matrix ...\n",
           "This will cause `glmnet()` to fail ...\n",
           "Please remove the feature or use `splyr::imputeNAs()`"),
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

test_that("`stability_selection()` trips the correct errors if kernel = pca.sd", {
  expect_error(
    stability_selection(x, y, kernel = "pca.sd",
                       lambda_seq = lambda_vec,
                       alpha = 0.8, Pw = 0.5,   # alpha must be PC dim
                       standardize = TRUE),
    "Bad alpha: 0.8\nPlease set to number of principal components to consider\\.",
  )
})

# Cox kernel ----
test_that("`stability_selection()` generates expected values for the Cox kernel", {
  xcox <- strip_meta(log_rfu(sim_adat))
  ycox <- survival::Surv(sim_adat$time, sim_adat$status)
  ss_cox <- stability_selection(xcox, ycox, kernel = "Cox", r_seed = 101)

  expect_equal(sum(ss_cox$lambda), 13.533816327717)
  expect_equal(sum(ss_cox$perm_lambda), 2.3110419787012)
  expect_length(ss_cox$perm_lambda, 20L)
  expect_equal(dim(ss_cox$stabpath_matrix), c(40, 120))
  expect_equal(rownames(ss_cox$stabpath_matrix),
               c("seq.2802.68", "seq.9251.29", "seq.1942.70", "seq.5751.80",
                 "seq.9608.12", "seq.3459.49", "seq.3865.56", "seq.3363.21",
                 "seq.4487.88", "seq.5994.84", "seq.9011.72", "seq.2902.23",
                 "seq.2260.48", "seq.4936.96", "seq.2277.95", "seq.2953.31",
                 "seq.3032.11", "seq.4330.4",  "seq.4914.10", "seq.3896.5",
                 "seq.5002.7",  "seq.3476.4",  "seq.1130.49", "seq.6356.60",
                 "seq.4579.40", "seq.8344.24", "seq.8441.53", "seq.9360.55",
                 "seq.7841.8",  "seq.8142.63", "seq.4461.56", "seq.9297.97",
                 "seq.9396.38", "seq.3300.26", "seq.2772.14", "seq.6615.18",
                 "seq.8797.98", "seq.9879.88", "seq.8993.16", "seq.9373.82")
  )
  expect_equal(diag(ss_cox$stabpath_matrix), rep(0, 40))
  expect_equal(sum(ss_cox$stabpath_matrix), 29.285)
  expect_equal(rowSums(ss_cox$stabpath_matrix),
               c(0.1, 0.325, 0.275, 0.61, 0.11, 3.195, 0.035, 0.12,
                 0.155, 0, 0.355, 7.37, 2.925, 1.7, 0.645, 0.015,
                 3.95, 0, 0.12, 1.815, 0.19, 0.005, 0.065, 0, 0.185,
                 0.035, 0.075, 0.045, 0.025, 0.03, 0.03, 4.165, 0.02,
                 0, 0.05, 0.05, 0.43, 0, 0.01, 0.055),
               ignore_attr = TRUE
  )
})
