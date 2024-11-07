
# Setup ----

# Testing ----
# No permutation ----
test_that("`get_stable_features()` without permutation check shape & dims", {
  s_feat <- get_stable_features(ss, 0.55)
  expect_s3_class(s_feat, "data.frame")
  expect_equal(dim(s_feat), c(n_feat, 2))       # all feats included
  expect_named(s_feat, c("MaxSelectProb", "FDRbound"))
  expect_setequal(rownames(s_feat), colnames(x))
})

# S3 method ----
test_that("the S3 methods work for `get_stable_features()`", {
  expect_equal(
    get_stable_features(ss, 0.55),
    get_stable_features(ss$stabpath_matrix, 0.55)
  )
  # default method; nothing for integer
  expect_error(
    get_stable_features(1:5L, 0.55),
    "Couldn't find a S3 method for this class object: 'integer'"
  )
})

test_that("`get_stable_features()` returns correct values at thresh = 0.55", {
  s_feat <- get_stable_features(ss, 0.55)
  expect_equal(
    s_feat,
    data.frame(
      row.names = c("feat_d", "feat_t", "feat_a",
                   "feat_j", "feat_s", "feat_m", "feat_l", "feat_f", "feat_q",
                   "feat_g", "feat_e", "feat_n", "feat_r", "feat_c",
                   "feat_k", "feat_b", "feat_o", "feat_h", "feat_p", "feat_i"),
      MaxSelectProb = c(0.945, 0.935, 0.915, 0.915, 0.905, 0.900, 0.890,
                        0.885, 0.880, 0.870, 0.860, 0.860, 0.860, 0.855, 0.850,
                        0.845, 0.840, 0.835, 0.830, 0.785),
      FDRbound = c(0.025, 0.05, 0.0749999999999999,
                   0.0999999999999999, 0.125, 0.15, 0.175, 0.2, 0.225,
                   0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45,
                   0.475, 0.5)
    )
  )
  expect_true(all(s_feat$MaxSelectProb <= 1L))
  expect_true(all(s_feat$MaxSelectProb >= 0L))
  expect_true(all(s_feat$FDRbound <= 1L))
  expect_true(all(s_feat$FDRbound >= 0L))
  expect_setequal(colnames(x), rownames(s_feat))
})

test_that("`get_stable_features()` returns correct values at thresh <= 0.5", {
  expect_warning(
    s_feat <- get_stable_features(ss, 0.49, warn = TRUE),
    "FDR upper bound not defined for `thresh <= 0.5`"
  )
  expect_equal(
    s_feat,
    data.frame(row.names = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s",
                             "feat_m", "feat_l", "feat_f", "feat_q", "feat_g",
                             "feat_e", "feat_n", "feat_r", "feat_c", "feat_k",
                             "feat_b", "feat_o", "feat_h", "feat_p", "feat_i"),
               MaxSelectProb = c(0.945, 0.935, 0.915, 0.915, 0.905, 0.900,
                                 0.890, 0.885, 0.880, 0.870, 0.860, 0.860,
                                 0.860, 0.855, 0.850, 0.845, 0.840, 0.835,
                                 0.830, 0.785),
               FDRbound = rep(NA_real_, n_feat)
    )
  )
})

test_that("`get_stable_features` returns correct values at thresh = 0.90", {
  s_feat <- get_stable_features(ss, 0.90)
  expect_equal(
    s_feat,
    data.frame(
      row.names = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s", "feat_m"),
      MaxSelectProb = c(0.945, 0.935, 0.915, 0.915, 0.905, 0.900),
      FDRbound = c(0.003125, 0.00625, 0.009375, 0.01250, 0.015625, 0.018750)
    )
  )
  expect_true(all(s_feat$MaxSelectProb <= 1L))
  expect_true(all(s_feat$MaxSelectProb >= 0L))
  expect_true(all(s_feat$FDRbound <= 1L))
  expect_true(all(s_feat$FDRbound >= 0L))
  expect_false(all(colnames(x) %in% rownames(s_feat)))
  expect_true(all(rownames(s_feat) %in% colnames(x)))
})

# Trips warnings ----
test_that("the `get_stable_features()` trips proper warnings", {
  expect_warning(
    get_stable_features(ss, thresh = 0.49, warn = TRUE),   # thresh < 0.5 warning
    "FDR upper bound not defined for `thresh <= 0.5`"
  )
  expect_warning(
    get_stable_features(ss, thresh = 0.25, warn = TRUE),   # thresh < 0.5 warning
    "FDR upper bound not defined for `thresh <= 0.5`"
  )
  expect_warning(
    get_stable_features(ss, thresh = 0.97, warn = TRUE),   # no stable features
    "No stable features at `thresh = 0.97`"
  )
  expect_warning(
    get_stable_features(ss, thresh = 0.98, warn = TRUE),   # no stable features
    "No stable features at `thresh = 0.98`"
  )
  expect_warning(
    get_stable_features(ss, thresh = 0.99, warn = TRUE),   # no stable features
    "No stable features at `thresh = 0.99`"
  )
})

# Empty when thresh high ----
test_that("`get_stable_features()` returns empty data frame when thresh too high", {
  expect_warning(
    s_feat <- get_stable_features(ss, 0.99, warn = TRUE),
    "No stable features at `thresh = 0.99`"
  )
  expect_equal(dim(s_feat), c(0, 2))
  expect_named(s_feat, c("MaxSelectProb", "FDRbound"))
})

# with permutation ----
test_that("`get_stable_features()` with permutation adds `EmpFDR` column", {
  s_feat_perm <- get_stable_features(ss_perm10, 0.85)
  expect_s3_class(s_feat_perm, "data.frame")
  expect_false(all(colnames(x) %in% rownames(s_feat_perm)))
  expect_true(all(rownames(s_feat_perm) %in% colnames(x)))

  # un-selected features
  expect_equal(setdiff(colnames(x), rownames(s_feat_perm)),
               c("feat_b", "feat_g", "feat_h", "feat_i", "feat_k",
                 "feat_l", "feat_n", "feat_o", "feat_p", "feat_q"))
  expect_equal(s_feat_perm,
               data.frame(row.names = c("feat_s", "feat_t", "feat_d", "feat_m",
                                        "feat_j", "feat_r", "feat_a", "feat_c",
                                        "feat_e", "feat_f"),
                         MaxSelectProb = c(0.955, 0.955, 0.925, 0.925, 0.89,
                                           0.87, 0.865, 0.865, 0.86, 0.86),
                              FDRbound = c(0.00357142857142857, 0.00714285714285714,
                                           0.0107142857142857, 0.0142857142857143,
                                           0.0178571428571429,
                                           0.0214285714285714, 0.025, 0.0285714285714286,
                                           0.0321428571428571, 0.0357142857142857),
                                EmpFDR = c(0.13, 0.13, 0.24, 0.24, 0.54, 0.85,
                                           0.95, 0.95, 1.02, 1.02)
                           )
     )
})
