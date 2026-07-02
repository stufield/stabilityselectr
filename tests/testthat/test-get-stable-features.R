
# Setup ----
pi_vec <- seq(0.6, 0.9, 0.1)

# Testing ----
# No permutation ----
test_that("`get_stable_features()` without permutation check shape & dims", {
  s_feat <- get_stable_features(ss, 0.55)$thresh_0.55
  expect_s3_class(s_feat, "tbl_df")
  expect_equal(dim(s_feat), c(n_feat, 3L))       # all feats included
  expect_named(s_feat, c("feature", "MaxSelectProb", "FDRbound"))
  expect_setequal(s_feat$feature, colnames(x))
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
  s_feat <- get_stable_features(ss, 0.55)$thresh_0.55

  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(s_feat, n = Inf) )
  )

  expect_true(all(s_feat$MaxSelectProb <= 1L))
  expect_true(all(s_feat$MaxSelectProb >= 0L))
  expect_true(all(s_feat$FDRbound <= 1L))
  expect_true(all(s_feat$FDRbound >= 0L))
  expect_setequal(colnames(x), s_feat$feature)
})

test_that("`get_stable_features()` returns correct values at thresh <= 0.5", {
  expect_warning(
    s_feat <- get_stable_features(ss, 0.49, warn = TRUE)$thresh_0.49,
    "FDR upper bound not defined for `thresh <= 0.5`"
  )
  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(s_feat, n = Inf) )
  )
})

test_that("`get_stable_features` returns correct values at thresh = 0.90", {
  s_feat <- get_stable_features(ss, 0.9)$thresh_0.9
  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(s_feat, n = Inf) )
  )
  expect_true(all(s_feat$MaxSelectProb <= 1L))
  expect_true(all(s_feat$MaxSelectProb >= 0L))
  expect_true(all(s_feat$FDRbound <= 1L))
  expect_true(all(s_feat$FDRbound >= 0L))
  expect_false(all(colnames(x) %in% s_feat$feature))
  expect_true(all(s_feat$feature %in% colnames(x)))
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
    s_feat <- get_stable_features(ss, 0.99, warn = TRUE)$thresh_0.99,
    "No stable features at `thresh = 0.99`"
  )
  expect_equal(dim(s_feat), c(0, 3))
  expect_named(s_feat, c("feature", "MaxSelectProb", "FDRbound"))
})

# with permutation ----
test_that("`get_stable_features()` with permutation adds `EmpFDR` column", {
  s_feat_perm <- get_stable_features(ss_perm, 0.88)$thresh_0.88

  expect_s3_class(s_feat_perm, "tbl_df")
  expect_false(all(colnames(x) %in% s_feat_perm$feature))

  # un-selected features
  expect_equal(setdiff(colnames(x), s_feat_perm$feature),
               c("feat_b", "feat_e", "feat_h", "feat_i", "feat_j",
                 "feat_k", "feat_n", "feat_o", "feat_p", "feat_q"))

  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(s_feat_perm, n = Inf) )
  )
})


test_that("`get_stable_features()` accepts a vector of thresholds", {
  thresh_feat_tbl <- get_stable_features(ss$stabpath_matrix, thresh = pi_vec)
  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(thresh_feat_tbl) )
  )
})

test_that("`get_stable_features()` trips warning when no stable features found", {
  expect_warning(
    thresh_feat_tbl <- get_stable_features(
      ss$stabpath_matrix, warn = TRUE, thresh = c(0.92, 0.97)),
    "No stable features at `thresh = 0.97`"
  )

  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(thresh_feat_tbl) )
  )
})
