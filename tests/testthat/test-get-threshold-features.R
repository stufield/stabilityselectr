
# Setup ----
pi_vec <- seq(0.6, 0.9, 0.1)


# Testing ----
test_that("`get_threshold_features()` returns the correct shape & dimensions", {
  feat <- get_threshold_features(ss$stabpath_matrix, thresh_vec = pi_vec)
  expect_type(feat, "list")
  expect_named(feat, paste0("thresh_", pi_vec))
  lapply(feat, expect_s3_class, "tbl_df")
  lapply(feat, expect_named, c("feature", "MaxSelectProb", "FDRbound"))
  expect_equal(lapply(feat, dim),
               list(thresh_0.6 = c(20, 3),
                    thresh_0.7 = c(20, 3),
                    thresh_0.8 = c(20, 3),
                    thresh_0.9 = c(5, 3)))
})

test_that("`get_threshold_features()` returns the correct values", {
  thresh_feat_tbl <- get_threshold_features(ss$stabpath_matrix,
                                            thresh_vec = pi_vec)
  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(thresh_feat_tbl) )
  )
})

test_that("`get_threshold_features()` returns NAs and trips warning when FDR <= 0.5", {
  expect_warning(
    thresh_feat_tbl <- get_threshold_features(
      ss$stabpath_matrix, warn = TRUE, thresh_vec = c(0.49, 0.51)),
    "FDR upper bound not defined for `thresh <= 0.5`"
  )
  expect_type(thresh_feat_tbl, "list")
  expect_true(all(vapply(thresh_feat_tbl, function(.x) nrow(.x) == 20, NA)))
  expect_true(all(is.na(thresh_feat_tbl$thresh_0.49$FDRbound)))
  expect_false(all(is.na(thresh_feat_tbl$thresh_0.51$FDRbound)))
  expect_equal(thresh_feat_tbl$thresh_0.51$FDRbound, seq(0.125, 2.5, by = 0.125))
})

test_that("`get_threshold_features()` trips warning when no stable features found", {
  expect_warning(
    thresh_feat_tbl <- get_threshold_features(
      ss$stabpath_matrix, warn = TRUE, thresh_vec = c(0.92, 0.97)),
    "No stable features at `thresh = 0.97`"
  )

  withr::with_options(list(pillar.sigfig = 4L),
    expect_snapshot( print(thresh_feat_tbl) )
  )

  # also test that same result if only 1 threshold passed
  tbl2 <- get_threshold_features(ss$stabpath_matrix, warn = FALSE,
                                 thresh_vec = 0.92)
  expect_equal(thresh_feat_tbl[-2L], tbl2)
})
