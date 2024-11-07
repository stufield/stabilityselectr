
# Setup ----
pi_vec <- seq(0.6, 0.9, 0.1)


# Testing ----
test_that("`get_threshold_features()` returns the correct shape & dimensions", {
  feat <- get_threshold_features(ss$stabpath_matrix, thresh_vec = pi_vec)
  expect_type(feat, "list")
  expect_named(feat, paste0("thresh_", pi_vec))
  lapply(feat, expect_s3_class, "data.frame")
  lapply(feat, expect_named, c("MaxSelectProb", "FDRbound"))
  expect_equal(lapply(feat, dim),
               list(thresh_0.6 = c(20, 2),
                    thresh_0.7 = c(20, 2),
                    thresh_0.8 = c(19, 2),
                    thresh_0.9 = c(6, 2)))
})

test_that("`get_threshold_features()` returns the correct values", {
  feat <- get_threshold_features(ss$stabpath_matrix, thresh_vec = pi_vec)
  expect_equal(lapply(feat, colSums),
               list(thresh_0.6 = c(MaxSelectProb = 17.460, FDRbound = 2.62500),
                    thresh_0.7 = c(MaxSelectProb = 17.460, FDRbound = 1.3125),
                    thresh_0.8 = c(MaxSelectProb = 16.675000, FDRbound = 0.7916666666),
                    thresh_0.9 = c(MaxSelectProb = 5.515000, FDRbound = 0.065625)))

  expect_equal(lapply(feat, rownames),
    list(thresh_0.6 = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s",
                        "feat_m", "feat_l", "feat_f", "feat_q", "feat_g",
                        "feat_e", "feat_n", "feat_r", "feat_c", "feat_k",
                        "feat_b", "feat_o", "feat_h", "feat_p", "feat_i"),
         thresh_0.7 = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s",
                        "feat_m", "feat_l", "feat_f", "feat_q", "feat_g",
                        "feat_e", "feat_n", "feat_r", "feat_c", "feat_k",
                        "feat_b", "feat_o", "feat_h", "feat_p", "feat_i"),
        thresh_0.8 = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s",
                       "feat_m", "feat_l", "feat_f", "feat_q", "feat_g",
                       "feat_e", "feat_n", "feat_r", "feat_c", "feat_k",
                       "feat_b", "feat_o", "feat_h", "feat_p"),
         thresh_0.9 = c("feat_d", "feat_t", "feat_a", "feat_j", "feat_s", "feat_m"))
  )
  # picked this one at random because it's small enough to check
  expect_equal(feat$thresh_0.9,
               data.frame(row.names = c("feat_d", "feat_t", "feat_a",
                                        "feat_j", "feat_s", "feat_m"),
                          MaxSelectProb = c(0.945, 0.935, 0.915, 0.915,
                                            0.905, 0.9),
                          FDRbound = c(0.003125, 0.00625, 0.009375,
                                       0.0125, 0.015625, 0.01875)
               )
   )
})

test_that("`get_threshold_features()` returns NAs and trips warning when FDR <= 0.5", {
  expect_warning(
    feat <- get_threshold_features(ss$stabpath_matrix, warn = TRUE,
                                  thresh_vec = c(0.49, 0.51)),
    "FDR upper bound not defined for `thresh <= 0.5`"
  )
  expect_type(feat, "list")
  expect_true(all(vapply(feat, function(.x) nrow(.x) == 20, NA)))
  expect_true(all(is.na(feat$thresh_0.49$FDRbound)))
  expect_false(all(is.na(feat$thresh_0.51$FDRbound)))
  expect_equal(feat$thresh_0.51$FDRbound, seq(0.125, 2.5, by = 0.125))
})

test_that("`get_threshold_features()` trips warning when no stable features found", {
  expect_warning(
    feat <- get_threshold_features(ss$stabpath_matrix, warn = TRUE,
                                  thresh_vec = c(0.93, 0.95)),
    "No stable features at `thresh = 0.95`"
  )
  expect_type(feat, "list")
  expect_equal(feat$thresh_0.93,
               data.frame(MaxSelectProb = c(0.945, 0.935),
                          FDRbound      = c(0.00290697674418605,
                                            0.00581395348837209),
                          row.names     = c("feat_d", "feat_t")
                          )
  )
  expect_equal(feat$thresh_0.95,
               data.frame(MaxSelectProb = numeric(0),
                          FDRbound      = numeric(0),
                          row.names     = character(0)))
})

test_that("`get_threshold_features()` returns same if only 1 threshold passed", {
  a <- get_threshold_features(ss$stabpath_matrix, 0.9)
  b <- get_threshold_features(ss$stabpath_matrix, 0.91)
  expect_type(a, "list")
  expect_type(b, "list")
  expect_equal(a$thresh_0.9,
               data.frame(row.names = c("feat_d", "feat_t",
                                        "feat_a", "feat_j",
                                        "feat_s", "feat_m"),
                              MaxSelectProb = c(0.945, 0.935, 0.915,
                                                0.915, 0.905, 0.9),
                          FDRbound = c(0.003125, 0.00625, 0.009375,
                                       0.0125, 0.015625, 0.01875)
                            )
  )
  expect_equal(b$thresh_0.91,
               data.frame(row.names = c("feat_d", "feat_t", "feat_a", "feat_j"),
                          MaxSelectProb = c(0.945, 0.935, 0.915, 0.915),
                          FDRbound = c(0.00304878048780488,
                                       0.00609756097560976,
                                       0.00914634146341463,
                                       0.0121951219512195)
               )
   )
})
