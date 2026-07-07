
# Setup ----
# ss_perm from setup.R
pi_vec <- seq(0.5, 0.9, 0.1)
fdr    <- calc_emp_fdr(ss_perm, thresh_seq = pi_vec)


# Testing ----
# Shapes ----
test_that("`calc_emp_fdr()` returns the correct shape & dimensions", {
  expect_null(dim(fdr))
  expect_type(fdr, "double")
  expect_length(fdr, length(pi_vec))
  expect_named(fdr, c("thresh_0.5", "thresh_0.6",
                      "thresh_0.7", "thresh_0.8", "thresh_0.9"))
})

# Values ----
test_that("`calc_emp_fdr()` returns the correct values", {
  expect_equal(
    fdr,
    c(thresh_0.5 = 20.00, thresh_0.6 = 20.00, thresh_0.7 = 20.0,
      thresh_0.8 = 19.2, thresh_0.9 = 4.80)
  )
})

# Threshold @ 1 value ----
test_that("`calc_emp_fdr()` returns correct values when thresh a double", {
  expect_equal(
    calc_emp_fdr(ss_perm, thresh_seq = 0.85), c("thresh_0.85" = 12.5)
  )
})

# Trips a warning < 5 perms ----
test_that("`calc_emp_fdr()` trips a warning when < 5 permutations", {
  ss2 <- ss_perm
  ss2$permpath_list <- head(ss2$permpath_list, 4L)
  expect_message(
    ssW <- calc_emp_fdr(ss2, thresh_seq = seq(0.1, 0.9, 0.1)),
    "Returning mean of < 5 permutations"
  )
  expect_equal(
    ssW, c(thresh_0.1 = 20, thresh_0.2 = 20, thresh_0.3 = 20,
           thresh_0.4 = 20, thresh_0.5 = 20, thresh_0.6 = 20,
           thresh_0.7 = 20, thresh_0.8 = 19.50, thresh_0.9 = 4.25)
  )
})
