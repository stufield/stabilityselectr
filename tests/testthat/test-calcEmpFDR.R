
# Setup ----
# ss_perm10 from helper.R
pi_vec <- seq(0.5, 0.9, 0.1)
fdr    <- calcEmpFDR(ss_perm10, thresh.seq = pi_vec)


# Testing ----
# Shapes ----
test_that("`calcEmpFDR()` returns the correct shape & dimensions", {
  expect_null(dim(fdr))
  expect_type(fdr, "double")
  expect_length(fdr, length(pi_vec))
  expect_named(fdr, c("thresh_0.5", "thresh_0.6",
                      "thresh_0.7", "thresh_0.8", "thresh_0.9"))
})

# Values ----
test_that("`calcEmpFDR()` returns the correct values", {
  expect_equal(
    fdr,
    c(thresh_0.5 = 20.00, thresh_0.6 = 20.00, thresh_0.7 = 19.9,
      thresh_0.8 = 18.7, thresh_0.9 = 3.90)
  )
})

# Threshold 1 value ----
test_that("`calcEmpFDR()` returns correct values when thresh is 1 value", {
  pi75 <- calcEmpFDR(ss_perm10, thresh.seq = 0.75)
  expect_type(pi75, "double")
  expect_length(pi75, 1L)
  expect_equal(pi75, c(thresh_0.75 = 19.90))
})

# Trips a warning < 5 perms ----
test_that("`calcEmpFDR()` trips a warning when < 5 permutations", {
  ss2 <- ss_perm10
  ss2$permpath.list <- head(ss2$permpath.list, 4L)
  expect_warning(
    ssW <- calcEmpFDR(ss2, thresh.seq = seq(0.1, 0.9, 0.1)),
    "Returning mean of < 5 permutations"
  )
  expect_type(ssW, "double")
  expect_length(ssW, 9L)
  expect_equal(ssW, c(thresh_0.1 = 20,
                      thresh_0.2 = 20,
                      thresh_0.3 = 20,
                      thresh_0.4 = 20,
                      thresh_0.5 = 20,
                      thresh_0.6 = 20,
                      thresh_0.7 = 19.75,
                      thresh_0.8 = 18.50,
                      thresh_0.9 = 4.75))
})
