
# Setup ----
empFDR <- calc_emp_fdr_breaks(ss_perm, thresh_seq = seq(1, 0.1, -0.1))


# Testing ----
test_that("the `calc_emp_fdr_breaks()` function generates correct output", {
  expect_type(empFDR, "list")
  expect_length(empFDR, 2L)
  expect_named(empFDR, c("fdr_data", "breaks"))
  lapply(empFDR, expect_type, "list")
})

test_that("`calc_emp_fdr_breaks()` function generates `fdr_data` entry", {
  expect_s3_class(empFDR$fdr_data, "tbl_df")
  expect_named(empFDR$fdr_data, c("MeanFPs", "n_selected", "piThresh"))
  expect_equal(vapply(empFDR$fdr_data, typeof, ""),
               c(MeanFPs = "double", n_selected = "integer", piThresh  = "double"))
  expect_equal(colSums(empFDR$fdr_data),
               c(MeanFPs = 164.0, n_selected = 167.0, piThresh = 5.5))
  expect_equal(vapply(empFDR$fdr_data, median, pi),
               c(MeanFPs = 20.00, n_selected = 20.00, piThresh = 0.55))
  expect_equal(empFDR$fdr_data$piThresh, seq(1, 0.1, -0.1))
})

test_that("`calc_emp_fdr_breaks()` function generates `breaks` entry", {
  expect_snapshot( print(empFDR$breaks, n = Inf) )
})
