
test_that("the `clust_data` object is correct", {
  target <- data.frame(F1 = c(0.4, 0.98, 0.35, 0.62, 0.48, 0.63, 0.22, 0.92, 0.38, 0.05,
                              0.91, 0.57, 0.77, 0.41, 0.57, 0.02, 0.19, 0.63, 0.95, 0.16),
                       F2 = c(0.85, 0.24, 0.98, 0.08, 0.79, 0.38, 0.88, 0.14, 0.75, 0.56,
                              0.33, 0.32, 0.22, 0.75, 0.35, 0.61, 0.63, 0.34, 0.09, 0.67))
  expect_equal(clust_data, target)
})

test_that("the `progeny_data` object is correct", {
  expect_true(is.matrix(progeny_data))
  expect_equal(dim(progeny_data), c(150L, 2L))
  expect_equal(colSums(progeny_data), c(F1 = -1.06695922710, F2 = -5.45820080843))
  expect_equal(median(rowSums(progeny_data)), 0.569218705589)
})
