
# Setup ----
# `n_iter =` and `size =` are passed to `progeny_k()` via the '...'
pclust <- withr::with_seed(101,
  progeny_cluster(progeny_data, clust_iter = 2:6L,
                  reps = 5L, n_iter = 10L, size = 6)
)


# Testing  ----
test_that("the object returned by `progeny_cluster()` is correct", {
  expect_s3_class(pclust, "pclust")
  expect_true(is_pclust(pclust))   # check the helper
  expect_type(pclust, "list")
  expect_named(pclust, c("scores", "mean_scores",
                         "ci95_scores", "random_scores",
                         "mean_random_scores", "D_max",
                         "D_gap", "clust_iter",
                         "repeats", "n_iter",
                         "size", "call"))
  expect_snapshot( lapply(pclust, class) )
})

test_that("the object elements are the numerically expected results", {
  skip_on_os("linux")
  expect_snapshot( unclass(pclust) )
})

test_that("`progeny_cluster()` can handle data frames as well as matrices", {
  p <- withr::with_seed(101,
    progeny_cluster(data.frame(progeny_data),
                    clust_iter = 2:6L, reps = 5L, n_iter = 10L, size = 6)
  )
  expect_equal(
    discard_it(p, names(p) == "call"),
    discard_it(pclust, names(pclust) == "call")
  )
})

# S3 print method ----
test_that("`pclust` S3 print method returns expected known output", {
  skip_on_os("linux")
  expect_snapshot_output(pclust)
})

# error handling ----
test_that("cluster comparisons is at least 3 ... throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 2:3L,
    "The number of clusters for comparison is too small: 2")
  )
})

test_that("clusters are all > 2 ... throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 1:5L),
    "The number of clusters must be greater than 1."
  )
})

test_that("the `reps =` argument is >= 1 ..., otherwise throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 2:6L, reps = 0L),
    paste("The number of repeats can't be zero or negative.",
          "Please use a positive value.")
  )
})

test_that("the `n_iter =` argument is passed via the `...`; throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 2:6L),
    "You must pass an `n_iter =` argument via the `...` to `progeny_cluster()`.",
    fixed = TRUE
  )
})

test_that("the `n_iter =` argument is >= 1 `...`; throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 2:6L, reps = 5L, n_iter = 0L),
    paste("The number of iterations can't be zero or negative.",
          "Please use a positive value.")
  )
})

test_that("`progeny_k()` trips appropriate errors as expected", {
  expect_error(
    progeny_k(head(progeny_data, 15), k = 2, size = 18, n_iter = 10),
    paste("You are probably progeny sampling with too many samples for",
          "these data ... try a smaller number")
  )
  expect_error(
    progeny_k(head(progeny_data, 10), k = 4, size = 4, n_iter = 10),
    paste("You are probably progeny sampling with too many",
          "samples ... perhaps try `size` < 4"),
  )
})
