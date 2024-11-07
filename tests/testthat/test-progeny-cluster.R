
# Setup ----
# iter= and size= are passed to progenyK() via the '...'
pclust <- withr::with_seed(101,
  progeny_cluster(progeny_data, clust_iter = 2:6L,
                  reps = 5L, iter = 10L, size = 6)
)


# Testing  ----
test_that("the object returned by `progeny_cluster()` is correct", {
  expect_s3_class(pclust, "pclust")
  expect_type(pclust, "list")
  expect_named(pclust, c("scores", "mean_scores",
                         "ci95_scores", "random_scores",
                         "mean_random_scores", "D_max", "D_gap",
                         "clust_iter", "repeats", "iter", "size", "call"))
})

test_that("the individual outputs are currect for each element of pclust", {
  expect_true(is.matrix(pclust$scores))
  expect_equal(dim(pclust$scores), c(5, 5))
  expect_equal(sum(pclust$scores), 393.9343373)
  expect_equal(pclust$mean_scores, apply(pclust$scores, 2, mean))

  expect_true(is.matrix(pclust$ci95_scores))
  expect_equal(dim(pclust$ci95_scores), c(2L, 5L))
  expect_equal(apply(pclust$ci95_scores, 1, mean),
               c("2.5%"  = 9.812698941, "97.5%" = 25.214726238))
})

test_that("`is_pclust()` returns correct logical", {
  expect_true(is_pclust(pclust))
  expect_false(is_pclust(unclass(pclust)))
})

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

test_that("the `iter =` argument is passed via the ..., otherwise throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 2:6L),
    "You must pass an `iter =` argument via the `...` to `progeny_cluster()`.",
    fixed = TRUE
  )
})

test_that("the `iter =` argument is >= 1 ..., otherwise throw error", {
  expect_error(
    progeny_cluster(progeny_data, clust_iter = 2:6L, reps = 5L, iter = 0L),
    paste("The number of iterations can't be zero or negative.",
          "Please use a positive value.")
  )
})

test_that("`progeny_cluster()` can handle data frames as well as matrices", {
  p <- withr::with_seed(101,
    progeny_cluster(data.frame(progeny_data),
                    clust_iter = 2:6L, reps = 5L, iter = 10L, size = 6)
  )
  expect_equal(discard_it(p, names(p) == "call"),
               discard_it(pclust, names(pclust) == "call"))
})

test_that("`progeny_k()` trips appropriate errors as expected", {
  expect_error(
    progeny_k(head(progeny_data, 15), k = 2, size = 18, iter = 10),
    paste("You are probably progeny sampling with too many samples for",
          "these data ... try a smaller number\\.")
  )
  expect_error(
    progeny_k(head(progeny_data, 10), k = 4, size = 4, iter = 10),
    paste("You are probably progeny sampling with too many",
          "samples ... perhaps try `size` < 4\\."),
  )
})

# S3 print method ----
test_that("`pclust` S3 print method returns expected known output", {
  expect_snapshot_output(pclust)
})
