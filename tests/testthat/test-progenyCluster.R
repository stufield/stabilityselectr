
# Setup ----
# iter= and size= are passed to progenyK() via the '...'
pclust <- withr::with_seed(101,
  progenyCluster(progeny_data, clust.iter = 2:6, reps = 5, iter = 10, size = 6)
)


# Testing  ----
test_that("the object returned by `progenyCluster()` is correct", {
  expect_s3_class(pclust, "pclust")
  expect_type(pclust, "list")
  expect_named(pclust, c("scores", "mean.scores",
                         "ci95.scores", "random.scores",
                         "mean.random.scores", "D.max", "D.gap",
                         "clust.iter", "repeats", "iter", "size", "call"))
})

test_that("the individual outputs are currect for each element of pclust", {
  expect_true(is.matrix(pclust$scores))
  expect_equal(dim(pclust$scores), c(5, 5))
  expect_equal(sum(pclust$scores), 393.9343373)
  expect_equal(pclust$mean.scores, apply(pclust$scores, 2, mean))

  expect_true(is.matrix(pclust$ci95.scores))
  expect_equal(dim(pclust$ci95.scores), c(2, 5))
  expect_equal(apply(pclust$ci95.scores, 1, mean),
               c("2.5%"  = 9.812698941, "97.5%" = 25.214726238))
})

test_that("`is.pclust()` returns correct logical", {
  expect_true(is.pclust(pclust))
  expect_false(is.pclust(unclass(pclust)))
})

test_that("cluster comparisons is at least 3 ... throw error", {
  expect_error(
    progenyCluster(progeny_data, clust.iter = 2:3),
    "The number of clusters for comparison is too small: 2"
  )
})

test_that("clusters are all > 2 ... throw error", {
  expect_error(
    progenyCluster(progeny_data, clust.iter = 1:5),
    "The number of clusters must be greater than 1."
  )
})

test_that("the `reps =` argument is >= 1 ..., otherwise throw error", {
  expect_error(
    progenyCluster(progeny_data, clust.iter = 2:6, reps = 0),
    paste("The number of repeats can't be zero or negative.",
          "Please use a positive value.")
  )
})

test_that("the `iter =` argument is passed via the ..., otherwise throw error", {
  expect_error(
    progenyCluster(progeny_data, clust.iter = 2:6),
    "You must pass an `iter =` argument via the `...` to `progenyCluster()`.",
    fixed = TRUE
  )
})

test_that("the `iter =` argument is >= 1 ..., otherwise throw error", {
  expect_error(
    progenyCluster(progeny_data, clust.iter = 2:6, reps = 5, iter = 0),
    paste("The number of iterations can't be zero or negative.",
          "Please use a positive value.")
  )
})

test_that("`progenyCluster()` can handle data frames as well as matrices", {
  p <- withr::with_seed(101,
    progenyCluster(data.frame(progeny_data),
                   clust.iter = 2:6, reps = 5, iter = 10, size = 6)
  )
  expect_equal(discard_it(p, names(p) == "call"),
               discard_it(pclust, names(pclust) == "call"))
})

test_that("`progenyK()` trips appropriate errors as expected", {
  expect_error(
    progenyK(head(progeny_data, 15), k = 2, size = 18, iter = 10),
    paste("You are probably progeny sampling with too many samples for",
          "these data ... try a smaller number\\.")
  )
  expect_error(
    progenyK(head(progeny_data, 10), k = 4, size = 4, iter = 10),
    paste("You are probably progeny sampling with too many",
          "samples ... perhaps try `size` < 4\\."),
  )
})

# S3 print method ----
test_that("`pclust` S3 print method returns expected known output", {
  expect_snapshot_output(pclust)
})
