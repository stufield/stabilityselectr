#' Perform Progeny Clustering
#'
#' Determine the most stable (optimal) number of clusters via
#'   Progeny Clustering algorithm.
#'
#' @family cluster
#'
#' @param data A (`n x p`) data matrix containing *n* samples and *p* features.
#'   Can also be a data frame where each row corresponds to a sample or
#'   observation, whereas each column corresponds to a feature or variable.
#' @param clust_iter `integer(n)`. Span of `k` clusters to interrogate
#' @param reps `integer(1)`. The number of repeat iterations to perform.
#'   Particularly useful if error bars during plotting are desired.
#' @param verbose `logical(1)`. Print the progress of the clustering repeats
#'   to the console. Defaults to [interactive()].
#' @param ... Additional parameters passed to the internal `progenyK()`,
#'   typically `iter =` and `size =`.
#' @return An object of class `pclust`, a list containing:
#'   \item{scores:}{A matrix of stability scores for each iteration in a matrix,
#'      with `k` columns}
#'   \item{mean_scores:}{The mean stability scores for each cluster `k`}
#'   \item{ci95_scores:}{95% confidence interval scores}
#'   \item{random_scores:}{The reference (random) scores for each iteration
#'     at each clustering level (`k`)}
#'   \item{mean_random_scores:}{The mean of the reference (random) data set, i.e.
#'     column means of `random.scores`}
#'   \item{D_max:}{The distance between the mean stability scores and the mean
#'     reference scores for each cluster `k`}
#'   \item{D_gap:}{The "gap" distance metric for neighboring cluster k differences.
#'     See original paper for reference.}
#'   \item{clust_iter:}{Integer Sequence of `k` clusters interrogated}
#'   \item{repeats:}{The number of repeat iterations to performed}
#'   \item{iter:}{The number of progeny iterations to performed}
#'   \item{size:}{The progeny size used in each iteration}
#'   \item{call:}{The call made to [progeny_cluster()]}
#'
#' @author Stu Field
#' @references Hu, C.W., Kornblau, S.M., Slater, J.H. and A.A. Qutub (2015).
#'   Progeny Clustering: A Method to Identify Biological Phenotypes.
#'   Scientific Reports, 5:12894. \url{http://www.nature.com/articles/srep12894}
#' @seealso [stats::kmeans()]
#'
#' @examples
#' # `iter=` and `size=` are passed to `progeny_k()`
#' pclust <- withr::with_seed(1234,
#'   progeny_cluster(progeny_data, clust_iter = 2:9L, iter = 20L, size = 6)
#' )
#'
#' # Test progeny clustering on Iris data set
#' # Doesn't work quite as well as the simulated data set
#' clustIris <- withr::with_seed(99,
#'   progeny_cluster(iris[, -5L], clust_iter = 2:5L, size = 6L, iter = 50)
#' )
#' clustIris    # true n clusters = 3
#'
#' @importFrom stats quantile runif sd
#' @export
progeny_cluster <- function(data, clust_iter = 2:10L, reps = 10L,
                            verbose = interactive(), ...) {

  if ( length(clust_iter) < 3L ) {
    stop(
      "The number of clusters for comparison is too small: ",
      value(length(clust_iter)), call. = FALSE
    )
  }
  if ( sum(clust_iter < 2L) != 0 ) {
    stop("The number of clusters must be greater than 1.", call. = FALSE)
  }
  if ( reps < 1L ) {
    stop(
      "The number of repeats can't be zero or negative. ",
      "Please use a positive value.", call. = FALSE
    )
  }
  if ( !"iter" %in% names(list(...)) ) {
    stop(
      "You must pass an `iter =` argument via the `...` to `progeny_cluster()`.",
      call. = FALSE
    )
  }
  if ( list(...)$iter < 1 ) {
    stop(
      "The number of iterations can't be zero or negative. ",
      "Please use a positive value.", call. = FALSE
    )
  }
  if ( inherits(data, "data.frame") ) {
    data <- data.matrix(data)
  }

  # internal: generate stability for a specific `k`
  .calc_k <- function(.data, .i) {
    calc_stability(progeny_k(.data, k = .i, ...))
  }

  scores <- lapply(1:reps, function(.r) {
       if ( verbose ) {
         signal_done("Performing iteration repeat ...", value(.r))
       }
       vapply(clust_iter, .calc_k, .data = data, FUN.VALUE = 0.1)
    }) |> do.call(what = rbind)

  rdata <- scramble_data(data)   # scrambled runif data
  random_scores <- lapply(1:reps, function(.r) {
      vapply(clust_iter, .calc_k, .data = rdata, FUN.VALUE = 0.1)
    }) |>
    do.call(what = rbind)

  stopifnot(identical(dim(scores), dim(random_scores)))   # check
  colnames(scores)        <- sprintf("k=%i", clust_iter)
  colnames(random_scores) <- colnames(scores)

  cm <- be_hard(colMeans, na.rm = TRUE)  # set default for use below

  ret <- list()
  ret$scores      <- scores
  ret$mean_scores <- cm(scores)
  ret$ci95_scores <- apply(scores, 2, function(.x) {
    stats::quantile(.x, probs = c(0.025, 0.975), na.rm = TRUE)
  })
  ret$random_scores      <- random_scores
  ret$mean_random_scores <- cm(random_scores)
  ret$D_max       <- cm(scores) - cm(random_scores)
  ret$D_gap       <- calc_gap_score(scores)
  ret$clust_iter  <- clust_iter
  ret$repeats     <- reps
  ret             <- c(ret, list(...))
  ret$call        <- match.call(expand.dots = TRUE)
  add_class(ret, "pclust")
}


#' Test for object type "pclust"
#'
#' The [is_pclust()] function checks whether
#' an object is of class `pclust`. See [inherits()].
#'
#' @rdname progeny_cluster
#' @return [is_pclust()] returns a logical boolean.
#'
#' @examples
#' # Test for class `pclust`
#' is_pclust(pclust)
#'
#' @export
is_pclust <- function(x) inherits(x, "pclust")


#' Perform Progeny Clustering 1x (internal)
#'
#' Perform progeny clustering for 1 iteration of K in k-means clustering.
#'
#' @param data The data matrix (n x p) containing n samples and p features.
#' @param k Integer. The number of clusters.
#' @param iter Integer. The number of progeny sampling iterations to perform.
#' @param size Integer. The number of progeny to sample in each cluster.
#'   Must be less than the number of samples in the data matrix,
#'   typically > 10, minimum = 5.
#' @return A list containing:
#' \item{Pij:}{A `k*size` x `k*size` matrix of the individual probabilities of
#'   repeat clustering for each progeny sample. See reference for full description}
#' \item{size:}{The progeny sample size}
#' \item{k:}{The k integer used in k-means}
#' @author Stu Field
#' @importFrom stats kmeans
#' @noRd
progeny_k <- function(data, k, size, iter) {

  data <- data.matrix(data, rownames.force = FALSE)

  if ( size > nrow(data) * 0.9 ) {
    stop(
      "You are probably progeny sampling with too many samples for ",
      "these data ... try a smaller number.", call. = FALSE
    )
  }

  kn <- k * size

  if ( kn > nrow(data) / 2 ) {
    stop(
      "You are probably progeny sampling with too many samples ... ",
      "perhaps try `size` < ", value(size), ".", call. = FALSE
    )
  }

  orig_clust   <- stats::kmeans(data, centers = k, nstart = 25, iter.max = 20)
  cluster_list <- lapply(1:k, function(.k) which(orig_clust$cluster == .k))
  Qij          <- matrix(0, ncol = kn, nrow = kn)

  # loop replaces Qij matrix in place via `+` at each iter
  for ( r in 1:iter ) {
    progeny <- lapply(data.frame(data), function(.p) {
                 # cat("* feature", .p, "\n")  # nolint: commented_code_linter.
                 lapply(cluster_list, function(.clust) {
                   # sample entries of .p according to .clust
                   .p[sample(.clust, size, replace = TRUE)]
                 }) |> unlist()
    }) |> data.frame() |> data.matrix()

    new_clust <- stats::kmeans(progeny, k, nstart = 25, iter.max = 20)$cluster
    q_ij      <- outer(new_clust, new_clust, function(x, y) as.numeric(x == y))
    q_ij[ lower.tri(q_ij) ] <- 0
    q_ij <- q_ij + t(q_ij)
    stopifnot(isSymmetric(q_ij))
    Qij <- Qij + q_ij
  }
  diag(Qij) <- iter
  list(Pij = Qij / iter, size = size, k = k)
}

#' Clustering Stability Metric (internal)
#'
#' This is an internal function calculating `S`, the progeny stability metric.
#'
#' @param x An object created via `progenyK()`, primarily of a `Pij`
#'   probability matrix.
#' @return A progeny stability metric, the ratio of the true
#'   classifications vs. the false classifications.
#'
#' @note `x` is a Pij symmetric matrix.
#' @noRd
calc_stability <- function(x) {
  pij  <- x$Pij
  size <- x$size
  k    <- x$k

  trueprob <- vapply(1:k, function(.c) {
                     from <- to <- ((.c - 1) * size + 1):(.c * size)
                     sum(pij[from, to])}, double(1)) |> sum()

  falseprob   <- sum(pij) - trueprob
  true_score  <- ((trueprob - size * k) / ((size - 1) * size * k))
  false_score <- (falseprob / (size * (k - 1) * size * k))
  return(true_score / ifelse(false_score == 0, NA, false_score))
}

#' Calculate the Gap Score
#'
#' @param x A scores matrix.
#' @noRd
calc_gap_score <- function(x) {
  gmat <- vapply(2:(ncol(x) - 1), FUN.VALUE = double(nrow(x)),
                 function(.i) 2 * x[, .i] - x[, .i - 1] - x[, .i + 1])

  g1 <- x[, 1] - 2 * x[, 2] + x[, 3]
  g2 <- x[, ncol(x)] - 2 * x[, ncol(x) - 1] + x[, ncol(x) - 2]
  g  <- colMeans(cbind(g1, gmat, g2), na.rm = TRUE)
  names(g) <- colnames(x)
  g
}

#' Generate Reference Dataset of random noise from original data
#'   Do we want runif()? What about rnorm()?
#'
#' @param data A data matrix.
#' @return A data matrix.
#'
#' @importFrom stats runif
#' @noRd
scramble_data <- function(data) {
  n <- nrow(data)
  vapply(data.frame(data), function(.p) {
    stats::runif(n, min = min(.p, na.rm = TRUE),
                 max = max(.p, na.rm = TRUE))
    }, double(n))
}
