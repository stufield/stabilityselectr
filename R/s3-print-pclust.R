#' S3 Print Method for Class `pclust`
#'
#' @rdname progenyCluster
#' @param x An object of class `pclust`.
#' @examples
#' # S3 print method
#' pclust
#' clustIris
#'
#' @importFrom withr local_options
#' @export
print.pclust <- function(x, ...) {

  local_options(c(digits = 3L))

  starMax <- function(v) {
    names(v)[which.max(v)] <- paste0(names(v)[which.max(v)], "*")
    v
  }

  signal_rule("Progeny Cluster Object", line_col = "blue", lty = "double")
  key <- c(
    "Call",
    "Progeny Size",
    "No. of Iterations",
    "K Iterations"
  ) |> pad(25)
  value <- c(
    gsub(" {2,}", " ", paste(deparse(x$call), collapse = "")),
    x$size,
    x$iter,
    paste0(x$clust.iter, collapse = " ")
  )
  writeLines(paste0("   ", key, value))

  cat("\n")
  signal_rule(
    "Mean & CI95 Stability Scores", line_col = "magenta", lty = "double"
  )
  print(rbind(starMax(x$mean.scores), x$ci95.scores)[c(2, 1, 3), ])

  cat("\n")
  signal_rule("Maximum Distance Scores", line_col = "magenta", lty = "double")
  print(starMax(x$D.max))

  cat("\n")
  signal_rule("Gap Distance Scores", line_col = "magenta", lty = "double")
  print(starMax(x$D.gap))
  signal_rule(line_col = "green", lty = "double")
  invisible(x)
}
