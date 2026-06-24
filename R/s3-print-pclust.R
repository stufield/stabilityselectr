#' S3 Print Method for `pclust`
#'
#' @rdname progeny_cluster
#'
#' @importFrom withr local_options
#'
#' @export
print.pclust <- function(x, ...) {

  local_options(c(digits = 3L))

  starMax <- function(v) {
    smb <- symbl$star
    names(v)[which.max(v)] <- paste0(smb, names(v)[which.max(v)], smb)
    v
  }

  signal_rule("Progeny Cluster Object", line_col = "blue", lty = "double")
  key <- c(
    "Call",
    "Progeny Size",
    "No. iterations",
    "K iterations",
    "No. repeats"
  ) |> pad(25)
  value <- c(
    gsub(" {2,}", " ", paste(deparse(x$call), collapse = "")),
    x$size,
    x$n_iter,
    paste(min(x$clust_iter), symbl$arrow_right, max(x$clust_iter)),
    x$repeats
  )
  writeLines(paste0("   ", key, value))

  cat("\n")
  signal_rule("Mean & CI95 Stability Scores", line_col = "cyan")
  print(rbind(mean = starMax(x$mean_scores),
              x$ci95_scores)[c(2L, 1L, 3L), ])

  cat("\n")
  signal_rule("Maximum Distance Scores", line_col = "cyan")
  print(starMax(x$D_max))

  cat("\n")
  signal_rule("Gap Distance Scores", line_col = "cyan")
  print(starMax(x$D_gap))

  signal_rule(line_col = "green", lty = "double")
  invisible(x)
}
