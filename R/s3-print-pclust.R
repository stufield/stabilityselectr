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

  signal_rule("Progeny Clustering", line_col = "blue", lty = "double")
  key <- c(
    "Call",
    "Progeny Size",
    "K iterations",
    "No. iterations",
    "No. repeats"
  ) |> pad(17)
  value <- list(
    gsub(" {2,}", " ", paste(deparse(x$call), collapse = "")),
    x$size,
    paste(min(x$clust_iter), symbl$arrow_right, max(x$clust_iter)),
    x$n_iter,
    x$repeats
  )
  liter(key, value, function(.x, .y) {
    writeLines(paste(add_color(symbl$bullet, "red"), .x, value(.y)))
  })

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
