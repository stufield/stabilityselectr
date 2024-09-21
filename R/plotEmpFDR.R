#' Plot Empirical FDR
#'
#' Plot the mean number of false positives (FPs) versus the number of selected
#'   features by a sequence of selection probability thresholds. For this to
#'   be possible, the `stab_sel` object _must_ have permuted
#'   data in order to calculate empirical false discovery rates. The area
#'   in the sub-diagonal represents where more features are added without
#'   a commensurate increase in false positives (**Good**). The inverse is true
#'   for the super-diagonal, false positives are being included faster
#'   than additional features (**Bad**). The legend highlights pre-defined
#'   empirical FDR breaks: `c(0.5, 1, 2, 3, 5)` evaluated to the nearest
#'   threshold cutoff.
#'
#' @family empirical FDR
#' @inheritParams calcEmpFDR
#' @author Stu Field
#' @seealso [stabilitySelection()], [getStableFeatures()]
#' @examples
#' # l1-logistic
#' withr::with_seed(101, {
#'   n_feat      <- 20
#'   n_samples   <- 100
#'   x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
#'   colnames(x) <- paste0("feat", "_", head(letters, n_feat))
#'   y           <- sample(1:2, n_samples, replace = TRUE)
#' })
#'
#' # typically set > 75 permutations
#' stab_sel <- stabilitySelection(x, y, "l1-logistic", num.perms = 25,
#'                                r.seed = 101, parallel = TRUE)
#' plotEmpFDR(stab_sel)
#'
#' @importFrom ggplot2 ggplot aes geom_abline geom_point
#' @importFrom ggplot2 scale_color_manual labs annotate
#' @importFrom SomaPlotr theme_soma
#' @export
plotEmpFDR <- function(x, thresh.seq = seq(1, 0.1, by = -0.01)) {

  stopifnot(is.stab_sel(x), x$perm.data)

  L <- length(x$permpath.list)
  emp_breaks <- calcEmpFDRbreaks(x, thresh.seq = thresh.seq)
  emp_breaks$breaks$thresh_mean <- paste0(emp_breaks$breaks$piThresh, " | ",
                                          emp_breaks$breaks$MeanFPs)

  emp_breaks$fdr_data |>
    ggplot2::ggplot(ggplot2::aes(x = n_selected, y = MeanFPs)) +
    ggplot2::geom_step(size = 0.5) +
    ggplot2::geom_point(alpha = 0.75, size = 3.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_point(data = emp_breaks$breaks, size = 2.5,
                        ggplot2::aes(x = n_selected,
                                     y = MeanFPs, colour = thresh_mean)) +
    ggplot2::scale_color_manual(name = "threshold | Mean FPs",
                                values = head(unlist(SomaPlotr::soma_colors),
                                              nrow(emp_breaks$breaks))) +
    ggplot2::labs(
      x = "Number Selected Features",
      y = sprintf("Mean FPs over %i permuted (null) datasets", L),
      title = "Mean FPs vs. Number Stable Features"
      ) +
    ggplot2::annotate("text", label = "Bad",
                      x = max(emp_breaks$fdr_data$n_selected) * 0.1,
                      y = max(emp_breaks$fdr_data$MeanFPs) * 2 / 3, size = 10) +
    ggplot2::annotate("text", label = "Good",
                      x = max(emp_breaks$fdr_data$n_selected) * 2 / 3,
                      y = max(emp_breaks$fdr_data$MeanFPs) * 0.1, size = 10) +
    theme_soma(legend.position = "right") +
    NULL
}
