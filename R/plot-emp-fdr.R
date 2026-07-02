#' Plot Empirical FDR
#'
#' @describeIn calc_emp_fdr
#'   plots the mean number of false positives (FPs)
#'   versus the numberof selected features by a sequence of
#'   selection probability thresholds. For this to be possible,
#'   the `stab_sel` object *must* have permuted data in order to
#'   calculate empirical false discovery rates. The area in the
#'   sub-diagonal represents where more features are added without
#'   a commensurate increase in false positives (**Good**).
#'   The inverse is true for the super-diagonal, false positives
#'   are being included faster than additional features (**Bad**).
#'   The legend highlights pre-defined empirical FDR breaks:
#'   `c(0.5, 1, 2, 3, 5)` evaluated to the nearest threshold cutoff.
#'
#' @examples
#' # plot the FDR
#' plot_emp_fdr(ss)  # typically set permutations > 75
#'
#' @importFrom ggplot2 ggplot aes geom_abline geom_point theme
#' @importFrom ggplot2 scale_color_manual labs annotate geom_step
#' @export
plot_emp_fdr <- function(x, thresh_seq = seq(1, 0.1, by = -0.01)) {

  stopifnot(is_stab_sel(x), x$perm_data)

  L <- length(x$permpath_list)
  emp_breaks <- calc_emp_fdr_breaks(x, thresh_seq = thresh_seq)
  emp_breaks$breaks$thresh_mean <- paste0(emp_breaks$breaks$piThresh, " | ",
                                          emp_breaks$breaks$MeanFPs)

  emp_breaks$fdr_data |>
    ggplot(aes(x = n_selected, y = MeanFPs)) +
    geom_step(size = 0.5) +
    geom_point(alpha = 0.75, size = 3.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(data = emp_breaks$breaks, size = 2.5,
               aes(x = n_selected, y = MeanFPs, colour = thresh_mean)) +
    scale_color_manual(name   = "threshold | Mean FPs",
                       values = head(unname(col_palette),
                                     nrow(emp_breaks$breaks))) +
    labs(
      x = "Number Selected Features",
      y = sprintf("Mean FPs over %i permuted (null) datasets", L),
      title = "Mean FPs vs. Number Stable Features"
      ) +
    annotate("text", label = "Bad",
             x = max(emp_breaks$fdr_data$n_selected) * 0.1,
             y = max(emp_breaks$fdr_data$MeanFPs) * 2 / 3, size = 10) +
    annotate("text", label = "Good",
             x = max(emp_breaks$fdr_data$n_selected) * 2 / 3,
             y = max(emp_breaks$fdr_data$MeanFPs) * 0.1, size = 10) +
    theme(legend.position = "right") +
    NULL
}
