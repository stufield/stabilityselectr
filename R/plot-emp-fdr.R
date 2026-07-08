#' Plot Empirical FDR
#'
#' @describeIn calc_emp_fdr
#'   plots the mean number of false positives (FPs)
#'   versus the number of selected features by a sequence of
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
#' @param ... Additional arguments passed to either
#'   [calc_emp_fdr_breaks()] (i.e. `thresh_seq` and `fdr_breaks`)
#'   or [plot.stab_sel()].
#'
#' @examples
#' # plot the FDR
#' plot_emp_fdr(ss)  # typically set permutations > 75
#'
#' @importFrom ggplot2 ggplot aes geom_abline geom_point theme
#' @importFrom ggplot2 scale_color_manual labs annotate geom_step
#' @importFrom ggplot2 scale_y_continuous scale_x_continuous
#' @export
plot_emp_fdr <- function(x, ...) {

  .check_perm(x)

  L <- length(x$permpath_list)
  emp_breaks <- calc_emp_fdr_breaks(x, ...)
  emp_breaks$breaks$thresh_mean <- paste(emp_breaks$breaks$piThresh, "|",
                                         emp_breaks$breaks$MeanFPs)

  emp_breaks$fdr_data |>
    ggplot(aes(x = n_selected, y = MeanFPs)) +
    geom_step(linewidth = 0.5) +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_continuous(breaks = scales::breaks_pretty()) +
    geom_point(alpha = 0.5, size = 3.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(data = emp_breaks$breaks, size = 2.5,
               aes(x = n_selected, y = MeanFPs, colour = thresh_mean)) +
    scale_color_manual(name   = "threshold | Mean FPs",
                       values = unname(col_palette)) +
    labs(
      x = "Number Selected Features",
      y = sprintf("Mean FPs over %i permuted (null) datasets", L),
      title = "False Positives vs. Number Stable Features"
      ) +
    annotate("text", label = "Bad", size = 10, alpha = 0.75,
             x = -Inf, y = Inf, hjust = -0.5, vjust = 2) +
    annotate("text", label = "Good", size = 10, alpha = 0.75,
             x = Inf,  y = -Inf, hjust = 1.25, vjust = -0.75) +
    theme(legend.position = "right") +
    NULL
}
