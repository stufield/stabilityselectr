#' Test for object type "stab_sel"
#'
#' The [is_stab_sel()] function checks whether
#'   an object is class `stab_sel`. See [inherits()].
#'
#' @rdname stability_selection
#'
#' @return The `is_stab_sel` function returns a logical boolean.
#'
#' @examples
#' # Test for class `stab_sel`
#' is_stab_sel(stab_sel)
#'
#' @export
is_stab_sel <- function(x) inherits(x, "stab_sel")


#' @describeIn stability_selection
#'   S3 `print` method for class `stab_sel`.
#'
#' @return The S3 print method returns:
#' \item{Stability Selection Kernel}{The kernel used in
#'       the stability selection algorithm.}
#'
#' \item{weakness}{the weakness used (\eqn{\alpha}).}
#'
#' \item{P(alpha)}{the probability of the \eqn{\alpha}
#'       being applied (`Pw =` argument).}
#'
#' \item{iter}{number of iterations in the stability selection algorithm.}
#'
#' \item{standardized}{was the data pre-standardized?}
#'
#' \item{imputed outliers}{were statistical outliers imputed via
#'       Gaussian approximation prior to stability selection?}
#'
#' \item{max lambda}{the maximum \eqn{\lambda} (tuning parameter)
#'       used.}
#'
#' \item{min ratio}{the minimum ratio used: \eqn{\lambda/max(\lambda)}.}
#'
#' \item{permuted data}{was permuted (null) data generated during
#'       stability selection?.}
#'
#' \item{random seed}{the seed passed to the random number generator
#'       during stability selection.}
#'
#' @examples
#' # S3 print method
#' stab_sel
#'
#' @export
print.stab_sel <- function(x, ...) {
  signal_rule(
    sprintf("Stability Selection (%s: %s)",
            add_color("Kernel", "cyan"),
            add_color(x$kernel, "green")),
    line_col = "blue", lty = "double"
  )
  key <- c(
    "weakness (\u03B1)",
    "P(\u03B1) -> (Pw)",
    "iter",
    "standardized",
    "imputed outliers",
    "max \u03BB",
    "min ratio (\u03BB)",
    "permuted data",
    "random Seed"
    ) |> pad(27)
  value <- list(
    round(x$alpha, 2L),
    round(x$Pw, 2L),
    round(x$n_iter),
    x$standardize,
    x$impute_outliers,
    round(x$lambda[1L], 4L),
    round(x$lambda_min_ratio, 1L),
    x$perm_data,
    x$r_seed
  )
  liter(key, value, function(.x, .y) {
    writeLines(paste(add_color(symbl$bullet, "red"), .x, value(.y)))
  })
  signal_rule(line_col = "green", lty = "double")
  invisible(x)
}


#' @describeIn stability_selection
#'   The S3 `summary` method for class `stab_sel`.
#'
#' @inheritParams get_stable_features
#'
#' @param object A `stab_sel` class object.
#'
#' @note Additional features can be passed as strings to the
#'   summary method via the `add_features` argument.
#'
#' @examples
#' # S3 summary method
#' summary(stab_sel, thresh = 0.6)
#'
#' summary(stab_sel, thresh = 0.8, add_features = "feat_c")   # force feat_c into table
#'
#' @importFrom dplyr select everything matches
#'
#' @export
summary.stab_sel <- function(object, thresh = NULL, ...) {
  if ( is.null(thresh) ) {
    stop(
      "Must pass a stability threshold to S3 summary method. ",
      "Please pass `thresh =`.", call. = FALSE
    )
  }

  lambda_norm <- object$lambda / max(object$lambda)
  df <- get_stable_features(object, thresh = thresh, ...)[[1L]]

  if ( nrow(df) > 0L ) { # reorder feats
    df$AUC <- object$stabpath_matrix[df$feature, , drop = FALSE] |>
      calc_path_auc(values = lambda_norm)
  }

  dplyr::select(df, -matches("FDR"), everything())
}


#' @describeIn stability_selection
#'   The S3 `plot` method plots the selection paths for the features.
#'   This plot closely resembles a lasso coefficient plot with
#'   the regularization parameter (lambda) plotted on x-axis and
#'   the feature selection probability (rather than the model
#'   coefficient) is plotted on the y-axis.
#'
#' Plots the regularization parameter (lambda) on the x-axis and the
#'   selection probability on the y-axis. The regularization
#'   parameter is plotted as lambda/max(lambda) so that it is in
#'   the range from 1 to 0. The selection probability corresponds
#'   to the number of times a particular marker was chosen at a
#'   given value of lambda. Each line in the plot is a marker and
#'   represents the stability selection path over the range of
#'   regularization parameter. All features that have a maximum
#'   selection probability greater than `thresh` (shown as a dotted
#'   horizontal line) are colored and labeled and the remaining
#'   features are colored gray and unlabeled. Additionally, you can
#'   provide a set of custom labels that will be colored and labeled
#'   regardless of their max selection probability. Each feature is
#'   labeled with a capital letter and the full name of the feature
#'   is indicated in the legend along with the AUC for its curve in
#'   parentheses.
#'
#' @param custom_labels `character(n)`. Character vector of
#'   additional features to label in the plot, see `Details`.
#'
#' @param main `character(1)`. Optional plot title (default depends on kernel).
#'
#' @param sort_by_AUC `logical(1)`. If `TRUE`, entries in
#'   the legend will be sorted by their curve AUC values
#'   which are in parentheses following the variable name
#'   in the legend.
#'
#' @param ln_cols `character(n)`. A vector of colors to be used as
#'   line colors in plotting. Recycled as necessary.
#'
#' @param add_perm `logical(1)`. Should empirical false discovery lines
#'   from the null permutation be added to the plot
#'  (only if permutation was performed)? This can be time
#'  consuming depending on the number of permutations.
#'
#' @param emp_thresh `numeric(n)`. A vector describing the empirical
#'   threshold values to be used (`default = seq(1, 0.1, by = 0.01)`).
#'
#' @return A `ggplot`.
#'
#' @examples
#' # S3 plot method
#' plot(stab_sel, thresh = 0.8)
#'
#' @importFrom utils head
#' @importFrom dplyr case_when arrange select mutate desc left_join
#' @importFrom ggplot2 theme geom_line scale_colour_manual aes
#' @importFrom ggplot2 scale_x_reverse geom_hline guides guide_legend
#' @importFrom ggplot2 scale_linetype_manual scale_y_continuous
#' @export
plot.stab_sel <- function(x, thresh = 0.60,
                          custom_labels = NULL,
                          main = NULL,
                          sort_by_AUC = TRUE,
                          ln_cols = col_palette,
                          add_perm = FALSE,
                          emp_thresh = seq(1, 0.1, by = -0.01), ...) {

  lambda_norm <- x$lambda / max(x$lambda)

  if ( is.null(main) ) {
    main <- bquote(
      "Stability Paths (" * alpha == .(x$alpha) *
        "," ~ italic(P)(alpha) == .(x$Pw) * ")")
  }

  custom_features <- intersect(custom_labels, rownames(x$stabpath_matrix))
  summary_df      <- summary(x, thresh = thresh, warn = TRUE,
                             add_features = custom_features)
  L <- nrow(summary_df)

  if ( L == 0L ) {
    summary_df$AUC <- NA_real_
  }

  if ( sort_by_AUC ) {
    summary_df <- dplyr::arrange(summary_df, dplyr::desc(AUC))
  }

  legend_cols <- dplyr::case_when(L > 12 & L < 24 ~ 2,
                                  L > 24 & L < 36 ~ 3,
                                  L > 36 ~ 4,
                                  TRUE ~ 1)

  # pre-process the x-axis values for dplyr::left_join to the prob. paths
  # The names will be X1 ... Xn; following conversion to tbl_df
  x_df <- data.frame(
    column      = paste0("X", seq_len(ncol(x$stabpath_matrix))),
    lambda_norm = lambda_norm)

  auc_df <- dplyr::select(summary_df, feature, AUC)

  # internal helper
  .strip_seq <- function(x) sub("\\.", "-", sub("^seq\\.", "", x))

  # get order of feature labels
  label_order <- sprintf("%s (%0.2f)", .strip_seq(auc_df$feature), auc_df$AUC)

  line_cols <- rep(unname(ln_cols), length.out = L)
  names(line_cols) <- label_order

  if ( L < nrow(x$stabpath_matrix) ) {  # if un-selected features, append color map
    line_cols <- c(line_cols, "Unselected (NA)" = "grey80")
  }

  data <- data.frame(x$stabpath_matrix) |>
    rn2col("feature") |>
    left_join(auc_df, by = "feature") |>
    dplyr::mutate(
      feat_sel = case_when(
        feature %in% summary_df$feature ~ .strip_seq(feature),
        TRUE ~ "Unselected"),
      label = factor(sprintf("%s (%0.2f)", feat_sel, AUC),
                     levels = names(line_cols))
    ) |>
    tidyr::pivot_longer(c(-feature, -feat_sel, -AUC, -label),
                        names_to = "column", values_to = "prob") |>
    left_join(x_df, by = "column")

  # sanity check for safety
  stopifnot(all(names(line_cols) %in% unique(data$label)))

  p <- data |>
    ggplot(aes(x = lambda_norm, y = prob)) +
    geom_line(aes(group = feature, colour = label)) +
    scale_colour_manual(values = line_cols) +
    scale_x_reverse(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = thresh, alpha = 0.75,
               colour = "black", linetype = "dashed") +
    guides(colour = guide_legend(title = "Feature", ncol = legend_cols)) +
    labs(x = bquote(lambda ~ "/" ~ max()(lambda)),
         y = "Selection Probability",
         title = main) +
    theme(legend.position = "right")

  if ( add_perm ) {
    if ( !x$perm_data ) {
      warning(
        "No permutation data present in `stab_sel` object. ",
        "Skipping permutation cutoffs.", call. = FALSE
      )
    } else {
      emp_breaks_list <- calc_emp_fdr_breaks(x, emp_thresh)
      emp_breaks_list$breaks$line <- "dotted"
      break_cols <- unname(col_palette) |> head(5L)
      p <- p +
        geom_hline(
          data = emp_breaks_list$breaks,
          mapping = aes(yintercept = piThresh, linetype = factor(MeanFPs)),
          colour = break_cols
          ) +
        scale_linetype_manual(
          name = "MeanFPs",
          values = rep("dashed", length(break_cols))
        ) +
        guides(linetype = guide_legend(
          override.aes = list(colour = break_cols))
        )
    }
  }
  p
}
