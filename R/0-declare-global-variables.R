# --------------------------- #
# Declaring Global Variables:
# This is mostly for passing R CMD checks
# global variables that come from other dependant
# packages, or objects in the data/ directory
# Reference: https://github.com/tidyverse/magrittr/issues/29
# -------------------------------------------------------- #
if ( getRversion() >= "2.15.1" ) {
  utils::globalVariables(
    c(".",
      "FDR_breaks",    # calc_emp_fdr_breaks() in `dplyr::mutate()`
      "null_score", "score", "progeny_score",         # `plot.pclust()`
      "2.5%", "97.5%", "cluster", "type", "gap_dist", # `plot.pclust()`
      "max_dist", "score type",        # `plot.pclust()`
      "MaxSelectProb", # in `get_stable_features()`
      "AUC",           # in `plot.stab_sel()`
      "feature",       # in `plot.stab_sel()`
      "feat_sel",      # in `plot.stab_sel()`
      "prob",          # in `plot.stab_sel()`
      "label",         # in `plot.stab_sel()`
      "thresh_mean",   # in `plot_emp_fdr()`
      "MeanFPs",       # in `plot_emp_fdr()`
      "piThresh",      # in `plot_emp_fdr()`
      "n_selected"     # in `plot_emp_fdr()`
    )
  )
}
