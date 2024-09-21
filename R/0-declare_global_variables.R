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
      "FdrBreaks",     # calcEmpFDRbreaks in `dplyr::mutate()`
      "MaxSelectProb", # in `getStableFeatures()`
      "AUC",           # in `plot.stab_sel()`
      "feature",       # in `plot.stab_sel()`
      "feat_sel",      # in `plot.stab_sel()`
      "prob",          # in `plot.stab_sel()`
      "label",         # in `plot.stab_sel()`
      "thresh_mean",   # in `plotEmpFDR()`
      "MeanFPs",       # in `plotEmpFDR()`
      "piThresh",      # in `plotEmpFDR()`
      "n_selected"     # in `plotEmpFDR()`
    )
  )
}
