# Calculate Stable Features

Returns a data frame object of all features with a maximum selection
probability greater than a minimum threshold.

`get_threshold_features()` calculates the features at or above a given
selection probability threshold. It is a thin wrapper around
`get_stable_features()` returning a list of thresholds corresponding to
`thresh_vec`.

## Usage

``` r
get_stable_features(x, thresh, add_features, warn)

get_threshold_features(x, thresh_vec, ...)
```

## Arguments

- x:

  An `stab_sel` class object OR a matrix containing selection
  probabilities, i.e. the `stabpath_matrix` entry of a `stab_sel` class
  object.

- thresh:

  `numeric(1)` in `[0, 1]`. Minimum selection probability threshold.

- add_features:

  `character(n)`. A string of additional features to *force* into the
  resulting table, irrespective of their threshold. Used mostly in the
  S3 plot method to see a given stability path of a feature not meeting
  a threshold cutoff. Must be exact string match.

- warn:

  `logical(1)`. Should warnings be triggered if no stable features were
  found at the specified threshold OR if the FDR upper bound is
  undefined at thresholds `<= 0.5`?

- thresh_vec:

  `numeric(1)`. A vector of threshold values.

- ...:

  Additional arguments passed to `get_stable_features()`, typically the
  `add_features =` and `warn =` arguments.

## Value

A two column data frame containing maximum selection probabilities and
FDR upper bounds when appropriate, see `Details`.

For `get_threshold_features()`, a *list* of data frames for various
thresholds corresponding to the entries of `thresh_vec`.

## Details

A "stable feature" is defined as a feature with a maximum selection
probability greater than a supplied threshold. This function returns a
`data.frame` of all features that satisfy this criterion along with the
maximum selection probability and the upper bound on the false discover
rate. This false discovery rate bound is only defined for
`thresh > 0.5`, it is otherwise undefined.

*IMPORTANT!* If you pass to the `matrix` method, there is no permutation
analysis performed, i.e. no `$EmpFDR` column in the returned data frame.
This calculation takes a long time and is not always desired, so this
method offers the user a control mechanism for the output behavior.

## See also

[`stability_selection()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_selection.md)

## Author

Stu Field

## Examples

``` r
# l1-logistic
withr::with_seed(101, {
  n_feat      <- 20
  n_samples   <- 100
  x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
  colnames(x) <- paste0("feat", "_", head(letters, n_feat))
  y           <- sample(1:2, n_samples, replace = TRUE)
})

stab_sel <- stability_selection(x, y)
#> ✓ Using kernel: 'l1-logistic' and 1 core (serial)

# Stable features at thresh = 0.55
get_stable_features(stab_sel, 0.55)
#>        MaxSelectProb FDRbound
#> feat_d         0.945    0.025
#> feat_j         0.945    0.050
#> feat_m         0.945    0.075
#> feat_s         0.935    0.100
#> feat_t         0.935    0.125
#> feat_c         0.905    0.150
#> feat_a         0.900    0.175
#> feat_r         0.900    0.200
#> feat_f         0.885    0.225
#> feat_e         0.880    0.250
#> feat_g         0.875    0.275
#> feat_h         0.875    0.300
#> feat_n         0.865    0.325
#> feat_o         0.860    0.350
#> feat_q         0.855    0.375
#> feat_l         0.850    0.400
#> feat_p         0.850    0.425
#> feat_b         0.835    0.450
#> feat_i         0.800    0.475
#> feat_k         0.800    0.500
# l1-logistic
withr::with_seed(101, {
  n_feat      <- 20
  n_samples   <- 100
  x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
  colnames(x) <- paste0("feat", "_", head(letters, n_feat))
  y           <- sample(1:2, n_samples, replace = TRUE)
})
stab_sel <- stability_selection(x, y, "l1-logistic", r_seed = 101)
#> ✓ Using kernel: 'l1-logistic' and 1 core (serial)
get_threshold_features(stab_sel, seq(0.6, 0.9, 0.1))
#> $thresh_0.6
#>        MaxSelectProb FDRbound
#> feat_d         0.945   0.0125
#> feat_t         0.935   0.0250
#> feat_a         0.915   0.0375
#> feat_j         0.915   0.0500
#> feat_s         0.905   0.0625
#> feat_m         0.900   0.0750
#> feat_l         0.890   0.0875
#> feat_f         0.885   0.1000
#> feat_q         0.880   0.1125
#> feat_g         0.870   0.1250
#> feat_e         0.860   0.1375
#> feat_n         0.860   0.1500
#> feat_r         0.860   0.1625
#> feat_c         0.855   0.1750
#> feat_k         0.850   0.1875
#> feat_b         0.845   0.2000
#> feat_o         0.840   0.2125
#> feat_h         0.835   0.2250
#> feat_p         0.830   0.2375
#> feat_i         0.785   0.2500
#> 
#> $thresh_0.7
#>        MaxSelectProb FDRbound
#> feat_d         0.945  0.00625
#> feat_t         0.935  0.01250
#> feat_a         0.915  0.01875
#> feat_j         0.915  0.02500
#> feat_s         0.905  0.03125
#> feat_m         0.900  0.03750
#> feat_l         0.890  0.04375
#> feat_f         0.885  0.05000
#> feat_q         0.880  0.05625
#> feat_g         0.870  0.06250
#> feat_e         0.860  0.06875
#> feat_n         0.860  0.07500
#> feat_r         0.860  0.08125
#> feat_c         0.855  0.08750
#> feat_k         0.850  0.09375
#> feat_b         0.845  0.10000
#> feat_o         0.840  0.10625
#> feat_h         0.835  0.11250
#> feat_p         0.830  0.11875
#> feat_i         0.785  0.12500
#> 
#> $thresh_0.8
#>        MaxSelectProb    FDRbound
#> feat_d         0.945 0.004166667
#> feat_t         0.935 0.008333333
#> feat_a         0.915 0.012500000
#> feat_j         0.915 0.016666667
#> feat_s         0.905 0.020833333
#> feat_m         0.900 0.025000000
#> feat_l         0.890 0.029166667
#> feat_f         0.885 0.033333333
#> feat_q         0.880 0.037500000
#> feat_g         0.870 0.041666667
#> feat_e         0.860 0.045833333
#> feat_n         0.860 0.050000000
#> feat_r         0.860 0.054166667
#> feat_c         0.855 0.058333333
#> feat_k         0.850 0.062500000
#> feat_b         0.845 0.066666667
#> feat_o         0.840 0.070833333
#> feat_h         0.835 0.075000000
#> feat_p         0.830 0.079166667
#> 
#> $thresh_0.9
#>        MaxSelectProb FDRbound
#> feat_d         0.945 0.003125
#> feat_t         0.935 0.006250
#> feat_a         0.915 0.009375
#> feat_j         0.915 0.012500
#> feat_s         0.905 0.015625
#> feat_m         0.900 0.018750
#> 
```
