# Calculate Stable Features

Returns a data frame object of all features with a maximum selection
probability greater than a minimum threshold.

## Usage

``` r
get_stable_features(x, thresh, add_features, warn)
```

## Arguments

- x:

  An `stab_sel` class object OR a matrix containing selection
  probabilities, i.e. the `stabpath_matrix` entry of a `stab_sel`
  object.

- thresh:

  `numeric(1)` in `[0, 1]`. the minimum selection probability threshold.
  In some instances this value can also be a vector, but is generally a
  scalar \> 0.50 for `get_stable_features()`.

- add_features:

  `character(n)`. A string of additional features to *force* into the
  resulting table, irrespective of their threshold. Used mostly in the
  S3 plot method to see a given stability path of a feature not meeting
  a threshold cutoff. Must be exact string match.

- warn:

  `logical(1)`. Should warnings be triggered if no stable features were
  found at the specified threshold OR if the FDR upper bound is
  undefined at thresholds `<= 0.5`?

## Value

A two column data frame containing maximum selection probabilities and
FDR upper bounds when appropriate, see `Details`.

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
# logistic regression
n_feat      <- 20
n_samples   <- 100
x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
colnames(x) <- paste0("feat", "_", head(letters, n_feat))
y           <- sample(1:2, n_samples, replace = TRUE)

stab_sel <- stability_selection(x, y, r_seed = 101)
#> ✓ Using kernel: 'binomial' and 1 core (serial)
#> ✓ Stablity path run time: 0.36s

# Stable features at `thresh =`
get_stable_features(stab_sel, 0.75)
#> $thresh_0.75
#> # A tibble: 18 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.98     0.005
#>  2 feat_e          0.95     0.01 
#>  3 feat_j          0.895    0.015
#>  4 feat_p          0.87     0.02 
#>  5 feat_h          0.865    0.025
#>  6 feat_i          0.85     0.03 
#>  7 feat_g          0.835    0.035
#>  8 feat_n          0.83     0.04 
#>  9 feat_s          0.83     0.045
#> 10 feat_a          0.825    0.05 
#> 11 feat_f          0.825    0.055
#> 12 feat_m          0.82     0.06 
#> 13 feat_t          0.815    0.065
#> 14 feat_o          0.8      0.07 
#> 15 feat_b          0.79     0.075
#> 16 feat_q          0.775    0.08 
#> 17 feat_r          0.77     0.085
#> 18 feat_l          0.76     0.09 
#> 

get_stable_features(stab_sel, c(0.75, 0.9))
#> $thresh_0.75
#> # A tibble: 18 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.98     0.005
#>  2 feat_e          0.95     0.01 
#>  3 feat_j          0.895    0.015
#>  4 feat_p          0.87     0.02 
#>  5 feat_h          0.865    0.025
#>  6 feat_i          0.85     0.03 
#>  7 feat_g          0.835    0.035
#>  8 feat_n          0.83     0.04 
#>  9 feat_s          0.83     0.045
#> 10 feat_a          0.825    0.05 
#> 11 feat_f          0.825    0.055
#> 12 feat_m          0.82     0.06 
#> 13 feat_t          0.815    0.065
#> 14 feat_o          0.8      0.07 
#> 15 feat_b          0.79     0.075
#> 16 feat_q          0.775    0.08 
#> 17 feat_r          0.77     0.085
#> 18 feat_l          0.76     0.09 
#> 
#> $thresh_0.9
#> # A tibble: 2 × 3
#>   feature MaxSelectProb FDRbound
#>   <chr>           <dbl>    <dbl>
#> 1 feat_d           0.98  0.00312
#> 2 feat_e           0.95  0.00625
#> 
```
