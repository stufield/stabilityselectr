# Calculate Empirical FDR Break Points

Calculates the stability selection threshold, the mean number of false
positive selected features (empirical), and the number of selected
features for specified FDR break points. Relies on
[`calc_emp_fdr()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr.md)
to calculate the mean false discovery based on permutations during the
stability selection algorithm.

## Usage

``` r
calc_emp_fdr_breaks(
  x,
  thresh_seq = seq(1, 0.1, by = -0.01),
  fdr_breaks = c(0.5, 1, 2, 3, 5)
)
```

## Arguments

- x:

  A `stab_sel` class object generated via
  [`stability_selection()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_selection.md).

- thresh_seq:

  `numeric(n)`. A sequence in `[0, 1]` specifying the thresholds to
  evaluate.

- fdr_breaks:

  `numeric(n)`. A vector specifying the desired mean number of empirical
  false positives at which to calculate various thresholds.

## Value

A list consisting of:

- n_selected:

  A vector of the number of features selected at each empirical
  stability selection threshold

- meanFPs:

  A vector of the mean number of false positive selected features at
  each empirical stability selection threshold

- breaks:

  A `tibble` of containing empirical false positive summary statistics
  at each FDR specified break point

## See also

[`get_stable_features()`](https://stufield.github.io/stabilityselectr/dev/reference/get_stable_features.md),
[`stability_selection()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_selection.md)

Other empirical FDR:
[`calc_emp_fdr()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr.md),
[`plot_emp_fdr()`](https://stufield.github.io/stabilityselectr/dev/reference/plot_emp_fdr.md)

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
  y        <- sample(1:2, n_samples, replace = TRUE)
})
stab_sel <- stability_selection(x, y, "l1-logistic", num_iter = 25,
                                num_perms = 25,
                                r_seed = 101, parallel = TRUE)
#> ✓ Using kernel: 'l1-logistic' and 1 core (serial)
calc_emp_fdr_breaks(stab_sel)
#> $fdr_data
#> # A tibble: 91 × 3
#>    MeanFPs n_selected piThresh
#>      <dbl>      <int>    <dbl>
#>  1    0.32          0     1   
#>  2    0.32          0     0.99
#>  3    0.8           0     0.98
#>  4    0.8           0     0.97
#>  5    1.6           0     0.96
#>  6    1.6           0     0.95
#>  7    2.24          2     0.94
#>  8    2.24          2     0.93
#>  9    3.68          3     0.92
#> 10    3.68          3     0.91
#> # ℹ 81 more rows
#> 
#> $breaks
#> # A tibble: 5 × 4
#>   FDR_breaks MeanFPs n_selected piThresh
#>        <dbl>   <dbl>      <int>    <dbl>
#> 1        0.5    0.8           0     0.98
#> 2        1      1.6           0     0.96
#> 3        2      2.24          2     0.94
#> 4        3      3.68          3     0.92
#> 5        5      6.48         11     0.88
#> 
```
