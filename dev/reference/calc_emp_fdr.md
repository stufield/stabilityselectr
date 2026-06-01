# Calculate Empirical Number of False Positives

Calculate the mean number of false positive features from a permutation
analysis performed during a stability selection run. This assumes a
permutation set was generated during stability selection, (i.e.
`num_perms > 0`).

## Usage

``` r
calc_emp_fdr(x, thresh_seq, warn = TRUE)
```

## Arguments

- x:

  A `stab_sel` class object generated via
  [`stability_selection()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_selection.md).

- thresh_seq:

  `numeric(n)`. A sequence in `[0, 1]` specifying the thresholds to
  evaluate.

- warn:

  `logical(1)`. Should warnings be triggered if mean of `< 5`
  permutations is being returned?

## Value

A named vector indicating the average number (counts) of false positive
features selected at the various thresholds specified by `thresh_seq`.

## See also

[`get_threshold_features()`](https://stufield.github.io/stabilityselectr/dev/reference/get_stable_features.md)

Other empirical FDR:
[`calc_emp_fdr_breaks()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr_breaks.md),
[`plot_emp_fdr()`](https://stufield.github.io/stabilityselectr/dev/reference/plot_emp_fdr.md)

## Author

Stu Field, Michael R. Mehan

## Examples

``` r
withr::with_seed(101, {
  n_feat      <- 20
  n_samples   <- 100
  x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
  colnames(x) <- paste0("feat", "_", head(letters, n_feat))
  y  <- sample(1:2, n_samples, replace = TRUE)
})
ss <- stability_selection(x, y, "l1-logistic", num_iter = 25,
                          num_perms = 25, r_seed = 101, parallel = TRUE)
#> ✓ Using kernel: 'l1-logistic' and 1 core (serial)
calc_emp_fdr(ss, seq(0.5, 0.9, 0.1))
#> thresh_0.5 thresh_0.6 thresh_0.7 thresh_0.8 thresh_0.9 
#>      20.00      20.00      19.56      15.20       4.80 
```
