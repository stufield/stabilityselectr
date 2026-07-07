# Introduction to Stability Selection

The `stabilityselectr` package performs stability selection with a
variety of kernels provided by the `glmnet` package, and provides simple
tools for plotting and extracting selected features. There is additional
functionality designed to facilitate various forms of permutation
clustering analyses.

------------------------------------------------------------------------

## Useful functions in `stabilityselectr`

- [`stability_selection()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_selection.md)
- [`is_stab_sel()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_selection.md)
- [`get_stable_features()`](https://stufield.github.io/stabilityselectr/dev/reference/get_stable_features.md)
- [`calc_emp_fdr()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr.md)
- [`calc_emp_fdr_breaks()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr.md)
- [`plot_emp_fdr()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr.md)
- [`plot_permuted_data()`](https://stufield.github.io/stabilityselectr/dev/reference/calc_emp_fdr.md)
- [`progeny_cluster()`](https://stufield.github.io/stabilityselectr/dev/reference/progeny_cluster.md)
- [`stability_cluster()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_cluster.md)

------------------------------------------------------------------------

## Examples

A typical stability selection analysis might be similar to the one
below.

#### Run `stability_selection()`

``` r

set.seed(101)
n_feat      <- 20L
n_samples   <- 100L
x           <- matrix(rnorm(n_feat * n_samples), n_samples, n_feat)
colnames(x) <- paste0("feat", "_", head(letters, n_feat))
y           <- sample(1:2, n_samples, replace = TRUE)
stab_sel    <- stability_selection(x, y, "binomial", n_iter = 500L)
#> ✓ Using kernel: 'binomial' and 1 core (serial)
#> ✓ Stablity path run time: 1.158s
is_stab_sel(stab_sel)
#> [1] TRUE
stab_sel
#> ══ Stability Selection (Kernel: binomial) ═════════════════════════════
#> • Weakness (alpha)            0.5
#> • Weakness Probability (Pw)   0.2
#> • Number of Iterations        500
#> • Standardized                'Yes'
#> • Imputed Outliers            'No'
#> • Lambda Max                  0.2261
#> • Lambda Min Ratio            0.1
#> • Permuted Data               'No'
#> • Random Seed                 1234
#> ═══════════════════════════════════════════════════════════════════════
```

#### Plot stability paths

``` r

plot(stab_sel, thresh = 0.85)
```

![](stabilityselectr_files/figure-html/plot-stab-sel-1.png)

#### Stable Features at a threshold

``` r

get_stable_features(stab_sel, thresh = 0.85)
#> $thresh_0.85
#> # A tibble: 4 × 3
#>   feature MaxSelectProb FDRbound
#>   <chr>           <dbl>    <dbl>
#> 1 feat_t          0.917  0.00357
#> 2 feat_d          0.892  0.00714
#> 3 feat_m          0.871  0.0107 
#> 4 feat_j          0.866  0.0143

# at multiple thresholds
get_stable_features(stab_sel, thresh = seq(0.7, 0.9, 0.05))
#> $thresh_0.7
#> # A tibble: 20 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_t          0.917  0.00625
#>  2 feat_d          0.892  0.0125 
#>  3 feat_m          0.871  0.0188 
#>  4 feat_j          0.866  0.025  
#>  5 feat_s          0.846  0.0313 
#>  6 feat_a          0.824  0.0375 
#>  7 feat_g          0.823  0.0438 
#>  8 feat_r          0.806  0.05   
#>  9 feat_c          0.795  0.0563 
#> 10 feat_q          0.795  0.0625 
#> 11 feat_l          0.79   0.0688 
#> 12 feat_f          0.789  0.075  
#> 13 feat_n          0.787  0.0813 
#> 14 feat_e          0.761  0.0875 
#> 15 feat_h          0.751  0.0938 
#> 16 feat_b          0.747  0.1    
#> 17 feat_p          0.747  0.106  
#> 18 feat_o          0.74   0.113  
#> 19 feat_i          0.739  0.119  
#> 20 feat_k          0.718  0.125  
#> 
#> $thresh_0.75
#> # A tibble: 15 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_t          0.917    0.005
#>  2 feat_d          0.892    0.01 
#>  3 feat_m          0.871    0.015
#>  4 feat_j          0.866    0.02 
#>  5 feat_s          0.846    0.025
#>  6 feat_a          0.824    0.03 
#>  7 feat_g          0.823    0.035
#>  8 feat_r          0.806    0.04 
#>  9 feat_c          0.795    0.045
#> 10 feat_q          0.795    0.05 
#> 11 feat_l          0.79     0.055
#> 12 feat_f          0.789    0.06 
#> 13 feat_n          0.787    0.065
#> 14 feat_e          0.761    0.07 
#> 15 feat_h          0.751    0.075
#> 
#> $thresh_0.8
#> # A tibble: 8 × 3
#>   feature MaxSelectProb FDRbound
#>   <chr>           <dbl>    <dbl>
#> 1 feat_t          0.917  0.00417
#> 2 feat_d          0.892  0.00833
#> 3 feat_m          0.871  0.0125 
#> 4 feat_j          0.866  0.0167 
#> 5 feat_s          0.846  0.0208 
#> 6 feat_a          0.824  0.025  
#> 7 feat_g          0.823  0.0292 
#> 8 feat_r          0.806  0.0333 
#> 
#> $thresh_0.85
#> # A tibble: 4 × 3
#>   feature MaxSelectProb FDRbound
#>   <chr>           <dbl>    <dbl>
#> 1 feat_t          0.917  0.00357
#> 2 feat_d          0.892  0.00714
#> 3 feat_m          0.871  0.0107 
#> 4 feat_j          0.866  0.0143 
#> 
#> $thresh_0.9
#> # A tibble: 1 × 3
#>   feature MaxSelectProb FDRbound
#>   <chr>           <dbl>    <dbl>
#> 1 feat_t          0.917  0.00313
```

------------------------------------------------------------------------

## Progeny and Stability Clustering

See separate vignette on clustering:
[`vignette("progeny-clustering")`](https://stufield.github.io/stabilityselectr/dev/articles/progeny-clustering.md).

------------------------------------------------------------------------

## References

Meinshausen, N and P Buhlmann. (2010). Stability selection. Journal of
the Royal Statistical Society: Series B (Statistical Methodology),
**72**: 417-473. doi: 10.1111/j.1467-9868.2010.00740.x
