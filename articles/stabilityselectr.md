# Introduction to Stability Selection

The `stabilityselectr` package performs stability selection with a
variety of kernels provided by the `glmnet` package, and provides simple
tools for plotting and extracting selected features. There is additional
functionality designed to facilitate various forms of permutation
clustering analyses.

------------------------------------------------------------------------

## Useful functions in `stabilityselectr`

- [`stability_selection()`](https://stufield.github.io/stabilityselectr/reference/stability_selection.md)
- [`is_stab_sel()`](https://stufield.github.io/stabilityselectr/reference/stability_selection.md)
- [`get_stable_features()`](https://stufield.github.io/stabilityselectr/reference/get_stable_features.md)
- [`calc_emp_fdr()`](https://stufield.github.io/stabilityselectr/reference/calc_emp_fdr.md)
- [`calc_emp_fdr_breaks()`](https://stufield.github.io/stabilityselectr/reference/calc_emp_fdr.md)
- [`plot_emp_fdr()`](https://stufield.github.io/stabilityselectr/reference/calc_emp_fdr.md)
- [`plot_permuted_data()`](https://stufield.github.io/stabilityselectr/reference/calc_emp_fdr.md)
- [`progeny_cluster()`](https://stufield.github.io/stabilityselectr/reference/progeny_cluster.md)
- [`stability_cluster()`](https://stufield.github.io/stabilityselectr/reference/stability_cluster.md)

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
#> ✓ Stablity path run time: 1.709
is_stab_sel(stab_sel)
#> [1] TRUE
stab_sel
#> ══ Stability Selection (Kernel: binomial) ═════════════════════════════
#> • Weakness (alpha)            0.8
#> • Weakness Probability (Pw)   0.5
#> • Number of Iterations        500
#> • Standardized                'Yes'
#> • Imputed Outliers            'No'
#> • Lambda Max                  0.1471
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
#> # A tibble: 14 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.941  0.00357
#>  2 feat_t          0.927  0.00714
#>  3 feat_j          0.922  0.0107 
#>  4 feat_m          0.914  0.0143 
#>  5 feat_a          0.901  0.0179 
#>  6 feat_s          0.889  0.0214 
#>  7 feat_g          0.88   0.025  
#>  8 feat_n          0.867  0.0286 
#>  9 feat_r          0.865  0.0321 
#> 10 feat_e          0.862  0.0357 
#> 11 feat_h          0.859  0.0393 
#> 12 feat_c          0.858  0.0429 
#> 13 feat_p          0.856  0.0464 
#> 14 feat_q          0.856  0.05

# at multiple thresholds
get_stable_features(stab_sel, thresh = seq(0.7, 0.9, 0.05))
#> $thresh_0.7
#> # A tibble: 20 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.941  0.00625
#>  2 feat_t          0.927  0.0125 
#>  3 feat_j          0.922  0.0188 
#>  4 feat_m          0.914  0.025  
#>  5 feat_a          0.901  0.0313 
#>  6 feat_s          0.889  0.0375 
#>  7 feat_g          0.88   0.0438 
#>  8 feat_n          0.867  0.05   
#>  9 feat_r          0.865  0.0563 
#> 10 feat_e          0.862  0.0625 
#> 11 feat_h          0.859  0.0688 
#> 12 feat_c          0.858  0.075  
#> 13 feat_p          0.856  0.0813 
#> 14 feat_q          0.856  0.0875 
#> 15 feat_f          0.845  0.0938 
#> 16 feat_l          0.842  0.1    
#> 17 feat_i          0.832  0.106  
#> 18 feat_b          0.823  0.113  
#> 19 feat_k          0.823  0.119  
#> 20 feat_o          0.82   0.125  
#> 
#> $thresh_0.75
#> # A tibble: 20 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.941    0.005
#>  2 feat_t          0.927    0.01 
#>  3 feat_j          0.922    0.015
#>  4 feat_m          0.914    0.02 
#>  5 feat_a          0.901    0.025
#>  6 feat_s          0.889    0.03 
#>  7 feat_g          0.88     0.035
#>  8 feat_n          0.867    0.04 
#>  9 feat_r          0.865    0.045
#> 10 feat_e          0.862    0.05 
#> 11 feat_h          0.859    0.055
#> 12 feat_c          0.858    0.06 
#> 13 feat_p          0.856    0.065
#> 14 feat_q          0.856    0.07 
#> 15 feat_f          0.845    0.075
#> 16 feat_l          0.842    0.08 
#> 17 feat_i          0.832    0.085
#> 18 feat_b          0.823    0.09 
#> 19 feat_k          0.823    0.095
#> 20 feat_o          0.82     0.1  
#> 
#> $thresh_0.8
#> # A tibble: 20 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.941  0.00417
#>  2 feat_t          0.927  0.00833
#>  3 feat_j          0.922  0.0125 
#>  4 feat_m          0.914  0.0167 
#>  5 feat_a          0.901  0.0208 
#>  6 feat_s          0.889  0.025  
#>  7 feat_g          0.88   0.0292 
#>  8 feat_n          0.867  0.0333 
#>  9 feat_r          0.865  0.0375 
#> 10 feat_e          0.862  0.0417 
#> 11 feat_h          0.859  0.0458 
#> 12 feat_c          0.858  0.05   
#> 13 feat_p          0.856  0.0542 
#> 14 feat_q          0.856  0.0583 
#> 15 feat_f          0.845  0.0625 
#> 16 feat_l          0.842  0.0667 
#> 17 feat_i          0.832  0.0708 
#> 18 feat_b          0.823  0.075  
#> 19 feat_k          0.823  0.0792 
#> 20 feat_o          0.82   0.0833 
#> 
#> $thresh_0.85
#> # A tibble: 14 × 3
#>    feature MaxSelectProb FDRbound
#>    <chr>           <dbl>    <dbl>
#>  1 feat_d          0.941  0.00357
#>  2 feat_t          0.927  0.00714
#>  3 feat_j          0.922  0.0107 
#>  4 feat_m          0.914  0.0143 
#>  5 feat_a          0.901  0.0179 
#>  6 feat_s          0.889  0.0214 
#>  7 feat_g          0.88   0.025  
#>  8 feat_n          0.867  0.0286 
#>  9 feat_r          0.865  0.0321 
#> 10 feat_e          0.862  0.0357 
#> 11 feat_h          0.859  0.0393 
#> 12 feat_c          0.858  0.0429 
#> 13 feat_p          0.856  0.0464 
#> 14 feat_q          0.856  0.05   
#> 
#> $thresh_0.9
#> # A tibble: 5 × 3
#>   feature MaxSelectProb FDRbound
#>   <chr>           <dbl>    <dbl>
#> 1 feat_d          0.941  0.00313
#> 2 feat_t          0.927  0.00625
#> 3 feat_j          0.922  0.00938
#> 4 feat_m          0.914  0.0125 
#> 5 feat_a          0.901  0.0156
```

------------------------------------------------------------------------

## Progeny and Stability Clustering

See separate vignette on clustering:
[`vignette("progeny-clustering")`](https://stufield.github.io/stabilityselectr/articles/progeny-clustering.md).
