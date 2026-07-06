# Perform Progeny Clustering

Determine the most stable (optimal) number of clusters via Progeny
Clustering algorithm.

The `is_pclust()` function checks whether an object is class `pclust`.
See [`inherits()`](https://rdrr.io/r/base/class.html).

## Usage

``` r
progeny_cluster(data, clust_iter = 2:10L, repeats = 10L, r_seed = 1234, ...)

is_pclust(x)

# S3 method for class 'pclust'
plot(x, ...)

# S3 method for class 'pclust'
print(x, ...)
```

## Arguments

- data:

  A \\n \times p\\ data matrix containing *n* samples and *p* features.
  Can also be a data frame where each row corresponds to a sample or
  observation, and each column corresponds to a feature or variable.

- clust_iter:

  `integer(n)`. Span of `k` clusters to interrogate over.

- repeats:

  `integer(1)`. The number of repeat iterations to perform. Particularly
  useful if error bars during plotting are desired.

- r_seed:

  `integer(1)`. Seed for the random number generator, for
  reproducibility.

- ...:

  Important! Parameters passed to the internal `progeny_k()`, i.e.
  `n_iter =` and `size =`. Or, for the
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method,
  arguments passed to the corresponding graphics device.

- x:

  A `pclust` class object (or an object to be tested for one).

## Value

A `pclust` class object, a list containing:

- scores:

  A matrix of stability scores for each iteration in a matrix, with `k`
  columns

- mean_scores:

  The mean stability scores for each cluster `k`

- ci95_scores:

  95% confidence interval scores

- random_scores:

  The reference (random) scores for each iteration at each clustering
  level (`k`)

- mean_random_scores:

  The mean of the reference (random) data set, i.e. column means of
  `random_scores`

- D_max:

  The distance between the mean stability scores and the mean reference
  scores for each cluster `k`

- D_gap:

  The "gap" distance metric for neighboring cluster k differences. See
  original paper for reference.

- clust_iter:

  Integer Sequence of `k` clusters interrogated

- repeats:

  The number of repeat iterations to performed

- n_iter:

  The number of progeny iterations to performed

- size:

  The progeny size used in each iteration

- call:

  The call made to `progeny_cluster()`

`is_pclust()` returns a logical boolean.

## References

Hu, C.W., Kornblau, S.M., Slater, J.H. and A.A. Qutub (2015). Progeny
Clustering: A Method to Identify Biological Phenotypes. Scientific
Reports, 5:12894. <http://www.nature.com/articles/srep12894>

## See also

[`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html)

Other cluster:
[`stability_cluster()`](https://stufield.github.io/stabilityselectr/reference/stability_cluster.md)

## Author

Stu Field

## Examples

``` r
# `n_iter =` and `size =` are passed to `progeny_k()`
pclust <- progeny_cluster(progeny_data, clust_iter = 2:9L,
                          n_iter = 20L, size = 6)
pclust
#> ══ Progeny Clustering ═════════════════════════════════════════════════
#> • Call              'progeny_cluster(data = progeny_data, clust_iter = 2:9L, n_iter = 20L, size = 6)'
#> • Progeny Size      6
#> • K iterations      '2 → 9'
#> • No. iterations    20
#> • No. repeats       10
#> 
#> ── Mean & CI95 Stability Scores ───────────────────────────────────────
#>        k=2  k=3  k=4  k=5  k=6  k=7  k=8 ★k=9★
#> 2.5%  2.40 13.4 10.2 12.6 10.8 15.0 16.2  21.8
#> mean  3.34 18.2 13.8 15.8 13.9 16.4 18.4  24.7
#> 97.5% 4.87 24.4 18.7 20.2 17.0 17.6 22.0  27.8
#> 
#> ── Maximum Distance Scores ────────────────────────────────────────────
#>    k=2  ★k=3★    k=4    k=5    k=6    k=7    k=8    k=9 
#>  0.666 11.403  4.447  2.235 -2.585  1.216  2.930  6.430 
#> 
#> ── Gap Distance Scores ────────────────────────────────────────────────
#>     k=2   ★k=3★     k=4     k=5     k=6     k=7     k=8     k=9 
#> -19.367  19.367  -6.464   3.913  -4.502   0.622  -4.341   4.341 
#> ═══════════════════════════════════════════════════════════════════════

# Test progeny clustering on iris data set
# Doesn't work as well as the simulated data set
clust_iris <- progeny_cluster(iris[, -5L], clust_iter = 2:5L,
                              size = 6L, n_iter = 250)  # pass `...`
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations
#> Warning: did not converge in 20 iterations

# true n_clusters = 3
clust_iris
#> ══ Progeny Clustering ═════════════════════════════════════════════════
#> • Call              'progeny_cluster(data = iris[, -5L], clust_iter = 2:5L, size = 6L, n_iter = 250)'
#> • Progeny Size      6
#> • K iterations      '2 → 5'
#> • No. iterations    250
#> • No. repeats       10
#> 
#> ── Mean & CI95 Stability Scores ───────────────────────────────────────
#>       ★k=2★  k=3  k=4  k=5
#> 2.5%    396 50.2 26.6 43.6
#> mean   1015 56.7 27.9 48.0
#> 97.5%  1499 65.1 30.3 52.5
#> 
#> ── Maximum Distance Scores ────────────────────────────────────────────
#> ★k=2★   k=3   k=4   k=5 
#> 969.3  48.7  19.6  38.3 
#> 
#> ── Gap Distance Scores ────────────────────────────────────────────────
#>  ★k=2★    k=3    k=4    k=5 
#>  927.7 -927.7  -48.9   48.9 
#> ═══════════════════════════════════════════════════════════════════════

# Test for class `pclust`
is_pclust(pclust)
#> [1] TRUE

# S3 plot method
plot(pclust)


plot(clust_iris)
```
