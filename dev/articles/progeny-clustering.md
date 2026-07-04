# Progeny and Stability Clustering

## Useful functions:

- [`progeny_cluster()`](https://stufield.github.io/stabilityselectr/dev/reference/progeny_cluster.md):
  performs progeny clustering
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`print()`](https://rdrr.io/r/base/print.html): S3 methods for class
  `pclust`
- [`stability_cluster()`](https://stufield.github.io/stabilityselectr/dev/reference/stability_cluster.md):
  performs stability clustering

------------------------------------------------------------------------

## Progeny Clustering via `progeny_cluster()`

Select the optimal number for clustering using Progeny Clustering. The
“true” number of clusters in the `progeny_data` object is 3.

``` r

pc <- progeny_cluster(progeny_data, clust_iter = 2:9L,
                      repeats = 10L, n_iter = 25L, size = 6)
pc
#> ══ Progeny Clustering ═════════════════════════════════════════════════
#> • Call              'progeny_cluster(data = progeny_data, clust_iter = 2:9L, repeats = 10L, n_iter = 25L, size = 6)'
#> • Progeny Size      6
#> • K iterations      '2 → 9'
#> • No. iterations    25
#> • No. repeats       10
#> ── Mean & CI95 Stability Scores ───────────────────────────────────────
#>        k=2  k=3   k=4  k=5  k=6  k=7  k=8 ★k=9★
#> 2.5%  2.80 14.6  9.52 13.3 11.2 15.3 17.3  19.5
#> mean  3.57 20.4 12.42 16.1 13.2 17.6 19.6  24.6
#> 97.5% 4.63 35.7 18.52 20.3 15.4 19.9 21.0  31.6
#> ── Maximum Distance Scores ────────────────────────────────────────────
#>    k=2  ★k=3★    k=4    k=5    k=6    k=7    k=8    k=9 
#>  0.401 13.009  4.108  4.697 -3.677 -7.018  2.029  5.034
#> ── Gap Distance Scores ────────────────────────────────────────────────
#>    k=2  ★k=3★    k=4    k=5    k=6    k=7    k=8    k=9 
#> -24.83  24.83 -11.63   6.49  -7.18   2.27  -2.96   2.96
#> ═══════════════════════════════════════════════════════════════════════
```

``` r

plot(pc)
```

![](progeny-clustering_files/figure-html/plot-pclust-1.png)

------------------------------------------------------------------------

## Stability Clustering via `stability_cluster()`

Partitioning Around Medoids (PAM) is used both because is uses a more
robust measurement of the cluster centers (medoids) and because this
implementation keeps the cluster labels consistent across runs, a key
feature in calculating the across run stability. This does not occur
using [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html) where the
initial cluster labels are arbitrarily assigned.

True clusters are:

- cluster 1 -\> samples `1:50`
- cluster 2 -\> samples `51:100`
- cluster 3 -\> samples `101:150`

``` r

stab_clust <- stability_cluster(progeny_data, k = 3L,
                                n_iter = 500L, r_seed = 999)

stab_clust$true_cluster <- rep(1:3L, each = 50L)

# view 3-way confusion matrix
stab_clust |>
  with(table(truth = true_cluster, predicted = prob_k))
#>      predicted
#> truth  1  2  3
#>     1 50  0  0
#>     2  0 47  3
#>     3  1  0 49

# identify false clusters
stab_clust <- stab_clust |>
  dplyr::mutate(
    sample = dplyr::row_number(),
    pch    = rep(16:18, each = 50),
    pch    = dplyr::case_when(
      sample <= 50 & prob_k != 1L ~ 13L,                 # cluster 1
      sample > 50 & sample <= 100 & prob_k != 2L ~ 13L,  # cluster 2
      sample > 100 & prob_k != 3L ~ 13L,                 # cluster 3
      TRUE ~ pch
    )
  )

# view incorrect clusters (n = 8)
stab_clust |>
  dplyr::filter(true_cluster != prob_k)
#> # A tibble: 4 × 7
#>   `k=1` `k=2` `k=3` prob_k true_cluster sample   pch
#>   <dbl> <dbl> <dbl>  <int>        <int>  <int> <int>
#> 1 0.182 0.318 0.5        3            2     58    13
#> 2 0.174 0.328 0.498      3            2     73    13
#> 3 0.182 0.364 0.454      3            2     86    13
#> 4 0.64  0.204 0.156      1            3    115    13
```

### Plotting Clusters

We can plot the `progeny_data` object, which has 3 main clusters, and
identify which samples that were “correctly” clustered via stability
clustering with an “X”.

``` r

withr::with_par(list(mgp = c(2.00, 0.75, 0.0),
                     mar = c(3, 4, 3, 1), mfrow = 1:2L), {
  plot(progeny_data, col = rep(2:4, each = 50L),
       pch = rep(16:18, each = 50), cex = 1.75,
       main = "Simulated 3 Cluster Data")
  plot(progeny_data, col = rep(2:4, each = 50),
       pch = stab_clust$pch, cex = 1.75,
       main = "Stability Clustering")
})
```

![](progeny-clustering_files/figure-html/plot-clusters-1.png)

------------------------------------------------------------------------

## References

Hu, C.W., Kornblau, S.M., Slater, J.H. and A.A. Qutub (2015). Progeny
Clustering: A Method to Identify Biological Phenotypes. Scientific
Reports, 5:12894. <http://www.nature.com/articles/srep12894>
