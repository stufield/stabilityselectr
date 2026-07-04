# Perform Stability Clustering

Partitioning Around Medoids (PAM) is used both because it uses a more
robust measurement of the cluster centers (medoids) and because this
implementation keeps the cluster labels consistent across runs, a key
feature in calculating the across run stability. This does not occur
using [`kmeans()`](https://rdrr.io/r/stats/kmeans.html) where the
initial cluster labels are arbitrarily assigned.

## Usage

``` r
stability_cluster(data, k, n_iter = 100L, r_seed = 1234)
```

## Arguments

- data:

  A \\n \times p\\ data matrix containing *n* samples and *p* features.
  Can also be a data frame where each row corresponds to a sample or
  observation, and each column corresponds to a feature or variable.

- k:

  `integer(1)`. The number of clusters.

- n_iter:

  `integer(1)`. Defining the number of sub-sampling iterations during
  each selection.

- r_seed:

  `integer(1)`. Seed for the random number generator, for
  reproducibility.

## Value

A \\n \times (k + 1)\\ dimensional `tibble` of clustering probabilities
for each `k`, plus a final column named `prob_k`, which indicates the
"most probable" cluster membership for that sample.

## Note

How do we make sure clusters are indexed the same as what comes out of
k-means? Could be susceptible to index errors (but seems ok for now).

## References

Hastie, et al. 2009.

## See also

[`cluster::pam()`](https://rdrr.io/pkg/cluster/man/pam.html)

Other cluster:
[`progeny_cluster()`](https://stufield.github.io/stabilityselectr/dev/reference/progeny_cluster.md)

## Author

Stu Field

## Examples

``` r
stab_clust <- stability_cluster(progeny_data, k = 3L, n_iter = 750L)

stab_clust$true_cluster <- rep(1:3L, each = 50L)

# View stable clusters
stab_clust
#> # A tibble: 150 × 5
#>    `k=1` `k=2` `k=3` prob_k true_cluster
#>    <dbl> <dbl> <dbl>  <int>        <int>
#>  1 0.68  0.159 0.161      1            1
#>  2 0.669 0.16  0.171      1            1
#>  3 0.657 0.172 0.171      1            1
#>  4 0.675 0.169 0.156      1            1
#>  5 0.647 0.168 0.185      1            1
#>  6 0.663 0.185 0.152      1            1
#>  7 0.673 0.175 0.152      1            1
#>  8 0.653 0.188 0.159      1            1
#>  9 0.665 0.177 0.157      1            1
#> 10 0.66  0.184 0.156      1            1
#> # ℹ 140 more rows

# confusion matrix
stab_clust |>
  with(table(truth = true_cluster, predicted = prob_k))
#>      predicted
#> truth  1  2  3
#>     1 49  1  0
#>     2  0 47  3
#>     3  1  1 48

# View the incorrectly clustered samples (n = 5)
filter(stab_clust, prob_k != true_cluster)
#> # A tibble: 6 × 5
#>   `k=1` `k=2` `k=3` prob_k true_cluster
#>   <dbl> <dbl> <dbl>  <int>        <int>
#> 1 0.355 0.457 0.188      2            1
#> 2 0.153 0.407 0.44       3            2
#> 3 0.141 0.393 0.465      3            2
#> 4 0.189 0.36  0.451      3            2
#> 5 0.639 0.201 0.16       1            3
#> 6 0.171 0.416 0.413      2            3

# Plot Stability Clusters
cols <- c("#24135F", "#00A499", "#840B55")
withr::with_par(list(mgp = c(2.00, 0.75, 0.0),
                     mar = c(3, 4, 3, 1), mfrow = 1:2L), {
  plot(progeny_data,
       col = cols[stab_clust$true_cluster],
       bg  = cols[stab_clust$true_cluster],
       pch = stab_clust$true_cluster + 20,
       lwd = 1, cex = 1.75, main = "Simulated 3 Cluster Data")
  plot(progeny_data,
       col = cols[stab_clust$prob_k],
       bg  = cols[stab_clust$true_cluster],
       pch = stab_clust$true_cluster + 20,
       lwd = 2.5, cex = 1.5, main = "Stability Clustering")
})
```
