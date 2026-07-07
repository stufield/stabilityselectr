# Package Objects

Package Objects

## Usage

``` r
clust_data

progeny_data
```

## Format

`clust_data`: A 2 column data matrix containing the 2 features in the
example, named "F1" and "F2", and containing 20 observations.

`progeny_data`: A 2 column data matrix containing 2 features, "F1" and
"F2", 150 observations, and 3 clusters in bivariate space.

## Source

`clust_data`: Hu, *et al*. (2015).

`progeny_data`: Stu Field

## Functions

- `clust_data`: The original clustering example data from the Progeny
  Clustering paper. A 2 feature data set example (see references) in the
  example illustration in Figure 1, page 3.

- `progeny_data`: A simulated clustering data set generated to contain
  **3 true clusters** in bivariate space. There are 2 features, named
  "F1" and "F2" and 150 observations.

## References

`clust_data`: Hu, CW, Kornblau, SM, Slater, JH and AA Qutub (2015).
Progeny Clustering: A Method to Identify Biological Phenotypes.
Scientific Reports, **5**:12894. pg. 3.
<http://www.nature.com/articles/srep12894>

## See also

[`MASS::mvrnorm()`](https://rdrr.io/pkg/MASS/man/mvrnorm.html)

## Examples

``` r
head(clust_data)
#>     F1   F2
#> 1 0.40 0.85
#> 2 0.98 0.24
#> 3 0.35 0.98
#> 4 0.62 0.08
#> 5 0.48 0.79
#> 6 0.63 0.38

plot(progeny_data, col = rep(2:4, each = 50L),
     pch = rep(16:18, each = 50L), cex = 1.75,
     main = "Original Simulated 3 Cluster Data")
```
