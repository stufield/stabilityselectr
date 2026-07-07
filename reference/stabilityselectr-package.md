# stabilityselectr: R Package to Perform Stability and Cluster Selection

The stabilityselectr package performs stability selection with a variety
of kernels provided by the 'glmnet' package, and provides simple tools
for plotting and extracting selected features. There is additional
functionality designed to facilitate various forms of permutation
clustering analyses.

## Details

Stability selection is performed using
[`stability_selection()`](https://stufield.github.io/stabilityselectr/reference/stability_selection.md),
which returns a `stab_sel` class object. The stability path can be
plotted using a S3
[`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)
method. A `tibble` with the highest selection probability can be created
using
[`get_stable_features()`](https://stufield.github.io/stabilityselectr/reference/get_stable_features.md).
See vignettes.

## References

Meinshausen, N and P Buhlmann (2010). Stability selection. Journal of
the Royal Statistical Society: Series B (Statistical Methodology),
**72**: 417-473. doi: 10.1111/j.1467-9868.2010.00740.x

Hu, CW, Kornblau, SM, Slater, JH and AA Qutub (*2015*). Progeny
Clustering: A Method to Identify Biological Phenotypes. Scientific
Reports, **5**:12894. <http://www.nature.com/articles/srep12894>

## See also

Useful links:

- <https://stufield.github.io/stabilityselectr>

[`glmnet::glmnet()`](https://glmnet.stanford.edu/reference/glmnet.html),
[`kmeans()`](https://rdrr.io/r/stats/kmeans.html)

## Author

**Maintainer**: Stu Field <stu.g.field@gmail.com>
([ORCID](https://orcid.org/0000-0002-1024-5859)) \[copyright holder\]
