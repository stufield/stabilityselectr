# Changelog

## stabilityselectr 0.2.0

- Bug fixes
  - fixed trapezoidal AUC calculation in the S3
    [`summary()`](https://rdrr.io/r/base/summary.html) method; previous
    AUCs were incorrect (also affected the plot routine)
  - fixed integer x and y axes in
    [`plot_emp_fdr()`](https://stufield.github.io/stabilityselectr/reference/calc_emp_fdr.md)
  - fixed unicode symbols in the S3 plot method via
    [`bquote()`](https://rdrr.io/r/base/bquote.html); previously failed
    during R CMD check
- Reproducibility
  - reworked seed handling for multi-core runs
  - [`stability_selection()`](https://stufield.github.io/stabilityselectr/reference/stability_selection.md)
    now sets the seed within `mclapply()`, as in
    [`stability_cluster()`](https://stufield.github.io/stabilityselectr/reference/stability_cluster.md)
    and
    [`progeny_cluster()`](https://stufield.github.io/stabilityselectr/reference/progeny_cluster.md);
    parallel and serial results now match
  - unit tests updated for the new seed behavior
- Changes
  - modified default `alpha` and `Pw` in
    [`stability_selection()`](https://stufield.github.io/stabilityselectr/reference/stability_selection.md)
    to match the docs
  - simplified `FDRbound` calculation in
    [`get_stable_features()`](https://stufield.github.io/stabilityselectr/reference/get_stable_features.md)
  - removed the scales dependency from the plot routine
  - cleaned up
    [`plot_emp_fdr()`](https://stufield.github.io/stabilityselectr/reference/calc_emp_fdr.md);
    fixed ggplot `linewidth` warning
  - cleaned up internal permutation checks
- Docs
  - cleaned up examples, references, vignettes, and spelling;
    re-rendered README
- CI
  - bumped GHA `actions/checkout` v3 -\> v7
  - corrected GHA runner to `macos-latest`

## stabilityselectr 0.1.0

- Development …
  - major structural changes
  - speed improvement
  - converge function parameters
  - general improved user experience

## stabilityselectr 0.0.2

- Release to public GitHub

## stabilityselectr 0.0.1.9000

- Structural package development …

## stabilityselectr 0.0.1 🎉

- Initial release! 🥳
