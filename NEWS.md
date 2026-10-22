# stabilityselectr 0.2.0

* Bug fixes
  - fixed trapezoidal AUC calculation in the S3
    `summary()` method; previous AUCs were incorrect
    (also affected the plot routine)
  - fixed integer x and y axes in `plot_emp_fdr()`
  - fixed unicode symbols in the S3 plot method via
    `bquote()`; previously failed during R CMD check
* Reproducibility
  - reworked seed handling for multi-core runs
  - `stability_selection()` now sets the seed within
    `mclapply()`, as in `stability_cluster()` and
    `progeny_cluster()`; parallel and serial results
    now match
  - unit tests updated for the new seed behavior
* Changes
  - modified default `alpha` and `Pw` in
    `stability_selection()` to match the docs
  - simplified `FDRbound` calculation in
    `get_stable_features()`
  - removed the scales dependency from the plot routine
  - cleaned up `plot_emp_fdr()`; fixed ggplot
    `linewidth` warning
  - cleaned up internal permutation checks
* Docs
  - cleaned up examples, references, vignettes,
    and spelling; re-rendered README
* CI
  - bumped GHA `actions/checkout` v3 -> v7
  - corrected GHA runner to `macos-latest`

# stabilityselectr 0.1.0

* Development ...
  - major structural changes
  - speed improvement
  - converge function parameters
  - general improved user experience

# stabilityselectr 0.0.2

* Release to public GitHub

# stabilityselectr 0.0.1.9000

* Structural package development ...

# stabilityselectr 0.0.1 :tada:

* Initial release! :partying_face:
