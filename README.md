
<!-- README.md is generated from README.Rmd. Please edit that file -->

    ## [1] "en_US.UTF-8"

# The `stabilityselectr` package

<!-- badges: start -->

![GitHub
version](https://img.shields.io/badge/Version-0.0.1.9000-success.svg?style=flat&logo=github)
[![CRAN
status](http://www.r-pkg.org/badges/version/stabilityselectr)](https://cran.r-project.org/package=stabilityselectr)
[![R-CMD-check](https://github.com/stufield/stabilityselectr/workflows/R-CMD-check/badge.svg)](https://github.com/stufield/stabilityselectr/actions)
[![](https://cranlogs.r-pkg.org/badges/grand-total/stabilityselectr)](https://cran.r-project.org/package=stabilityselectr)
[![test
coverage](https://codecov.io/gh/stufield/stabilityselectr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/stufield/stabilityselectr?branch=main)
[![lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://choosealicense.com/licenses/mit/)
<!-- badges: end -->

## Overview

The `stabilityselectr` package performs stability selection with a
variety of kernels provided by the `glmnet` package, and provides simple
tools for plotting and extracting selected features. There is additional
functionality designed to facilitate various forms of permutation
clustering analyses.

------------------------------------------------------------------------

## Installation

``` r
# current dev version
remotes::install_github("stufield/stabilityselectr")

# or a specific version
remotes::install_github("stufield/stabilityselectr@v0.0.1")
```

------------------------------------------------------------------------

## Usage

To load `stabilityselectr` simply make a call to `library()` as usual:

``` r
library(stabilityselectr)
```

## Help

``` r
library(help = stabilityselectr)
```

## Package Notes

- The `stabilityselectr` package is easy to use.
- It is my first go-to when looking at dimensionality reduction and
  upstream feature selection.
