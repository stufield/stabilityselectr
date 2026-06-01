# The `stabilityselectr` package

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
remotes::install_github("stufield/stabilityselectr@v0.0.2")
```

------------------------------------------------------------------------

## Usage

To load `stabilityselectr` simply make a call to
[`library()`](https://rdrr.io/r/base/library.html) as usual:

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
