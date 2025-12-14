# Changelog

## HTLR 0.4-3

CRAN release: 2020-09-09

### New Features

- Added option for users to keep samples of warmup iterations.

### Improvements

- Bug fix \[[\#7](https://github.com/longhaiSK/HTLR/issues/7)\].

## HTLR 0.4-2

CRAN release: 2020-01-17

### New Features

- Added new function
  [`std()`](https://longhaisk.github.io/HTLR/reference/std.md) for
  feature standardization.

- Added new dataset `diabetes392`.

### Improvements

- [`htlr()`](https://longhaisk.github.io/HTLR/reference/htlr.md) and
  [`predict.htlr.fit()`](https://longhaisk.github.io/HTLR/reference/predict.htlr.fit.md)
  now handles non-matrix input, i.e. data.frame.

- Minor speed improvement on
  [`htlr()`](https://longhaisk.github.io/HTLR/reference/htlr.md) and
  [`gendata_FAM()`](https://longhaisk.github.io/HTLR/reference/gendata_FAM.md).

- Updated documentation of
  [`htlr()`](https://longhaisk.github.io/HTLR/reference/htlr.md).

### Note

- Changed package license from GPLv2 to GPLv3.

## HTLR 0.4-1

CRAN release: 2019-10-08

### Bug Fixes

- Fixed potential memory leak issue in ARS module.

## HTLR 0.4

CRAN release: 2019-10-06

### New Features

- This is the first released version of revamped HTLR.

- The Gibbs sampling routine is completely refactored using
  RcppArmadillo, which leads to a significant performance gain on
  multi-core/distributed machines.

- The fitted model object is registered to S3 class `htlrfit`, coming
  with a set of useful S3 methods
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html), and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html).

- New model fitting function
  [`htlr()`](https://longhaisk.github.io/HTLR/reference/htlr.md) has a
  more accessible interface, while
  [`htlr_fit()`](https://longhaisk.github.io/HTLR/reference/htlr_fit.md)
  and
  [`htlr_predict()`](https://longhaisk.github.io/HTLR/reference/htlr_predict.md)
  are still keeped for the best possible backward compatibility.

- Better cohesion with `bayesplot` and other packages of RStan
  toolchain.

- Added new dataset `colon`.

## HTLR 0.3

### Note

- This is HTLR originally created by Longhai Li with legacy version
  number 3.1-1.

- Not compatible with macOS.
