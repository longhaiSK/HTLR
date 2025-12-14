# Get Indices of Non-Zero Coefficients

Get the indices of non-zero coefficients from fitted HTLR model objects.

## Usage

``` r
nzero_idx(fit, cut = 0.1)
```

## Arguments

- fit:

  An object of S3 class `htlr.fit`.

- cut:

  Threshold on relative SDB to distinguish zero coefficients.

## Value

Indices vector of non-zero coefficients in the model.

## Examples

``` r
set.seed(12345)
data("colon")

fit <- htlr(X = colon$X, y = colon$y, fsel = 1:100, iter = 20)
nzero_idx(fit)
#> [1] 14 43
```
