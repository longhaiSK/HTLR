# Order features by F-statistic

This function orders all features in terms of ANOVA F-statistic.

## Usage

``` r
order_ftest(X, y)
```

## Arguments

- X:

  Input matrix, of dimension `nobs` by `nvars`; each row is an
  observation vector.

- y:

  Vector of response variables.

## Value

Order of all features of length `nvars`.

## Examples

``` r
data("diabetes392")
order_ftest(diabetes392$X, diabetes392$y)
#> [1] 2 8 5 6 1 4 7 3
```
