# Order features by Kruskal-Wallis test

This function orders all features in terms of Kruskal-Wallis test
p-value.

## Usage

``` r
order_kruskal(X, y)
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
order_kruskal(diabetes392$X, diabetes392$y)
#> [1] 2 8 5 6 4 1 7 3
```
