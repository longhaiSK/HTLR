# Standardizes a Design Matrix

This function accepts a design matrix and returns a standardized version
of that matrix, the statistics of each column such as `median` and `sd`
are also provided.

## Usage

``` r
std(X, tol = 1e-06)
```

## Arguments

- X:

  Design matrix, of dimension `nobs` by `nvars`; each row is an
  observation vector; can also be an object that can be coerced to a
  matrix, e.g. a data.frame.

- tol:

  The tolerance value; a column of `X` is considered as singular if the
  `sd` of its entries (observations) is less than `tol`. Singular
  columns will be dropped by the end.

## Value

The standardized design matrix with the following attributes:

- nonsingular:

  Indices of non-singular columns.

- center:

  Median of each non-singular column which is used for standardization.

- scale:

  Standard deviation of each non-singular column which is used for
  standardization.

## Details

For each column of `X`, the standardization is done by first subtracting
its median, then dividing by its sample standard deviation, while the
original version in `ncvreg` uses mean and population standard
deviation. Its speed is slower than `ncvreg` because of the complexity
of median finding, but still substantially faster than
[`scale()`](https://rdrr.io/r/base/scale.html) provided by R base.

## See also

[http://pbreheny.github.io/ncvreg/reference/std.html](http://pbreheny.github.io/ncvreg/reference/std.md)

## Author

Patrick Breheny (original)  
Steven Liu (modification)

## Examples

``` r
set.seed(123)
mat <- matrix(rnorm(n = 80 * 90, mean = 100, sd = 50), 80, 90)
mat %>% as.numeric() %>% ggplot2::qplot(bins = 30, xlab = '')

mat %>% std() %>% as.numeric() %>% ggplot2::qplot(bins = 30, xlab = '')

 
```
