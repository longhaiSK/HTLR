# Generate Simulated Data with Multinomial Logistic Regression Model

This function generates the response variables `y` given optional
supplied `X` using a multinomial logistic regression model.

## Usage

``` r
gendata_MLR(n, p, NC = 3, nu = 2, w = 1, X = NULL, betas = NULL)
```

## Arguments

- n:

  Number of observations.

- p:

  Number of features.

- NC:

  Number of classes for response variables.

- nu, w:

  If `betas` is not supplied (default), the regression coefficients are
  generated with t prior with df = `nu`, scale = `sqrt(w)`; will be
  ignored if `betas` is supplied.

- X:

  The design matrix; will be generated from standard normal distribution
  if not supplied.

- betas:

  User supplied regression coefficients.

## Value

A list contains input matrix `X`, response variables `y`, and regression
coefficients `deltas`.

## See also

[`gendata_FAM`](https://longhaisk.github.io/HTLR/reference/gendata_FAM.md)

## Examples

``` r
set.seed(12345)
dat <- gendata_MLR(n = 100, p = 10)
ggplot2::qplot(dat$y, bins = 6)

corrplot::corrplot(cor(dat$X))

```
