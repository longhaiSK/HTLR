# Split Data into Train and Test Partitions

This function splits the input data and response variables into training
and testing parts.

## Usage

``` r
split_data(X, y, p.train = 0.7, n.train = round(nrow(X) * p.train))
```

## Arguments

- X:

  Input matrix, of dimension `nobs` by `nvars`; each row is an
  observation vector.

- y:

  Vector of response variables.

- p.train:

  Percentage of training set.

- n.train:

  Number of cases for training; will override `p.train` if specified.

## Value

List of training data `x.tr`, `y.tr` and testing data `x.te`, `y.te`.

## Examples

``` r
dat <- gendata_MLR(n = 100, p = 10)
dat <- split_data(dat$X, dat$y, p.train = 0.7)
dim(dat$x.tr)
#> [1] 70 10
dim(dat$x.te)
#> [1] 30 10
   
```
