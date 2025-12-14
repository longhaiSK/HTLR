# Create a Matrix of Markov Chain Samples

The Markov chain samples (without warmup) included in a `htlr.fit`
object will be coerced to a matrix.

## Usage

``` r
# S3 method for class 'htlr.fit'
as.matrix(x, k = NULL, include.warmup = FALSE, ...)
```

## Arguments

- x:

  An object of S3 class `htlr.fit`.

- k:

  Coefficients associated with class `k` will be drawn. Must be a
  positive integer in 1,2,...,C-1 for C-class traning labels (base class
  0 can not be chosen). By default the last class is selected. For
  binary logistic model this argument can be ignored.

- include.warmup:

  Whether or not to include warmup samples

- ...:

  Not used.

## Value

A matrix with `(p + 1)` columns and `i` rows, where `p` is the number of
features excluding intercept, and `i` is the number of iterations after
burnin.

## Examples

``` r
## No. of features used: 100; No. of iterations after burnin: 15 
fit <- htlr(X = colon$X, y = colon$y, fsel = 1:100, iter = 20, warmup = 5)

dim(as.matrix(fit))
#> [1]  15 101
  
```
