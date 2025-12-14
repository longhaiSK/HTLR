# Posterior Summaries

This function gives a summary of posterior of parameters.

## Usage

``` r
# S3 method for class 'htlr.fit'
summary(
  object,
  features = 1L:object$p,
  method = median,
  usedmc = get_sample_indice(dim(object$mcdeltas)[3], object$mc.param$iter.rmc),
  ...
)
```

## Arguments

- object:

  An object of S3 class `htlr.fit`.

- features:

  A vector of indices (int) or names (char) that specify the parameters
  we will look at. By default all parameters are selected.

- method:

  A function that is used to aggregate the MCMC samples. The default is
  `median`, other built-in/customized statistical functions such as
  `mean`, `sd`, and `mad` can also be used.

- usedmc:

  Indices of Markov chain iterations used for inference. By default all
  iterations are used.

- ...:

  Not used.

## Value

A point summary of MCMC samples.

## Examples

``` r
set.seed(12345)
data("colon")

fit <- htlr(X = colon$X, y = colon$y, fsel = 1:100, iter = 20)
summary(fit, features = 1:16)
#>                class 1
#> Intercept  0.938177120
#> V1         0.000000000
#> V2         0.000000000
#> V3         0.000000000
#> V4         0.000000000
#> V5         0.000000000
#> V6         0.000000000
#> V7         0.000000000
#> V8         0.000000000
#> V9         0.000000000
#> V10        0.000000000
#> V11        0.000000000
#> V12        0.000000000
#> V13        0.000000000
#> V14       -2.045512555
#> V15        0.002115235
#> V16       -0.019430491
#> attr(,"stats")
#> [1] "median"
  
```
