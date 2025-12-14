# Make Prediction on New Data

Similar to other predict methods, this function returns predictions from
a fitted `htlrfit` object.

## Usage

``` r
# S3 method for class 'htlr.fit'
predict(object, newx, type = c("response", "class"), ...)
```

## Arguments

- object:

  A fitted model object with S3 class `htlrfit`.

- newx:

  A Matrix of values at which predictions are to be made.

- type:

  Type of prediction required. Type "response" gives the fitted
  probabilities. Type "class" produces the class label corresponding to
  the maximum probability.

- ...:

  Advanced options to specify the Markov chain iterations used for
  inference. See
  [`htlr_predict`](https://longhaisk.github.io/HTLR/reference/htlr_predict.md).

## Value

The object returned depends on type.
