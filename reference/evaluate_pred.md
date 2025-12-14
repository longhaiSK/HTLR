# Evaluate Prediction Results

This function compares the prediction results returned by a classifier
with ground truth, and finally gives a summary of the evaluation.

## Usage

``` r
evaluate_pred(y.pred, y.true, caseid = names(y.true), showplot = TRUE, ...)
```

## Arguments

- y.pred:

  A matrix of predicted probabilities, as returned by a classifier.

- y.true:

  Ground truth labels vector.

- caseid:

  The names of test cases which we take account of. By default all test
  cases are used for evaluation.

- showplot:

  Logical; if `TRUE`, a summary plot will be generated.

- ...:

  Not used.

## Value

A summary of evaluation result.
