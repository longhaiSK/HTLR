# Bias-corrected Bayesian classification initial state

Generate initial Markov chain state with Bias-corrected Bayesian
classification.

## Usage

``` r
bcbcsf_deltas(X, y, alpha = 0)
```

## Arguments

- X:

  Design matrix of traning data; rows should be for the cases, and
  columns for different features.

- y:

  Vector of class labels in training or test data set. Must be coded as
  non-negative integers, e.g., 1,2,...,C for C classes.

- alpha:

  The regularization proportion (between 0 and 1) for mixing the
  diagonal covariance estimates and the sample covariance estimated with
  the training samples. The default is 0, the covariance matrix is
  assumed to be diagonal, which is the most robust.

## Value

A matrix - the initial state of Markov Chain for HTLR model fitting.

## Details

Caveat: This method can be used only for continuous predictors such as
gene expression profiles, and it does not make sense for categorical
predictors such as SNP profiles.

## References

Longhai Li (2012). Bias-corrected hierarchical Bayesian classification
with a selected subset of high-dimensional features. *Journal of the
American Statistical Association*, 107(497), 120-134.

## See also

[`lasso_deltas`](https://longhaisk.github.io/HTLR/reference/lasso_deltas.md)
