# Make Prediction on New Data (Advanced)

This function uses MCMC samples from fitted `htlrfit` object OR user
supplied regression coefficient to predict the class labels of test
cases.

## Usage

``` r
htlr_predict(
  X_ts,
  fithtlr = NULL,
  deltas = NULL,
  burn = NULL,
  thin = 1,
  usedmc = NULL,
  rep.legacy = TRUE
)
```

## Arguments

- X_ts:

  Matrix of values at which predictions are to be made.

- fithtlr:

  Fitted HTLR model object.

- deltas:

  The values of deltas (for example true deltas) used to make
  prediction; will override `fithtlr` if provided.

- burn, thin:

  `burn` of Markov chain (super)iterations will be discarded for
  prediction, and only every `thin` are used.

- usedmc:

  Indices of Markov chain iterations used for inference. If supplied,
  `burn` and `thin` will be ignored.

- rep.legacy:

  To reproduce (actually incorrect) results in legacy version. See
  <https://github.com/longhaiSK/HTLR/issues/7>.

## Value

A matrix of predictive probabilities, with rows for cases, cols for
classes.
