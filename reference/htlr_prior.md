# Generate Prior Configuration

Configure prior hyper-parameters for HTLR model fitting

## Usage

``` r
htlr_prior(
  ptype = c("t", "ghs", "neg"),
  df = 1,
  logw = -(1/df) * 10,
  eta = ifelse(df > 1, 3, 0),
  sigmab0 = 2000
)
```

## Arguments

- ptype:

  The prior to be applied to the model. Either "t" (student-t, default),
  "ghs" (horseshoe), or "neg" (normal-exponential-gamma).

- df:

  The degree freedom (aka alpha) of t/ghs/neg prior for coefficients.

- logw:

  The log scale of priors for coefficients.

- eta:

  The `sd` of the normal prior for logw. When it is set to 0, logw is
  fixed. Otherwise, logw is assigned with a normal prior and it will be
  updated during sampling.

- sigmab0:

  The `sd` of the normal prior for the intercept.

## Value

A configuration list containing `ptype`, `alpha`, `logw`, `eta`, and
`sigmab0`.

## Details

The output is a configuration list which is to be passed to `prior`
argument of `htlr`. For naive users, you only need to specify the prior
type and degree freedom, then the other hyper-parameters will be chosen
automatically. For advanced users, you can supply each prior
hyper-parameters by yourself. For suggestion of picking
hyper-parameters, see `references`.

## References

Longhai Li and Weixin Yao. (2018). Fully Bayesian Logistic Regression
with Hyper-Lasso Priors for High-dimensional Feature Selection. *Journal
of Statistical Computation and Simulation* 2018, 88:14, 2827-2851.
