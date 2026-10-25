# Repository guide

HTLR is an R package for Bayesian logistic regression, with Rcpp/RcppArmadillo C++ code.

- R code lives in `R/`, native code in `src/`, and tests in `tests/testthat/`.
- Follow existing style and keep changes focused.
- Edit roxygen comments in `R/` and regenerate documentation with `devtools::document()`; do not hand-edit generated `man/` files or `NAMESPACE`.
- Edit `README.Rmd`, the source for `README.md`.
- Run tests from the repository root: `Rscript -e 'devtools::test()'`.
- For package-wide validation, run `Rscript -e 'devtools::check()'`.
- Add regression tests for bug fixes; use fixed seeds for stochastic tests.
