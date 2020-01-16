
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HTLR: Bayesian Logistic Regression with Heavy-tailed Priors

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/HTLR)](https://CRAN.R-project.org/package=HTLR)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://cranlogs.r-pkg.org/badges/HTLR)](https://cran.r-project.org/package=HTLR)
[![Build
Status](https://travis-ci.org/longhaiSK/HTLR.svg?branch=master)](https://travis-ci.org/longhaiSK/HTLR)
<!-- badges: end -->

`HTLR` performs classification and feature selection by fitting Bayesian
polychotomous (multiclass, multinomial) logistic regression models based
on heavy-tailed priors with small degree freedom. This package is
suitable for classification with high-dimensional features, such as gene
expression profiles. Heavy-tailed priors can impose stronger shrinkage
(compared to Guassian and Laplace priors) to the coefficients associated
with a large number of useless features, but still allow coefficients of
a small number of useful features to stand out with little punishment.
Heavy-tailed priors can also automatically make selection within a large
number of correlated features. The posterior of coefficients and
hyperparameters is sampled with resitricted Gibbs sampling for
leveraging high-dimensionality and Hamiltonian Monte Carlo for handling
high-correlations among coefficients.

## Installation

[CRAN](https://CRAN.R-project.org) version (recommended):

``` r
install.packages("HTLR")
```

Development version on [GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("longhaiSK/HTLR")
```

This package uses library [Armadillo](http://arma.sourceforge.net/) for
carrying out most of matrix operations, you may get speed benefits from
using an alternative BLAS library such as
[ATLAS](http://math-atlas.sourceforge.net/),
[OpenBLAS](http://www.openblas.net/) or [Intel
MKL](https://software.intel.com/en-us/mkl). Check this
[post](http://brettklamer.com/diversions/statistical/faster-blas-in-r/)
for the comparision and the installation guide. Windows users may
consider installing [Microsoft R Open](https://mran.microsoft.com/open).

## Reference

Longhai Li and Weixin Yao (2018). Fully Bayesian Logistic Regression
with Hyper-Lasso Priors for High-dimensional Feature Selection.  2018,
88:14, 2827-2851, [the published
version](https://www.tandfonline.com/doi/full/10.1080/00949655.2018.1490418),
or [arXiv version](https://arxiv.org/pdf/1405.3319.pdf).
