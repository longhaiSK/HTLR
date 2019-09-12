
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HTLR

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/longhaiSK/HTLR.svg?branch=master)](https://travis-ci.org/longhaiSK/HTLR)
[![License: GPL
v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`HTLR` performs classification and feature selection by fitting Baeysian
polychotomous (multiclass, multinomial) logistic regression models based
on heavy-tailed priors with small degree freedom. This package is
suitable for classification with high-dimensional features, such as gene
expression profiles. Heavy-tailed priors can impose stronger shrinkage
(compared to Guassian and Laplace priors) to the coefficients associated
with a large number of useless features, but still allow coefficients of
a small number of useful features to stand out without punishment. It
can also automatically make selection within a large number of
correlated features. The posterior of coefficients and hyperparameters
is sampled with resitricted Gibbs sampling for leveraging
high-dimensionality and Hamiltonian Monte Carlo for handling
high-correlations among coefficients.

This site focuses mainly on illustrating the usage and syntax of `HTLR`.
For more details on the algorithm, see the original article:
\<[DOI:10.1080/00949655.2018.1490418](https://www.tandfonline.com/doi/full/10.1080/00949655.2018.1490418)\>
([PDF](https://arxiv.org/pdf/1405.3319.pdf)).

## Installation

You can install the released version of HTLR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("HTLR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("longhaiSK/HTLR")
```
