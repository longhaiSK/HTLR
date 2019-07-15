
<!-- README.md is generated from README.Rmd. Please edit that file -->
HTLR
====

<!-- badges: start -->
<!-- badges: end -->
This package performs classification and feature selection by fitting Bayesian polychotomous (multiclass) logistic regression models based on heavy-tailed priors with small degree freedom. The software is suitable for classification with high-dimensional features, such as gene expression profiles. Heavy-tailed priors can impose stronger shrinkage (compared to Guassian and Laplace priors) to the coefficients associated with a large number of useless features, but still allow coefficients of a small number of useful features to stand out without punishment. Heavy-tailed priors can also automatically make selection within a large number of correlated features. The posterior of coefficients and hyperparameters is sampled with resitricted Gibbs sampling for leveraging high-dimensionality and Hamiltonian Monte Carlo for handling high-correlations among coefficients. The core computation in this software is carried out with fast C++ code.

Installation
------------

~~You can install the released version of HTLR from [CRAN](https://CRAN.R-project.org) with:~~

``` r
install.packages("HTLR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("longhaiSK/HTLR")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(HTLR)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.svg" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
