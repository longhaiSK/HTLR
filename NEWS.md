# HTLR 1.0

* Changed the package maintainer.

* Changed the arXiv reference in the package description to DOI format.

# HTLR 0.5

* Set the minimum R version to 3.6.2 and the minimum Rcpp version to 1.0.0.

* Removed the explicit C++11 build setting and system requirement.

# HTLR 0.4-4

* Corrected the default ranking function in `lasso_deltas()` to `order_plain`.

# HTLR 0.4-3

* Added option for users to keep samples of warmup iterations.

* Bug fix [[#7](https://github.com/longhaiSK/HTLR/issues/7)].

# HTLR 0.4-2

* Added new function `std()` for feature standardization. 

* Added new dataset `diabetes392`.

* `htlr()` and `predict.htlr.fit()` now handles non-matrix input, i.e. data.frame.

* Minor speed improvement on `htlr()` and `gendata_FAM()`. 

* Updated documentation of `htlr()`.

## Note

* Changed package license from GPLv2 to GPLv3.

# HTLR 0.4-1

* Fixed potential memory leak issue in ARS module. 

# HTLR 0.4

* This is the first released version of revamped HTLR.

* The Gibbs sampling routine is completely refactored using RcppArmadillo, which leads to a significant performance gain on multi-core/distributed machines.  

* The fitted model object is registered to S3 class `htlrfit`, coming with a set of useful S3 methods `print()`, `summary()`, `predict()`, `as.matrix()`, and `nobs()`.   

* New model fitting function `htlr()` has a more accessible interface, while `htlr_fit()` and `htlr_predict()` are still keeped for the best possible backward compatibility.

* Better cohesion with `bayesplot` and other packages of RStan toolchain.

* Added new dataset `colon`. 

# HTLR 0.3

## Note

* This is HTLR originally created by Longhai Li with legacy version number 3.1-1.

* Not compatible with macOS. 
