# HTLR 0.4

## New Features

* This is the first released version of revamped HTLR.

* The Gibbs sampling routine is completely refactored using RcppArmadillo, which leads to a significant performance gain on multi-core/distributed machines.  

* The fitted model object is registered to S3 class `htlrfit`, coming with a set of useful S3 methods `print`, `summary`, `predict`, `as.matrix`, and `nobs`.   

* New model fitting function `htlr` has a more accessible interface, while `htlr_fit` and `htlr_predict` are still keeped for the best possible backward compatibility.

* Better cohesion with `bayesplot` and other packages of RStan toolchain.

* Added new dataset `colon`. 

# HTLR 0.3.1

* This is HTLR originally created by Longhai Li with legacy version number 3.1-1.

* Not compatible with macOS. 
