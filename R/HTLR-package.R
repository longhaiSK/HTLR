#' Bayesian Logistic Regression with Heavy-Tailed Priors
#'
#' Efficient Bayesian multinomial logistic regression based on heavy-tailed priors. 
#' This package is suitable for classification and feature selection with high-dimensional 
#' features, such as gene expression profiles. Heavy-tailed priors can impose stronger 
#' shrinkage (compared to Gaussian and Laplace priors) to the coefficients associated with 
#' a large number of useless features, but still allow coefficients of a small number of 
#' useful features to stand out without punishment. It can also automatically make selection 
#' within a large number of correlated features. The posterior of coefficients and hyper-
#' parameters is sampled with restricted Gibbs sampling for leveraging high-dimensionality 
#' and Hamiltonian Monte Carlo for handling high-correlations among coefficients.
#'
#' @references
#' Longhai Li and Weixin Yao (2018). Fully Bayesian Logistic Regression 
#' with Hyper-Lasso Priors for High-dimensional Feature Selection.
#' \emph{Journal of Statistical Computation and Simulation} 2018, 88:14, 2827-2851.
#'
#' @name HTLR-package
#' 
#' @docType package
#' 
NULL
