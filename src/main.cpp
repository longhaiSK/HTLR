#include "RcppArmadillo.h"
#include "gibbs.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List htlr_fit_helper(
  int p, int K, int n,
  arma::mat &X, arma::mat &ymat, arma::uvec &ybase,
  std::string ptype, double alpha, double s, double eta,
  int iters_rmc, int iters_h, int thin, 
  int leap_L, int leap_L_h, double leap_step,
  double hmc_sgmcut, arma::mat &deltas, double logw, 
  arma::vec &sigmasbt, int silence, bool legacy)
{
  auto f = Fit(
    p, K, n, X, ymat, ybase,
    ptype, alpha, s, eta,
    iters_rmc, iters_h, thin, 
    leap_L, leap_L_h, leap_step, hmc_sgmcut, 
    deltas, logw, sigmasbt,
    silence, legacy);

  f.StartSampling();
  return f.OutputR(); 
}