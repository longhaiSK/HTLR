#include "RcppArmadillo.h"
#include "gibbs.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List HtlrFit(
  int p, int K, int n,
  arma::mat &X, arma::mat &ymat, arma::uvec &ybase,
  std::string ptype, double alpha, double s, double eta, double sigmab0,
  int iters_rmc, int iters_h, int thin, 
  int leap_L, int leap_L_h, double leap_step,
  double hmc_sgmcut, arma::vec &DDNloglike,
  arma::mat &deltas, double logw, arma::vec &sigmasbt,
  int silence, int looklf, bool legacy)
{
  auto f = Fit(
    p, K, n, X, ymat, ybase,
    ptype, alpha, s, eta, sigmab0,
    iters_rmc, iters_h, thin, 
    leap_L, leap_L_h, leap_step, hmc_sgmcut, DDNloglike, 
    deltas, logw, sigmasbt,
    silence, looklf, legacy);

  f.StartSampling();
  return f.OutputR(); 
}