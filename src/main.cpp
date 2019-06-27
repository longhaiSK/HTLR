#include "RcppArmadillo.h"
#include "gibbs.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
Rcpp::List HtlrFit(
  int p, int K, int n,
  arma::mat X, arma::mat ymat, arma::vec ybase,
  std::string ptype, double alpha, double s, double eta, double sigmab0,
  int iters_rmc, int iters_h, int thin, 
  int leap_L, int leap_L_h, double leap_step,
  double hmc_sgmcut, arma::vec DDNloglike,
  arma::cube mcdeltas, arma::vec mclogw, arma::mat mcsigmasbt,
  int silence, int looklf)
{
  auto f = Fit(
    p, K, n, X, ymat, ybase,
    ptype, alpha, s, eta, sigmab0,
    iters_rmc, iters_h, thin, 
    leap_L, leap_L_h, leap_step, hmc_sgmcut, DDNloglike, 
    mcdeltas, mclogw, mcsigmasbt,
    silence, looklf);

  f.StartSampling();
  return f.OutputR(); 
}