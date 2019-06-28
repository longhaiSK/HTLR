#ifndef GIBBS_H
#define GIBBS_H 

#include "RcppArmadillo.h"
#include "utils.h"
#include "sampler.h"

class Fit 
{
  private:

  // data
  const int p, K, C, n;
  arma::mat X; arma::mat ymat; arma::vec ybase;

  // prior
  std::string ptype; 
  const double alpha, s, eta, sigmab0;

  // sampling
  const int iters_rmc, iters_h, thin, leap_L, leap_L_h; 
  const double leap_step, hmc_sgmcut; 
  arma::vec DDNloglike;

  // fit result
  arma::cube mcdeltas;
  arma::mat mcsigmasbt, mcvardeltas;
  arma::vec mclogw, mcloglike, mcuvar, mchmcrej;

  // other control or result
  const int silence, looklf;

  // internal
  int nvar, nuvar, nfvar;
  int *ids_update, *ids_fix;

  arma::mat 
      lv, lv_old, lv_fix, normlv,
      predprob, predprob_old,
      DNloglike, DNloglike_old, 
      deltas, deltas_old, momt,
      DNlogprior, DNlogprior_old, DNlogpost;  
  
  arma::vec 
      SUMsqdeltas, SUMsqdeltas_old,
      SUMdeltas, SUMdeltas_old,
      vardeltas, vardeltas_old,
      stepsizes, sigmasbt; 
  double loglike, loglike_old, sgmsqcut;
 
  void whichupdate(double cut);  
  void updatepredprob();
  void detach_fixlv();
  void updateDNloglike();
  void updateloglike();
  void updateDNlogprior();
  void updateDNlogpost();
  void updatevardeltas();
  double comp_nenergy();
  void gen_momt();
  void cache_oldvalues();
  void restore_oldvalues();
  void updatestepsizes();


  public:

  Fit(int p, int K, int n,
      arma::mat X, arma::mat ymat, arma::vec ybase,
      std::string ptype, double alpha, double s, double eta, double sigmab0,
      int iters_rmc, int iters_h, int thin, 
      int leap_L, int leap_L_h, double leap_step,
      double hmc_sgmcut, arma::vec DDNloglike,
      arma::cube mcdeltas, arma::vec mclogw, arma::mat mcsigmasbt,
      int silence, int looklf);

  void StartSampling();
  Rcpp::List OutputR(); 

};

#endif
