#ifndef GIBBS_H
#define GIBBS_H 
#define ARMA_NO_DEBUG

#include "RcppArmadillo.h"
#include "utils.h"
#include "sampler.h"

class Fit 
{
  private:

  // data
  const int p_, K_, C_, n_;
  arma::mat X_; arma::mat ymat_; arma::uvec ybase_;

  // prior
  std::string ptype_; 
  const double alpha_, s_, eta_;

  // sampling
  const int iters_rmc_, iters_h_, thin_, leap_L_, leap_L_h_; 
  const double leap_step_; 
  double sgmsq_cut_;
  arma::vec DDNloglike_;

  // fit result
  arma::cube mc_deltas_;
  arma::mat mc_sigmasbt_, mc_var_deltas_;
  arma::vec mc_logw_, mc_loglike_, mc_uvar_, mc_hmcrej_;

  // other control or result
  const int silence_;
  const bool legacy_;

  // internal
  const int nvar_;
  double logw_;
  int nuvar_, nfvar_;
  arma::uvec ids_update_, ids_fix_, iup_;

  arma::mat 
      lv_, lv_old_, lv_fix_, norm_lv_,
      pred_prob_, pred_prob_old_,
      DNloglike_, DNloglike_old_, 
      deltas_, deltas_old_, momt_,
      DNlogprior_, DNlogprior_old_, DNlogpost_;  
  
  arma::vec 
      sumsq_deltas_, sumsq_deltas_old_,
      sum_deltas_, sum_deltas_old_,
      var_deltas_, var_deltas_old_,
      step_sizes_, sigmasbt_; 
  double loglike_, loglike_old_;
 
  void WhichUpdate(bool init = false);  
  arma::uvec GetIdsUpdate();
  arma::uvec GetIdsFix();
  void DetachFixlv();
  void UpdatePredProb();
  void UpdateDNlogLike();
  void UpdateLogLike();
  void UpdateDNlogPrior();
  void UpdateDNlogPost();
  void UpdateMomtAndDeltas();
  void UpdateVarDeltas();
  void UpdateSigmas();
  void UpdateSigmasT();
  void UpdateSigmasSgm(SamplerSgm *target);
  void UpdateSigmasNeg();
  void UpdateSigmasGhs();
  void UpdateLogw();
  double CompNegEnergy();
  void GenMomt();
  void MoveMomt();
  void UpdateStepSizes();
  void CacheOldValues();
  void RestoreOldValues();
  bool IsFault(double cri = 20);
  void Initialize();
  void Traject(int i_mc);

  public:

  Fit(int p, int K, int n,
      arma::mat &X, arma::mat &ymat, arma::uvec &ybase,
      std::string ptype, double alpha, double s, double eta,
      int iters_rmc, int iters_h, int thin, 
      int leap_L, int leap_L_h, double leap_step,
      double hmc_sgmcut, arma::vec &DDNloglike,
      arma::mat &deltas, double logw, arma::vec &sigmasbt,
      int silence, bool legacy);

  void StartSampling();
  Rcpp::List OutputR(); 

};

#endif
