#ifndef SAMPLER_H
#define SAMPLER_H 

#include "RcppArmadillo.h"
#include "ars.h"

class SamplerSgm : public SampleTarget
{
  protected:

  int idx_;
  const int p_, K_;
  const double alpha_, log_aw_;
  const arma::vec vardeltas_;

  public:

  SamplerSgm(int p, const arma::vec &vardeltas, int K, double alpha, double log_aw);
  void set_idx(int i);
};

class SamplerSgmNeg : public SamplerSgm
{
  public:

  SamplerSgmNeg(int p, const arma::vec &vardeltas, int K, double alpha, double log_aw);
  void eval_logf(const double x, double &logf, double &dlogf) override; 
};

class SamplerSgmGhs : public SamplerSgm
{
  public:

  SamplerSgmGhs(int p, const arma::vec &vardeltas, int K, double alpha, double log_aw);
  void eval_logf(const double x, double &logf, double &dlogf) override;
};

class SamplerLogw : public SampleTarget
{
  protected:

  const int p_, K_;
  const double nu_, s_, eta_;
  const arma::vec vardeltas_;

  public:

  SamplerLogw(int p, const arma::vec &vardeltas, int K,
              double nu, double s, double eta);

  void eval_logf(const double x, double &logf, double &dlogf) override;
};

#endif
