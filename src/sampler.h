#ifndef SAMPLER_H
#define SAMPLER_H 

#include "RcppArmadillo.h"
#include "ars.h"

class SamplerSgm : public SampleTarget
{
  protected:

  int idx;
  const int p, K;
  const double alpha, log_aw;
  arma::vec vardeltas;

  public:

  SamplerSgm(int p, arma::vec vardeltas, int K, double alpha, double log_aw);
  void set_idx(int i);

};

class SamplerSgmNeg : public SamplerSgm
{
  public:

  SamplerSgmNeg(int p, arma::vec vardeltas, int K, double alpha, double log_aw);
  void eval_logf(const double x, double &logf, double &dlogf) override; 
};

class SamplerSgmGhs : public SamplerSgm
{
  public:

  SamplerSgmGhs(int p, arma::vec vardeltas, int K, double alpha, double log_aw);
  void eval_logf(const double x, double &logf, double &dlogf) override;
};

class SamplerLogw : public SampleTarget
{
  protected:

  const int p, K;
  const double nu, s, eta;
  const arma::vec vardeltas;

  public:

  SamplerLogw(int p, arma::vec vardeltas, int K,
              double nu, double s, double eta);

  void eval_logf(const double x, double &logf, double &dlogf) override;
};

#endif