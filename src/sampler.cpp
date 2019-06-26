#include "RcppArmadillo.h"
#include "ars.h"

class SamplerSgm : public SampleTarget
{
  protected:

  int idx;
  const int p, K;
  const double alpha, log_aw;
  const arma::vec vardeltas;

  public:

  SamplerSgm(int p, arma::vec vardeltas, int K, double alpha, double log_aw)
      : p(p), vardeltas(vardeltas), K(K), alpha(alpha), log_aw(log_aw)
  {}

  void set_idx(int i)
  {
    idx = i;
  }

};

class SamplerSgmNeg : public SamplerSgm
{
  public:

  SamplerSgmNeg(int p, arma::vec vardeltas, int K, double alpha, double log_aw)
      : SamplerSgm(p, vardeltas, K, alpha, log_aw)   
  {}

  void eval_logf(const double x, double &logf, double &dlogf)
  {
    logf = -(K / 2.0 - 1.0) * x;
    dlogf = -(K / 2.0 - 1.0);

    double vexi = vardeltas[idx] / 2.0 / exp(x);
    logf -= vexi;
    dlogf += vexi;

    double eximinlogaw = exp(x - log_aw);
    logf -= (alpha / 2.0 + 1.0) * log(1.0 + 2.0 * eximinlogaw);
    dlogf -= (alpha + 2.0) * eximinlogaw / (1.0 + 2.0 * eximinlogaw);
  }
};

class SamplerSgmGhs : public SamplerSgm
{
  public:

  SamplerSgmGhs(int p, arma::vec vardeltas, int K, double alpha, double log_aw)
      : SamplerSgm(p, vardeltas, K, alpha, log_aw)   
  {}

  void eval_logf(const double x, double &logf, double &dlogf)
  {
    logf  = - (K  - 1.0) / 2.0 * x ;
    dlogf = - (K  - 1.0) / 2.0 ;

    double vexi = vardeltas[idx] / 2.0 / exp(x);
    logf -= vexi;
    dlogf += vexi;

    double eximinlogaw = exp(x - log_aw);
    logf -= (alpha + 1.0) / 2.0 * log (1.0 + eximinlogaw);
    dlogf -= (alpha + 1.0) / 2.0 * eximinlogaw / (1.0 + eximinlogaw);
  }
};

class SamplerLogw : public SampleTarget
{
  protected:

  const int p, K;
  const double nu, s, eta;
  const arma::vec vardeltas;

  public:

  SamplerLogw(int p, arma::vec vardeltas, int K,
              double nu, double s, double eta)
      : p(p), vardeltas(vardeltas), K(K), nu(nu), s(s), eta(eta)
  {
  }

  void eval_logf(const double x, double &logf, double &dlogf)
  {
    double w = exp(x);
    double sdu = (x - s) / eta;
    double wnu = w * nu;
    // loglikelihood
    dlogf = 0;
    logf = 0;
    for (int j = 0; j < p; j++)
    {
      logf += log(vardeltas[j] + wnu);
      dlogf += wnu / (vardeltas[j] + wnu);
    }
    logf *= (-(nu + K) / 2);
    dlogf *= (-(nu + K) / 2);
    logf += p * nu / 2 * x;
    dlogf += p * nu / 2;
    logf += -R_pow_di(sdu, 2) / 2 - log(eta);
    dlogf += -sdu / eta;
  }

};

