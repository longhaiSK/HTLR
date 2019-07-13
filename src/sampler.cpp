#include "sampler.h"

SamplerSgm::SamplerSgm(int p, const arma::vec &vardeltas, int K, double alpha, double log_aw)
    : p_(p), K_(K), alpha_(alpha), log_aw_(log_aw), vardeltas_(vardeltas)
{
}

void SamplerSgm::set_idx(int i)
{
  idx_ = i;
}

SamplerSgmNeg::SamplerSgmNeg(int p, const arma::vec &vardeltas, int K, double alpha, double log_aw)
    : SamplerSgm(p, vardeltas, K, alpha, log_aw)
{
}

void SamplerSgmNeg::eval_logf(const double x, double &logf, double &dlogf)
{
  logf = -(K_ / 2.0 - 1.0) * x;
  dlogf = -(K_ / 2.0 - 1.0);

  double vexi = vardeltas_[idx_] / 2.0 / exp(x);
  logf -= vexi;
  dlogf += vexi;

  double eximinlogaw = exp(x - log_aw_);
  logf -= (alpha_ / 2.0 + 1.0) * log(1.0 + 2.0 * eximinlogaw);
  dlogf -= (alpha_ + 2.0) * eximinlogaw / (1.0 + 2.0 * eximinlogaw);
}

SamplerSgmGhs::SamplerSgmGhs(int p, const arma::vec &vardeltas, int K, double alpha, double log_aw)
    : SamplerSgm(p, vardeltas, K, alpha, log_aw)
{
}

void SamplerSgmGhs::eval_logf(const double x, double &logf, double &dlogf)
{
  logf = -(K_ - 1.0) / 2.0 * x;
  dlogf = -(K_ - 1.0) / 2.0;

  double vexi = vardeltas_[idx_] / 2.0 / exp(x);
  logf -= vexi;
  dlogf += vexi;

  double eximinlogaw = exp(x - log_aw_);
  logf -= (alpha_ + 1.0) / 2.0 * log(1.0 + eximinlogaw);
  dlogf -= (alpha_ + 1.0) / 2.0 * eximinlogaw / (1.0 + eximinlogaw);
}

SamplerLogw::SamplerLogw(int p, const arma::vec &vardeltas, int K,
                         double nu, double s, double eta)
    : p_(p), K_(K), nu_(nu), s_(s), eta_(eta), vardeltas_(vardeltas)
{
}

void SamplerLogw::eval_logf(const double x, double &logf, double &dlogf)
{
  double w = exp(x);
  double sdu = (x - s_) / eta_;
  double wnu = w * nu_;

  dlogf = arma::accu(wnu / (vardeltas_ + wnu));
  logf = arma::accu(arma::log(vardeltas_ + wnu));

  logf *= (-(nu_ + K_) / 2);
  dlogf *= (-(nu_ + K_) / 2);
  logf += p_ * nu_ / 2 * x;
  dlogf += p_ * nu_ / 2;
  logf += -R_pow_di(sdu, 2) / 2 - log(eta_);
  dlogf += -sdu / eta_;
}
