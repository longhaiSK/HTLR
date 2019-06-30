#include "ars.h"

// [[Rcpp::export]]
Rcpp::NumericVector sample_tnorm_ars(
    const int n, const double lb, const double ub, const bool verbose = false)
{
  class TruncNormTarget : public SampleTarget
  {
    private: 
    
    const double lb_, ub_;

    public:

    TruncNormTarget(double lb, double ub) : lb_(lb), ub_(ub) {}

    void eval_logf(const double x, double &logf, double &dlogf) override
    {
      // if ub and lb cannot be found explicitely, this checking condition can be
      // replaced by another checking expression, eg. power (x,2) < 1
      if ((x > ub_) || (x < lb_))
      {
        logf = R_NegInf;
        dlogf = R_NaN;
        return;
      }
      logf = -pow(x, 2) / 2;
      dlogf = -x;
    }
  };

  double init_tpoint = 0;
  if (R_FINITE(lb) & R_FINITE(ub))
    init_tpoint = (lb + ub) / 2;
  if (R_FINITE(lb) & !R_FINITE(ub))
    init_tpoint = lb + 1;
  if (!R_FINITE(lb) & R_FINITE(ub))
    init_tpoint = ub - 1;
  if (!R_FINITE(lb) & !R_FINITE(ub))
    init_tpoint = 0;

  // Do adaptive rejection sampling here
  // Note that the bounds are set to -Inf and +Inf, so the actual bounds
  // are to be determined by the ARS sampler.
  auto target = TruncNormTarget(lb, ub);
  auto sampler = ARS(n, &target, init_tpoint, R_NegInf, R_PosInf, verbose);
  return sampler.Sample(); 
}

// [[Rcpp::export]]
Rcpp::NumericVector sample_post_ichi(
    const int n, const Rcpp::NumericVector &sigmasq, const double alpha1, 
    const double alpha0 = 1E-5, const double w0 = 1E-5, 
    const bool verbose = false)
{
  class IchiTarget : public SampleTarget
  {
    private: 
    
    double lambda_p_, alpha_p_, lambda0_;

    public:

    IchiTarget(Rcpp::NumericVector sigmasq, double alpha1, double alpha0, double w0)
    {
      lambda_p_ = 0;
      int p = sigmasq.length();
      for (int i = 0; i < p; i++)
      {
        lambda_p_ += 1 / sigmasq[i];
      }
      lambda_p_ *= alpha1 / 2.0;

      alpha_p_ = (p * alpha1 - alpha0) / 2.0;
      lambda0_ = alpha0 * w0 / 2.0;

      if (alpha_p_ < 1.0)
      {
        Rcpp::stop(
            "Error in 'R_sample_post_ichi:\n'"
            "Posterior alpha is less than 1, not log-concave\n");
      }
    }

    void eval_logf(const double x, double &logf, double &dlogf) override
    {
      double exps = exp(x);
      double iexps = 1 / exps;
      logf = alpha_p_ * x - lambda_p_ * exp(x) - lambda0_ * iexps;
      dlogf = alpha_p_ - lambda_p_ * exps + lambda0_ * iexps;
    }
  };

  // Do adaptive rejection sampling here
  auto target = IchiTarget(sigmasq, alpha1, alpha0, w0);
  auto sampler = ARS(n, &target, 0, R_NegInf, R_PosInf, verbose);
  return sampler.Sample(); 
}

