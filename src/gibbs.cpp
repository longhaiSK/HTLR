#include <string>
#include "gibbs.h"

Fit::Fit(int p, int K, int n, const arma::mat &X, const arma::mat &ymat, const arma::uvec &ybase,
         std::string ptype, double alpha, double s, double eta,
         int iters_rmc, int iters_h, int thin,
         int leap_L, int leap_L_h, double leap_step,
         double hmc_sgmcut, const arma::mat &deltas, const arma::vec &sigmasbt,
         bool keep_warmup_hist, int silence, bool legacy)
    : p_(p), K_(K), C_(K + 1), n_(n), X_(X), ymat_(ymat), ybase_(ybase),
      ptype_(ptype), alpha_(alpha), s_(s), eta_(eta),
      iters_rmc_(iters_rmc), iters_h_(iters_h), thin_(thin),
      leap_L_(leap_L), leap_L_h_(leap_L_h), leap_step_(leap_step),
      sgmsq_cut_(hmc_sgmcut > 0 ? R_pow_di(hmc_sgmcut, 2) : hmc_sgmcut),
      DDNloglike_(col_sum(arma::square(X)) / 4), keep_warmup_hist_(keep_warmup_hist),
      silence_(silence), legacy_(legacy), nvar_(p + 1), logw_(s), sigmasbt_(sigmasbt)
{
  ids_update_ = arma::uvec(nvar_, arma::fill::zeros);
  ids_fix_ = arma::uvec(nvar_, arma::fill::zeros);

  int hist_len = keep_warmup_hist_ ? (iters_rmc_ + iters_h_ + 1) : (iters_rmc_ + 1);

  mc_logw_ = arma::vec(hist_len, arma::fill::zeros);
  mc_logw_[0] = logw_;
  
  mc_sigmasbt_ = arma::mat(nvar_, hist_len, arma::fill::zeros);
  mc_sigmasbt_.col(0) = sigmasbt_;

  deltas_ = deltas;
  mc_deltas_ = arma::cube(nvar_, K, hist_len, arma::fill::zeros);
  mc_deltas_.slice(0) = deltas;

  mc_var_deltas_ = arma::mat(nvar_, hist_len, arma::fill::zeros);
  mc_loglike_ = arma::vec(hist_len, arma::fill::zeros);
  mc_uvar_ = arma::vec(hist_len, arma::fill::zeros);
  mc_hmcrej_ = arma::vec(hist_len, arma::fill::zeros);

  lv_ = arma::mat(n, C_, arma::fill::zeros);
  lv_fix_ = arma::mat(n, C_, arma::fill::zeros);

  DNloglike_ = arma::mat(nvar_, K, arma::fill::zeros);
  momt_ = arma::mat(nvar_, K, arma::fill::zeros);
  DNlogprior_ = arma::mat(nvar_, K, arma::fill::zeros);
  DNlogpost_ = arma::mat(nvar_, K, arma::fill::zeros);

  sumsq_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  sum_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  var_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  step_sizes_ = arma::vec(nvar_, arma::fill::zeros);
}

void Fit::StartSampling()
{
  Initialize();

  /************************ start gibbs sampling **************************/
  for (int i_mc = 0; i_mc < iters_h_ + iters_rmc_; i_mc++)
  {
    /***************** thin iterations of Gibbs sampling ******************/
    double no_uvar = 0;
    double rej = 0;
    for (int i_thin = 0; i_thin < thin_; i_thin++)
    {
      /*********************** HMC Metropolis Update ********************/
      
      // initialize HMC
      WhichUpdate();
      no_uvar += nuvar_;

      GenMomt();
      UpdateStepSizes();

      DetachFixlv();
      CacheOldValues();
      
      double nenergy_old = CompNegEnergy();

      // start trajectory
      UpdateDNlogPrior(); 
      UpdateDNlogLike(); 
      UpdateDNlogPost(); 
      Traject(i_mc);

      // decide whether to accept it
      UpdateLogLike();
      UpdateVarDeltas();
      double nenergy = CompNegEnergy();

      GetRNGstate();
      if (log(R::runif(0, 1)) > (nenergy - nenergy_old) || IsFault())
      {
        RestoreOldValues();
        rej++;
      }
      PutRNGstate();

      UpdateSigmas();
    }

    no_uvar /= thin_;
    rej /= thin_;

    /****************** record the markov chain state ********************/
    int i_rmc = keep_warmup_hist_ ? (i_mc + 1) : (i_mc - iters_h_ + 1);
    if (i_rmc > 0)
    {
      mc_deltas_.slice(i_rmc) = deltas_;
      mc_sigmasbt_.col(i_rmc) = sigmasbt_;
      mc_var_deltas_.col(i_rmc) = var_deltas_;
      mc_logw_[i_rmc] = logw_;
      mc_loglike_[i_rmc] = loglike_;
      mc_uvar_[i_rmc] = no_uvar;
      mc_hmcrej_[i_rmc] = rej;
    }

    // print some results on screen
    if (silence_ == 0)
    {
      Rprintf(
          "Iter%4d: deviance=%5.3f, logw=%6.2f, nuvar=%3.0f, hmcrej=%4.2f\n",
          i_mc - iters_h_, -loglike_ / n_, logw_, no_uvar, rej);
    }

    if (i_mc % 256 == 0) R_CheckUserInterrupt();
  }
}

Rcpp::List Fit::OutputR()
{
  auto mc_param = Rcpp::List::create(
      Rcpp::Named("iter.rmc") = iters_rmc_,
      Rcpp::Named("iter.warm") = iters_h_,
      Rcpp::Named("thin") = thin_,
      Rcpp::Named("leap") = leap_L_,
      Rcpp::Named("leap.warm") = leap_L_h_,
      Rcpp::Named("leap.step") = leap_step_,
      Rcpp::Named("sgmsq.cut") = sgmsq_cut_,
      Rcpp::Named("DDNloglike") = DDNloglike_);

  return Rcpp::List::create(
      Rcpp::Named("p") = p_,
      Rcpp::Named("n") = n_,
      Rcpp::Named("K") = K_,
      Rcpp::Named("mc.param") = mc_param,
      Rcpp::Named("mcdeltas") = mc_deltas_,
      Rcpp::Named("mclogw") = mc_logw_,
      Rcpp::Named("mcsigmasbt") = mc_sigmasbt_,
      Rcpp::Named("mcvardeltas") = mc_var_deltas_,
      Rcpp::Named("mcloglike") = mc_loglike_,
      Rcpp::Named("mcuvar") = mc_uvar_,
      Rcpp::Named("mchmcrej") = mc_hmcrej_);
}

// This function determines which features to be updated.
// Modified: nuvar, nfvar, ids_update, ids_fix_  
void Fit::WhichUpdate(bool init)
{
  nuvar_ = 0;
  nfvar_ = 0;
  double cut = init ? -1 : sgmsq_cut_;

  for (int j = 0; j < nvar_; j++)
  {
    if (sigmasbt_(j) > cut)
      ids_update_(nuvar_++) = j;
    else
      ids_fix_(nfvar_++) = j;
  }
  iup_ = ids_update_.head(nuvar_); // save a quick reference 
}

// X: n * nvar
// deltas: nvar * K
// lv: n * (1 + K)
// Modified: lv, norm_lv, pred_prob
void Fit::UpdatePredProb()
{
  lv_.tail_cols(K_) = lv_fix_.tail_cols(K_);
  for (int j : iup_)
  {
    for (int k = 0; k < K_; k++)
    {
      for (int i = 0; i < n_; i++)
      {
        lv_(i,k + 1) += X_(i, j) * deltas_(j, k);
      }
    }
  }
  norm_lv_ = find_normlv(lv_);
  pred_prob_ = arma::exp(norm_lv_);
}

// lv: n * (1 + K)
// deltas: nvar * K
// X: n * nvar
// Modified: lv_fix
void Fit::DetachFixlv()
{
  if (nuvar_ <= nvar_ / 2)
  {
    lv_fix_.tail_cols(K_) = lv_.tail_cols(K_); 
    // remove updated part
    for (int j : iup_)
    {
      for (int k = 0; k < K_; k++)
      {
        for (int i = 0; i < n_; i++)
        {
          lv_fix_(i, k + 1) -= X_(i, j) * deltas_(j, k);
        }
      }
    }
  }
  else
  {
    lv_fix_.tail_cols(K_) = arma::mat(n_, K_, arma::fill::zeros);
    // add fixed part
    for (int j : GetIdsFix())
    {
      for (int k = 0; k < K_; k++)
      {
        for (int i = 0; i < n_; i++)
        {
          lv_fix_(i, k + 1) += X_(i, j) * deltas_(j, k);
        }
      }
    }
  }
}

// DNloglike: nvar * K
// X: n * nvar
// pred_prob: n * (1 + K)
// ymat: n * K
// Modified: DNloglike
void Fit::UpdateDNlogLike()
{
  arma::mat tmp = pred_prob_.tail_cols(K_) - ymat_;
  for (int j : iup_)
  {
    for (int k = 0; k < K_; k++)
    {
      DNloglike_(j, k) = 0;
      for (int i = 0; i < n_; i++)
      {
        DNloglike_(j, k) += X_(i, j) * tmp(i, k);
      }
    }
  }
}

// norm_lv: n * C
// Modified: loglike
void Fit::UpdateLogLike()
{
  loglike_ = 0;
  for (int i = 0; i < n_; i++)
  {
    loglike_ += norm_lv_(i, ybase_(i));
  }
}

// deltas: nvar * K
// DNlogprior: nvar * K
// sum_deltas: nvar
// Modified: sum_deltas, DNlogprior:  
void Fit::UpdateDNlogPrior()
{
  arma::mat deltas_tmp = deltas_.rows(iup_);
  sum_deltas_(iup_) = row_sum(deltas_tmp); 
  DNlogprior_.rows(iup_) = deltas_tmp.each_col() - sum_deltas_(iup_) / C_;
}

// DNloglike: nvar * K
// DNlogprior: nvar * K
// DNlogpost: nvar * K
// sigmasbt: nvar
// Modified: DNlogpost 
void Fit::UpdateDNlogPost()
{
  arma::mat DNlogprior_tmp = DNlogprior_.rows(iup_);
  DNlogpost_.rows(iup_) = DNloglike_.rows(iup_) + DNlogprior_tmp.each_col() / sigmasbt_(iup_);
}

// This function is called at the beginning of the trajectory loop.
// momt: nvar * K
// step_sizes: nvar
// DNlogpost: nvar * K
// deltas: nvar * K
// Modified: momt, deltas
void Fit::UpdateMomtAndDeltas()
{
  arma::mat DNlogpost_tmp = DNlogpost_.rows(iup_);
  momt_.rows(iup_) -= step_sizes_(iup_) / 2 % DNlogpost_tmp.each_col();
  arma::mat momt_tmp = momt_.rows(iup_);
  deltas_.rows(iup_) += step_sizes_(iup_) % momt_tmp.each_col();
}

void Fit::UpdateSigmas()
{
  if (ptype_.compare("t") == 0)
    UpdateSigmasT();
  else if (ptype_.compare("ghs") == 0)
    UpdateSigmasGhs();
  else if (ptype_.compare("neg") == 0)
    UpdateSigmasNeg();
  else
    Rcpp::stop("Unsupported prior type %s", ptype_);
}

void Fit::UpdateSigmasT()
{
  double alpha_post = (alpha_ + K_) / 2;
  if (legacy_)
  {
    for (int j = 1; j < nvar_; j++)
    {
      GetRNGstate();
      sigmasbt_(j) =
          1.0 / R::rgamma(alpha_post, 1.0) * (alpha_ * exp(logw_) + var_deltas_[j]) / 2.0;
      PutRNGstate();
    }
  }
  else
  {
    arma::vec var_deltas_p = var_deltas_.tail(p_); 
    sigmasbt_.tail(p_) = spl_sgm_ig(alpha_, K_, exp(logw_), var_deltas_p);
  }

  UpdateLogw();
}

void Fit::UpdateLogw()
{
  // logw Update
  if (eta_ > 1E-10)
  {
    if (eta_ < 0.01)
      logw_ = s_;
    else
    {
      arma::vec tmp = var_deltas_.tail(p_);
      auto target = SamplerLogw(p_, tmp, K_, alpha_, s_, eta_);
      auto spl = ARS(1, &target, logw_);
      logw_ = spl.Sample()[0];
    }
  }
}

// Helper function for UpdateSigmasGhs and UpdateSigmasNeg  
void Fit::UpdateSigmasSgm(SamplerSgm *target)
{
  if (legacy_)
  {
    for (int j = 1; j < nvar_; j++)
    {
      target->set_idx(j);    
      auto spl = ARS(1, target, log(var_deltas_(j) / K_));
      sigmasbt_(j) = exp(spl.Sample()[0]); // perform ars on log(sigma_j)
    }
  }
  else
  {
    arma::vec tmp = arma::linspace(1, p_, p_);
    tmp.for_each([this, &target](arma::vec::elem_type &val) {
      target->set_idx(val);
      auto spl = ARS(1, target, log(var_deltas_(val) / K_));
      val = exp(spl.Sample()[0]); // perform ars on log(sigma_j)
    });
    sigmasbt_.tail(p_) = tmp;
  }
}

void Fit::UpdateSigmasGhs()
{
  auto *target = 
    new SamplerSgmGhs(nvar_, var_deltas_, K_, alpha_, logw_ + log(alpha_));
  UpdateSigmasSgm(target);
  delete target;
}

void Fit::UpdateSigmasNeg()
{
  auto *target = 
    new SamplerSgmNeg(nvar_, var_deltas_, K_, alpha_, logw_ + log(alpha_));
  UpdateSigmasSgm(target);
  delete target;
}

void Fit::Traject(int i_mc)
{
  int L;

  if (i_mc < iters_h_ / 2.0)
  {
    L = leap_L_h_;
    logw_ = -10;
  }
  else if (i_mc < iters_h_)
  {
    L = leap_L_h_;
    logw_ = s_;
  }
  else
  {
    L = leap_L_;
    logw_ = s_;
  }

  for (int i_trj = 0; i_trj < L; i_trj++)
  {
    UpdateMomtAndDeltas();
    UpdatePredProb();
    UpdateDNlogPrior();
    UpdateDNlogLike();
    UpdateDNlogPost();
    MoveMomt();
  }
}

// deltas: nvar * K
// sumsq_deltas: nvar
// var_deltas: nvar
// Modified: sumsq_deltas, var_deltas  
void Fit::UpdateVarDeltas()
{
  sumsq_deltas_(iup_) = row_sum(arma::square(deltas_.rows(iup_)));
  var_deltas_(iup_) = sumsq_deltas_(iup_) - arma::square(sum_deltas_(iup_)) / C_;
}

// momt: nvar * K
double Fit::CompNegEnergy()
{
  double logprior = arma::sum(var_deltas_(iup_) / sigmasbt_(iup_));
  double logprior_momt = arma::accu(arma::square(momt_.rows(iup_)));
  return (loglike_ - logprior / 2 - logprior_momt / 2);
}

// momt: nvar * K
// Modified: momt
void Fit::GenMomt()
{
  if (true)
  {
    for (int j : iup_)
    {
      for (int k = 0; k < K_; k++)
      {
        GetRNGstate();
        momt_(j, k) = R::rnorm(0, 1);
        PutRNGstate();
      }
    }
  }
  else // might have problem
  {
    arma::vec rn = Rcpp::rnorm(nuvar_ * K_);
    momt_.rows(iup_) = arma::reshape(rn, nuvar_, K_);
  }
}

// This function moves momonton with new derivatives.
// step_sizes: nvar
// DNlogpost: nvar * K
// momt: nvar * K
// Modified: momt
void Fit::MoveMomt()
{
  arma::mat DNlogpost_tmp = DNlogpost_.rows(iup_);
  momt_.rows(iup_) -= step_sizes_(iup_) / 2 % DNlogpost_tmp.each_col();
}

// step_sizes: nvar
// DDNloglike_: nvar
// sigmasbt: nvar
// Modified: step_sizes
void Fit::UpdateStepSizes()
{
  step_sizes_(iup_) =
      leap_step_ / arma::sqrt(DDNloglike_(iup_) + K_ / sigmasbt_(iup_) / C_);
}

void Fit::CacheOldValues()
{
  lv_old_ = lv_;
  pred_prob_old_ = pred_prob_;
  deltas_old_ = deltas_;
  DNlogprior_old_ = DNlogprior_;
  var_deltas_old_ = var_deltas_;
  loglike_old_ = loglike_;
}

void Fit::RestoreOldValues()
{
  lv_ = lv_old_;
  pred_prob_ = pred_prob_old_;
  deltas_ = deltas_old_;
  DNlogprior_ = DNlogprior_old_;
  var_deltas_ = var_deltas_old_;
  loglike_ = loglike_old_;
}

bool Fit::IsFault(double cri)
{
  for (int j : iup_)
  {
    for (int k = 0; k < K_; k++)
    {
      if (fabs(deltas_(j, k)) > cri)
      {
        return true;
      }
    }
  }
  return false;
}

// This function is called once at the beginning of the sampling process.
void Fit::Initialize()
{
  WhichUpdate(true); // set to update all
  UpdatePredProb();  // lv is computed here

  UpdateLogLike();
  mc_loglike_[0] = loglike_;

  UpdateDNlogPrior();
  UpdateVarDeltas(); 
  mc_var_deltas_.col(0) = var_deltas_;
}
