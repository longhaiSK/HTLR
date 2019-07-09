#include <string>
#include "gibbs.h"

Fit::Fit(int p, int K, int n,
         arma::mat &X, arma::mat &ymat, arma::uvec &ybase,
         std::string ptype, double alpha, double s, double eta, double sigmab0,
         int iters_rmc, int iters_h, int thin,
         int leap_L, int leap_L_h, double leap_step,
         double hmc_sgmcut, arma::vec &DDNloglike_,
         arma::mat &deltas, double logw, arma::vec &sigmasbt,
         int silence, int looklf)
    : p_(p), K_(K), C_(K + 1), n_(n), X_(X), ymat_(ymat), ybase_(ybase),
      ptype_(ptype), alpha_(alpha), s_(s), eta_(eta), sigmab0_(sigmab0),
      iters_rmc_(iters_rmc), iters_h_(iters_h), thin_(thin),
      leap_L_(leap_L), leap_L_h_(leap_L_h), leap_step_(leap_step),
      DDNloglike_(DDNloglike_), silence_(silence), looklf_(looklf), nvar_(p + 1), logw_(logw)
{
  ids_update_ = arma::uvec(nvar_, arma::fill::zeros);
  ids_fix_ = arma::uvec(nvar_, arma::fill::zeros);

  mc_logw_ = arma::vec(iters_rmc + 1, arma::fill::zeros);
  mc_logw_[0] = logw;
  
  sigmasbt_ = sigmasbt;
  mc_sigmasbt_ = arma::mat(nvar_, iters_rmc + 1, arma::fill::zeros);
  mc_sigmasbt_.col(0) = sigmasbt;

  deltas_ = deltas;
  mc_deltas_ = arma::cube(nvar_, K, iters_rmc + 1, arma::fill::zeros);
  mc_deltas_.slice(0) = deltas;

  mc_var_deltas_ = arma::mat(nvar_, iters_rmc + 1, arma::fill::zeros);
  mc_loglike_ = arma::vec(iters_rmc + 1, arma::fill::zeros);
  mc_uvar_ = arma::vec(iters_rmc + 1, arma::fill::zeros);
  mc_hmcrej_ = arma::vec(iters_rmc + 1, arma::fill::zeros);

  lv_ = arma::mat(n, C_, arma::fill::zeros);
  lv_old_ = arma::mat(n, C_, arma::fill::zeros);
  lv_fix_ = arma::mat(n, C_, arma::fill::zeros);
  norm_lv_ = arma::mat(n, C_, arma::fill::zeros);
  pred_prob_ = arma::mat(n, C_, arma::fill::zeros);
  pred_prob_old_ = arma::mat(n, C_, arma::fill::zeros);

  DNloglike_ = arma::mat(nvar_, K, arma::fill::zeros);
  //DNloglike_old_ = arma::mat(nvar_, K, arma::fill::zeros);
  deltas_old_ = arma::mat(nvar_, K, arma::fill::zeros);
  momt_ = arma::mat(nvar_, K, arma::fill::zeros);
  DNlogprior_ = arma::mat(nvar_, K, arma::fill::zeros);
  DNlogprior_old_ = arma::mat(nvar_, K, arma::fill::zeros);
  DNlogpost_ = arma::mat(nvar_, K, arma::fill::zeros);

  sumsq_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  sumsq_deltas_old_ = arma::vec(nvar_, arma::fill::zeros);
  sum_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  sum_deltas_old_ = arma::vec(nvar_, arma::fill::zeros);
  var_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  var_deltas_old_ = arma::vec(nvar_, arma::fill::zeros);
  step_sizes_ = arma::vec(nvar_, arma::fill::zeros);

  sgmsq_cut_ = hmc_sgmcut > 0 ? R_pow_di(hmc_sgmcut, 2) : hmc_sgmcut;
}

void Fit::StartSampling()
{
  Initialize();

  /************************ start gibbs sampling **************************/
  for (int i_mc = 0; i_mc < iters_h_ + iters_rmc_; i_mc++)
  {
    int i_rmc = i_mc - iters_h_;
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

      for (int i_trj = 0; i_trj < L; i_trj++)
      {
        UpdateMomtAndDeltas();
        // compute derivative of minus log joint distribution
        UpdatePredProb();
        UpdateDNlogPrior();
        UpdateDNlogLike();
        UpdateDNlogPost();
        MoveMomt();
      }

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

      /*********************** Sigmas Update  ***********************/
      
      if (ptype_.compare("t") == 0)
      {
        double alpha_post = (alpha_ + K_) / 2;
        
        for (int j = 1; j < nvar_; j++)
        {
          GetRNGstate();
          // TODO: use for_each
          sigmasbt_[j] =
              1.0 / R::rgamma(alpha_post, 1.0) * (alpha_ * exp(logw_) + var_deltas_[j]) / 2.0;
          PutRNGstate();
        }
        /********************** logw Update  **********************/
        if (eta_ > 1E-10)
        {
          if (eta_ < 0.01)
          {
            logw_ = s_;
          }
          else
          {
            auto target = SamplerLogw(p_, var_deltas_, K_, alpha_, s_, eta_);
            auto spl = ARS(1, &target, logw_);
            logw_ = spl.Sample()[0];
          }
        }
      }

      double log_aw = logw_ + log(alpha_);

      if (ptype_.compare("ghs") == 0)
      {
        //GetRNGstate();
        auto target = SamplerSgmGhs(p_, var_deltas_, K_, alpha_, log_aw);
        for (int idx = 0; idx < p_; idx++)
        {
          // performa ars on log(sigma_j), which is still saved in sigma_j
          //if (log(var_deltas_[j]) > -20)
          target.set_idx(idx);
          auto spl = ARS(1, &target, log(var_deltas_[idx] / K_));
          //convert xi to sigma_j
          sigmasbt_[idx] = exp(spl.Sample()[0]);
        }
        //PutRNGstate();
      }

      if (ptype_.compare("neg") == 0)
      {
        //GetRNGstate();
        auto target = SamplerSgmNeg(p_, var_deltas_, K_, alpha_, log_aw);
        for (int idx = 0; idx < p_; idx++)
        {
          // performa ars on log(sigma_j), which is still saved in sigma_j
          //if (log(var_deltas_[j]) > -20)
          target.set_idx(idx);
          ARS spl = ARS(1, &target, log(var_deltas_[idx] / K_));
          //convert xi to sigma_j
          sigmasbt_[idx] = exp(spl.Sample()[0]);
        }
        //PutRNGstate();
      }
    }

    no_uvar /= thin_;
    rej /= thin_;

    /****************** record the markov chain state ********************/
    if (i_rmc >= 0)
    {
      mc_deltas_.slice(i_rmc + 1) = deltas_;
      mc_sigmasbt_.col(i_rmc + 1) = sigmasbt_;
      mc_var_deltas_.col(i_rmc + 1) = var_deltas_;
      mc_logw_(i_rmc + 1) = logw_;
      mc_loglike_(i_rmc + 1) = loglike_;
      mc_uvar_(i_rmc + 1) = no_uvar;
      mc_hmcrej_(i_rmc + 1) = rej;
    }

    // print some results on screen
    if (silence_ == 0)
    {
      Rprintf(
          "Iter%4d: deviance=%5.3f, logw=%6.2f, nuvar=%3.0f, hmcrej=%4.2f\n\n",
          i_rmc, -loglike_ / n_, logw_, no_uvar, rej);
    }
    if (i_rmc % 256 == 0) R_CheckUserInterrupt();
  }
}

Rcpp::List Fit::OutputR()
{
  auto prior_param = Rcpp::List::create(
      Rcpp::Named("ptype") = ptype_,
      Rcpp::Named("alpha") = alpha_,
      Rcpp::Named("s") = s_,
      Rcpp::Named("eta") = eta_,
      Rcpp::Named("sigmab0") = sigmab0_);

  auto mc_param = Rcpp::List::create(
      Rcpp::Named("iters.rmc") = iters_rmc_,
      Rcpp::Named("iters.h") = iters_h_,
      Rcpp::Named("thin") = thin_,
      Rcpp::Named("leap.L") = leap_L_,
      Rcpp::Named("leap.L.h") = leap_L_h_,
      Rcpp::Named("leap.step") = leap_step_,
      Rcpp::Named("sgmsq.cut") = sgmsq_cut_,
      Rcpp::Named("DDNloglike") = DDNloglike_);

  return Rcpp::List::create(
      Rcpp::Named("p") = p_,
      Rcpp::Named("n") = n_,
      Rcpp::Named("K") = K_,
      Rcpp::Named("prior") = prior_param,
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
}

// Get int vector ids_update of length nuvar. 
arma::uvec Fit::GetIdsUpdate()
{
  return ids_update_.head(nuvar_);
}

// Get int vector ids_fix of length nfvar.
arma::uvec Fit::GetIdsFix()
{
  return ids_fix_.head(nfvar_);
}

// X: n * nvar
// deltas: nvar * K
// lv: n * (1 + K)
// Modified: lv, norm_lv, pred_prob
void Fit::UpdatePredProb()
{
  lv_.tail_cols(K_) = lv_fix_.tail_cols(K_);
  for (int j : GetIdsUpdate())
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
    for (int j : GetIdsUpdate())
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
  for (int j : GetIdsUpdate())
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
  auto ids = GetIdsUpdate();
  sum_deltas_(ids) = row_sum(deltas_.rows(ids)); 
  DNlogprior_.rows(ids) = deltas_.rows(ids) - sum_deltas_(ids) / C_;
}

// DNloglike: nvar * K
// DNlogprior: nvar * K
// DNlogpost: nvar * K
// sigmasbt: nvar
// Modified: DNlogpost 
void Fit::UpdateDNlogPost()
{
  auto ids = GetIdsUpdate();
  DNlogpost_.rows(ids) = DNloglike_.rows(ids) + DNlogprior_.rows(ids) / sigmasbt_(ids);
}

// This function is called at the beginning of the trajectory loop.
// momt: nvar * K
// step_sizes: nvar
// DNlogpost: nvar * K
// deltas: nvar * K
// Modified: momt, deltas
void Fit::UpdateMomtAndDeltas()
{
  for (int j : GetIdsUpdate())
  {
    momt_.row(j) -= step_sizes_(j) / 2 * DNlogpost_.row(j);
    deltas_.row(j) += step_sizes_(j) * momt_.row(j);
  }
}

// deltas: nvar * K
// sumsq_deltas: nvar
// var_deltas: nvar
// Modified: sumsq_deltas, var_deltas  
void Fit::UpdateVarDeltas()
{
  auto ids = GetIdsUpdate();
  sumsq_deltas_(ids) = row_sum(arma::square(deltas_.rows(ids)));
  var_deltas_(ids) = sumsq_deltas_(ids) - arma::square(sum_deltas_(ids)) / C_;
}

// momt: nvar * K
double Fit::CompNegEnergy()
{
  double logprior = 0;
  double logprior_momt = 0;
  for (int j : GetIdsUpdate())
  {
    logprior += var_deltas_[j] / sigmasbt_[j];
    for (int k = 0; k < K_; k++)
      logprior_momt += R_pow_di(momt_(j, k), 2);
  }
  return (loglike_ - logprior / 2 - logprior_momt / 2);
}

// momt: nvar * K
// Modified: momt
void Fit::GenMomt()
{
  for (int j : GetIdsUpdate())
  {
    for (int k = 0; k < K_; k++)
    {
      GetRNGstate();
      momt_(j, k) = R::rnorm(0, 1);
      PutRNGstate();
    }
  }
}

// This function moves momonton with new derivatives.
// step_sizes: nvar
// DNlogpost: nvar * K
// momt: nvar * K
// Modified: momt
void Fit::MoveMomt()
{
  for (int j : GetIdsUpdate())
  {
    momt_.row(j) -= step_sizes_(j) / 2 * DNlogpost_.row(j); 
  }
}

// step_sizes: nvar
// DDNloglike_: nvar
// sigmasbt: nvar
// Modified: step_sizes
void Fit::UpdateStepSizes()
{
  auto ids = GetIdsUpdate();
  step_sizes_(ids) = leap_step_ /
                     arma::sqrt(DDNloglike_(ids) + K_ / sigmasbt_(ids) / C_);
}

void Fit::CacheOldValues()
{
  lv_old_ = copy(lv_);
  pred_prob_old_ = copy(pred_prob_);
  deltas_old_ = copy(deltas_);
  DNlogprior_old_ = copy(DNlogprior_);
  var_deltas_old_ = copy(var_deltas_);
  loglike_old_ = loglike_;
}

void Fit::RestoreOldValues()
{
  lv_ = copy(lv_old_);
  pred_prob_ = copy(pred_prob_old_);
  deltas_ = copy(deltas_old_);
  DNlogprior_ = copy(DNlogprior_old_);
  var_deltas_ = copy(var_deltas_old_);
  loglike_ = loglike_old_;
}

bool Fit::IsFault(double cri)
{
  for (int j : GetIdsUpdate())
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