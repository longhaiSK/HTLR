#include <string>
#include "gibbs.h"

Fit::Fit(int p, int K, int n,
         arma::mat &X, arma::mat &ymat, arma::vec &ybase,
         std::string ptype, double alpha, double s, double eta, double sigmab0,
         int iters_rmc, int iters_h, int thin,
         int leap_L, int leap_L_h, double leap_step,
         double hmc_sgmcut, arma::vec &DDNloglike_,
         arma::cube &mcdeltas, arma::vec &mclogw, arma::mat &mcsigmasbt,
         int silence, int looklf)
    : p_(p), K_(K), C_(K + 1), n_(n), X_(X.t()), ymat_(ymat), ybase_(ybase),
      ptype_(ptype), alpha_(alpha), s_(s), eta_(eta), sigmab0_(sigmab0),
      iters_rmc_(iters_rmc), iters_h_(iters_h), thin_(thin),
      leap_L_(leap_L), leap_L_h_(leap_L_h), leap_step_(leap_step),
      hmc_sgmcut_(hmc_sgmcut), DDNloglike_(DDNloglike_),
      silence_(silence), looklf_(looklf)
{
  this->nvar_ = p + 1;
  this->nuvar_ = 0;
  this->nfvar_ = 0;
  this->ids_update_ = arma::uvec(nvar_, arma::fill::zeros); 
  this->ids_fix_ = arma::uvec(nvar_, arma::fill::zeros);

  this->mc_deltas_ = mcdeltas;
  this->mc_logw_ = mclogw;
  this->mc_sigmasbt_ = mcsigmasbt;
  this->mc_var_deltas_ = arma::mat(nvar_, iters_rmc + 1, arma::fill::zeros);
  this->mc_loglike_ = arma::vec(iters_rmc + 1, arma::fill::zeros);
  this->mc_uvar_ = arma::vec(iters_rmc + 1, arma::fill::zeros);
  this->mc_hmcrej_ = arma::vec(iters_rmc + 1, arma::fill::zeros);

  this->lv_ = arma::mat(n, C_, arma::fill::zeros);
  this->lv_old_ = arma::mat(n, C_, arma::fill::zeros);
  this->lv_fix_ = arma::mat(n, C_, arma::fill::zeros);
  this->norm_lv_ = arma::mat(n, C_, arma::fill::zeros);
  this->pred_prob_ = arma::mat(n, C_, arma::fill::zeros);
  this->pred_prob_old_ = arma::mat(n, C_, arma::fill::zeros);

  this->DNloglike_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->DNloglike_old_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->deltas_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->deltas_old_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->momt_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->DNlogprior_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->DNlogprior_old_ = arma::mat(K, nvar_, arma::fill::zeros);
  this->DNlogpost_ = arma::mat(K, nvar_, arma::fill::zeros);

  this->sumsq_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  this->sumsq_deltas_old_ = arma::vec(nvar_, arma::fill::zeros);
  this->sum_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  this->sum_deltas_old_ = arma::vec(nvar_, arma::fill::zeros);
  this->var_deltas_ = arma::vec(nvar_, arma::fill::zeros);
  this->var_deltas_old_ = arma::vec(nvar_, arma::fill::zeros);
  this->step_sizes_ = arma::vec(nvar_, arma::fill::zeros);
  this->sigmasbt_ = arma::vec(nvar_, arma::fill::zeros);

  this->sgmsq_cut_ = hmc_sgmcut > 0 ? R_pow_di(hmc_sgmcut, 2) : hmc_sgmcut;
}

void Fit::StartSampling()
{
  /*********************** getting initial values ************************/
  deltas_ = mc_deltas_.slice(0);
  double logw = mc_logw_(0);
  sigmasbt_ = mc_sigmasbt_.col(0);

  // compute some initial values
  WhichUpdate(-1);  // set to update all
  UpdatePredProb(); // lv is computed here
  Rcpp::Rcerr << "here2.2\n";
  UpdateLogLike();
  Rcpp::Rcerr << "here2.3\n";
  UpdateDNlogPrior();
  Rcpp::Rcerr << "here2.4\n";
  UpdateVarDeltas();
  Rcpp::Rcerr << "here2.5\n";
  mc_loglike_[0] = loglike_;
  mc_var_deltas_.col(0) = var_deltas_;
  Rcpp::Rcerr << "here3\n";

  /************************ start gibbs sampling **************************/
  for (int i_mc = 0; i_mc < iters_h_ + iters_rmc_; i_mc++)
  {
    int i_rmc = i_mc;
    int L;

    if (i_mc < iters_h_ / 2.0)
    {
      L = leap_L_h_;
      logw = -10;
    }
    else if (i_mc < iters_h_)
    {
      L = leap_L_h_;
      logw = s_;
    }
    else
    {
      L = leap_L_;
      logw = s_;
    }

    /***************** thin iterations of Gibbs sampling ******************/
    double no_uvar = 0;
    double rej = 0;
    for (int i_thin = 0; i_thin < thin_; i_thin++)
    {
      /*********************** HMC Metropolis Update ********************/
      // initialize HMC
      Rcpp::Rcerr << "here3.1\n";
      WhichUpdate(sgmsq_cut_);
      no_uvar += nuvar_;
      Rcpp::Rcerr << "here3.2\n";

      gen_momt();
      updatestepsizes();
      Rcpp::Rcerr << "here3.3\n";
      DetachFixlv();
      Rcpp::Rcerr << "here3.4\n";
      cache_oldvalues();
      Rcpp::Rcerr << "here3.5\n";
      double nenergy_old = comp_nenergy();
      Rcpp::Rcerr << "here3.6\n";

      // start trajectory
      UpdateDNlogPrior();
      Rcpp::Rcerr << "here3.7\n";
      UpdateDNlogLike(); // recompute derivatives of log likelihood
      Rcpp::Rcerr << "here3.8\n";
      UpdateDNlogPost(); // recompute derivatives of log prior
      Rcpp::Rcerr << "here3.9\n";
      for (int i_trj = 0; i_trj < L; i_trj++)
      {
        for (int uj = 0; uj < nuvar_; uj++)
        {
          int j = ids_update_[uj];
          for (int k = 0; k < K_; k++)
          {
            momt_(k, j) -= step_sizes_(j) / 2 * DNlogpost_(k, j);
            deltas_(k, j) += step_sizes_(j) * momt_(k, j);
          }
        }
        Rcpp::Rcerr << "here4\n";
        // compute derivative of minus log joint distribution
        UpdatePredProb();
        Rcpp::Rcerr << "here4.1\n";
        UpdateDNlogPrior();
        Rcpp::Rcerr << "here4.2\n";
        UpdateDNlogLike();
        Rcpp::Rcerr << "here4.3\n";
        UpdateDNlogPost();
        Rcpp::Rcerr << "here4.4\n";
        // move momonton with new derivatives
        for (int uj = 0; uj < nuvar_; uj++)
        {
          int j = ids_update_[uj];
          for (int k = 0; k < K_; k++)
          {
            momt_(k, j) -= step_sizes_(j) / 2 * DNlogpost_(k, j);
          }
        }
        Rcpp::Rcerr << "here4.5\n";
      }

      bool isfault = false;
      for (int uj = 0; uj < nuvar_; uj++)
      {
        int j = ids_update_[uj];
        for (int k = 0; k < K_; k++)
        {
          if (fabs(deltas_(k, j)) > 20)
          {
            isfault = true;
            break;
          }
        }
        if (isfault)
          break;
      }

      // decide whether to accept it
      UpdateLogLike();
      UpdateVarDeltas();
      double nenergy = comp_nenergy();

      GetRNGstate();
      if (log(R::runif(0, 1)) > nenergy - nenergy_old || isfault)
      {
        restore_oldvalues();
        rej++;
      }
      PutRNGstate();

      /*********************** Sigmas Update  ***********************/
      double log_aw = logw + log(alpha_);

      if (ptype_.compare("t") == 0)
      {
        for (int j = 1; j < nvar_; j++)
        {
          GetRNGstate();
          double alpha_post = (alpha_ + K_) / 2;
          sigmasbt_[j] =
              1.0 / R::rgamma(alpha_post, 1.0) * (alpha_ * exp(logw) + var_deltas_(j)) / 2.0;
          PutRNGstate();
        }
        /********************** logw Update  **********************/
        if (eta_ > 1E-10)
        {
          if (eta_ < 0.01)
          {
            logw = s_;
          }
          else
          {
            auto target = SamplerLogw(p_, var_deltas_, K_, alpha_, s_, eta_);
            auto spl = ARS(1, &target, logw);
            logw = spl.Sample()[0];
          }
        }
      }

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
      mc_logw_(i_rmc + 1) = logw;
      mc_loglike_(i_rmc + 1) = loglike_;
      mc_uvar_(i_rmc + 1) = no_uvar;
      mc_hmcrej_(i_rmc + 1) = rej;
    }

    // print some results on screen
    if (silence_ == 0)
    {
      Rprintf(
          "Iter%4d: deviance=%5.3f, logw=%6.2f, nuvar=%3.0f, hmcrej=%4.2f\n",
          i_rmc, -loglike_ / n_, logw, no_uvar, rej);
    }
  }
}

Rcpp::List Fit::OutputR()
{
    auto prior_param = Rcpp::List::create
    (
      Rcpp::Named("ptype") = ptype_,
      Rcpp::Named("alpha") = alpha_,
      Rcpp::Named("s") = s_,
      Rcpp::Named("eta") = eta_,
      Rcpp::Named("sigmab0") = sigmab0_
    ); 

    auto mc_param = Rcpp::List::create
    (
      Rcpp::Named("iters.rmc") = iters_rmc_,
      Rcpp::Named("iters.h") = iters_h_,
      Rcpp::Named("thin") = thin_,
      Rcpp::Named("leap.L") = leap_L_,
      Rcpp::Named("leap.L.h") = leap_L_h_,
      Rcpp::Named("leap.step") = leap_step_,
      Rcpp::Named("hmc.sgmcut") = hmc_sgmcut_,
      Rcpp::Named("DDNloglike_") = DDNloglike_
    );

    return Rcpp::List::create
    (
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
      Rcpp::Named("mchmcrej") = mc_hmcrej_
    );
  }

// this function determines which features to be updated 
void Fit::WhichUpdate(double cut)
{
  nuvar_ = 0;
  nfvar_ = 0;
  for (int j = 0; j < nvar_; j++) 
  {
    if (sigmasbt_[j] > cut) 
      ids_update_[nuvar_++] = j;
    else 
      ids_fix_[nfvar_++] = j;
  }
}

// This function updates lv and pred_prob_.
// deltas: 
void Fit::UpdatePredProb()
{
  arma::mat deltas_update = deltas_.cols(ids_update_);
  arma::mat X_update = X_.rows(ids_update_);
  arma::mat updated_part = deltas_update * X_update;
  lv_.tail_cols(K_) = lv_fix_.tail_cols(K_) + updated_part.t();
  norm_lv_ = find_normlv(lv_);
  pred_prob_ = arma::exp(norm_lv_);
}

void Fit::DetachFixlv()
{
  if (nuvar_ <= nvar_ / 2)
  {
    lv_fix_ = cpvec(lv_);
    // remove updated part
    for (int uj = 0; uj < nuvar_; uj++)
    {
      int j = ids_update_[uj];
      for (int k = 0; k < K_; k++)
      {
        for (int i = 0; i < n_; i++)
        {
          lv_fix_(i, k + 1) -= deltas_(k, j) * X_(j, i);
        }
      }
    }
  }
  else
  {
    lv_fix_.tail_cols(K_) = arma::mat(n_, K_, arma::fill::zeros);
    // add fixed part
    for (int fj = 0; fj < nfvar_; fj++)
    {
      int j = ids_fix_[fj];
      for (int k = 0; k < K_; k++)
      {
        for (int i = 0; i < n_; i++)
        {
          lv_fix_(i, k + 1) += deltas_(k, j) * X_(j, i);
        }
      }
    }
  }
}

// This function computes derivatives of negative log-likelihood. 
// 
void Fit::UpdateDNlogLike()
{
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    for (int k = 0; k < K_; k++)
    {
      DNloglike_(k, j) = 0; 
      for (int i = 0; i < n_; i++) 
      {
        DNloglike_(k, j) += X_(j, i) * (pred_prob_(i, k + 1) - ymat_(k, i));  
      }
    }
  }
}

// this function computes log likelihood given normlv
void Fit::UpdateLogLike()
{
  loglike_ = 0;
  for (int i = 0; i < n_; i++) 
  {
    loglike_ += norm_lv_(i, ybase_[i]);
  }
}
  
// this function update sum_deltas_ and DNlogprior
void Fit::UpdateDNlogPrior()
{
  //arma::mat tmp1 = deltas.cols(ids_update);
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    //sum deltas for each feature
    sum_deltas_[j] = 0;
    for (int k = 0; k < K_; k++) 
    {
      sum_deltas_[j] += deltas_(k, j);
    //}
    //for (int k = 0; k < K; k++)
    //{    
      DNlogprior_(k, j) = deltas_(k, j) - sum_deltas_[j] / C_;
    }
  }
}

void Fit::UpdateDNlogPost()
{       
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    for (int k = 0; k < K_; k++)
    {    
      DNlogpost_(k, j) = DNloglike_(k, j) + DNlogprior_(k, j) / sigmasbt_[j];
    }
  }
}

// this function updates sumsq_deltas_, var_deltas_
void Fit::UpdateVarDeltas()
{
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    // prior deltas
    sumsq_deltas_[j] = 0;
    for (int k = 0; k < K_; k++) 
    {
      sumsq_deltas_[j] += R_pow_di(deltas_(k, j),2);
    }
    var_deltas_[j] = sumsq_deltas_[j] - R_pow_di(sum_deltas_[j],2) / C_;
  }
}

// this function compute neg. energy
double Fit::comp_nenergy()
{
  double logprior = 0;
  double logprior_momt_ = 0;
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    logprior += var_deltas_(j) / sigmasbt_[j];
    for (int k = 0; k < K_; k++)  
      logprior_momt_ += R_pow_di(momt_(k, j),2);
  }
  return (loglike_- logprior/2 - logprior_momt_/2);
}
  
void Fit::gen_momt()
{
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    for (int k = 0; k < K_; k++)
    {
      GetRNGstate();
      momt_(k, j) = R::rnorm(0,1);
      PutRNGstate();
    }
  } 
}
  
void Fit::cache_oldvalues()
{
  lv_old_ = cpvec(lv_);
  pred_prob_old_ = cpvec(pred_prob_);
  loglike_old_ = loglike_;

  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    for (int k = 0; k < K_; k++)
    {
      deltas_old_(k, j) = deltas_(k, j);
      DNlogprior_old_(k, j) = DNlogprior_(k, j);
    }
    var_deltas_old_(j) = var_deltas_(j);
  }
}

void Fit::restore_oldvalues()
{
  lv_ = cpvec(lv_old_);
  pred_prob_ = cpvec(pred_prob_old_);
  loglike_ = loglike_old_;

  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    for (int k = 0; k < K_; k++)
    {
        deltas_(k, j) = deltas_old_(k, j);
        DNlogprior_(k, j) = DNlogprior_old_(k, j);
    }
    var_deltas_(j) = var_deltas_old_(j);
  }
}

void Fit::updatestepsizes()
{
  for (int uj = 0; uj < nuvar_; uj++) 
  {
    int j = ids_update_[uj];
    step_sizes_(j) = leap_step_;
    step_sizes_(j) /= sqrt(DDNloglike_(j) + K_ / sigmasbt_[j] / C_);
  }
}
