#include <string>
#include "gibbs.h"

Fit::Fit(int p, int K, int n,
         arma::mat &X, arma::mat &ymat, arma::vec &ybase,
         std::string ptype, double alpha, double s, double eta, double sigmab0,
         int iters_rmc, int iters_h, int thin,
         int leap_L, int leap_L_h, double leap_step,
         double hmc_sgmcut, arma::vec &DDNloglike,
         arma::cube &mcdeltas, arma::vec &mclogw, arma::mat &mcsigmasbt,
         int silence, int looklf)
    : p(p), K(K), C(K + 1), n(n), X(X.t()), ymat(ymat), ybase(ybase),
      ptype(ptype), alpha(alpha), s(s), eta(eta), sigmab0(sigmab0),
      iters_rmc(iters_rmc), iters_h(iters_h), thin(thin),
      leap_L(leap_L), leap_L_h(leap_L_h), leap_step(leap_step),
      hmc_sgmcut(hmc_sgmcut), DDNloglike(DDNloglike),
      silence(silence), looklf(looklf)
{
  this->nvar = p + 1;
  this->nuvar = 0;
  this->nfvar = 0;
  this->ids_update = arma::uvec(nvar, arma::fill::zeros); 
  this->ids_fix = arma::uvec(nvar, arma::fill::zeros);

  this->mcdeltas = mcdeltas;
  this->mclogw = mclogw;
  this->mcsigmasbt = mcsigmasbt;
  this->mcvardeltas = arma::mat(nvar, iters_rmc + 1, arma::fill::zeros);
  this->mcloglike = arma::vec(iters_rmc + 1, arma::fill::zeros);
  this->mcuvar = arma::vec(iters_rmc + 1, arma::fill::zeros);
  this->mchmcrej = arma::vec(iters_rmc + 1, arma::fill::zeros);

  this->lv = arma::mat(n, C, arma::fill::zeros);
  this->lv_old = arma::mat(n, C, arma::fill::zeros);
  this->lv_fix = arma::mat(n, C, arma::fill::zeros);
  this->normlv = arma::mat(n, C, arma::fill::zeros);
  this->predprob = arma::mat(n, C, arma::fill::zeros);
  this->predprob_old = arma::mat(n, C, arma::fill::zeros);

  this->DNloglike = arma::mat(K, nvar, arma::fill::zeros);
  this->DNloglike_old = arma::mat(K, nvar, arma::fill::zeros);
  this->deltas = arma::mat(K, nvar, arma::fill::zeros);
  this->deltas_old = arma::mat(K, nvar, arma::fill::zeros);
  this->momt = arma::mat(K, nvar, arma::fill::zeros);
  this->DNlogprior = arma::mat(K, nvar, arma::fill::zeros);
  this->DNlogprior_old = arma::mat(K, nvar, arma::fill::zeros);
  this->DNlogpost = arma::mat(K, nvar, arma::fill::zeros);

  this->SUMsqdeltas = arma::vec(nvar, arma::fill::zeros);
  this->SUMsqdeltas_old = arma::vec(nvar, arma::fill::zeros);
  this->SUMdeltas = arma::vec(nvar, arma::fill::zeros);
  this->SUMdeltas_old = arma::vec(nvar, arma::fill::zeros);
  this->vardeltas = arma::vec(nvar, arma::fill::zeros);
  this->vardeltas_old = arma::vec(nvar, arma::fill::zeros);
  this->stepsizes = arma::vec(nvar, arma::fill::zeros);
  this->sigmasbt = arma::vec(nvar, arma::fill::zeros);

  this->sgmsqcut = hmc_sgmcut > 0 ? R_pow_di(hmc_sgmcut, 2) : hmc_sgmcut;
}

void Fit::StartSampling()
{
  /*********************** getting initial values ************************/
  deltas = mcdeltas.slice(0);
  double logw = mclogw(0);
  sigmasbt = mcsigmasbt.col(0);

  // compute some initial values
  whichupdate(-1);  // set to update all
  updatepredprob(); // lv is computed here
  Rcpp::Rcerr << "here2.2\n";
  updateloglike();
  Rcpp::Rcerr << "here2.3\n";
  updateDNlogprior();
  Rcpp::Rcerr << "here2.4\n";
  updatevardeltas();
  Rcpp::Rcerr << "here2.5\n";
  mcloglike[0] = loglike;
  mcvardeltas.col(0) = vardeltas;
  Rcpp::Rcerr << "here3\n";

  /************************ start gibbs sampling **************************/
  for (int i_mc = 0; i_mc < iters_h + iters_rmc; i_mc++)
  {
    int i_rmc = i_mc;
    int L;

    if (i_mc < iters_h / 2.0)
    {
      L = leap_L_h;
      logw = -10;
    }
    else if (i_mc < iters_h)
    {
      L = leap_L_h;
      logw = s;
    }
    else
    {
      L = leap_L;
      logw = s;
    }

    /***************** thin iterations of Gibbs sampling ******************/
    double no_uvar = 0;
    double rej = 0;
    for (int i_thin = 0; i_thin < thin; i_thin++)
    {
      /*********************** HMC Metropolis Update ********************/
      // initialize HMC
      Rcpp::Rcerr << "here3.1\n";
      whichupdate(sgmsqcut);
      no_uvar += nuvar;
      Rcpp::Rcerr << "here3.2\n";

      gen_momt();
      updatestepsizes();
      Rcpp::Rcerr << "here3.3\n";
      detach_fixlv();
      Rcpp::Rcerr << "here3.4\n";
      cache_oldvalues();
      Rcpp::Rcerr << "here3.5\n";
      double nenergy_old = comp_nenergy();
      Rcpp::Rcerr << "here3.6\n";

      // start trajectory
      updateDNlogprior();
      Rcpp::Rcerr << "here3.7\n";
      updateDNloglike(); // recompute derivatives of log likelihood
      Rcpp::Rcerr << "here3.8\n";
      updateDNlogpost(); // recompute derivatives of log prior
      Rcpp::Rcerr << "here3.9\n";
      for (int i_trj = 0; i_trj < L; i_trj++)
      {
        for (int uj = 0; uj < nuvar; uj++)
        {
          int j = ids_update[uj];
          for (int k = 0; k < K; k++)
          {
            momt(k, j) -= stepsizes(j) / 2 * DNlogpost(k, j);
            deltas(k, j) += stepsizes(j) * momt(k, j);
          }
        }
        Rcpp::Rcerr << "here4\n";
        // compute derivative of minus log joint distribution
        updatepredprob();
        Rcpp::Rcerr << "here4.1\n";
        updateDNlogprior();
        Rcpp::Rcerr << "here4.2\n";
        updateDNloglike();
        Rcpp::Rcerr << "here4.3\n";
        updateDNlogpost();
        Rcpp::Rcerr << "here4.4\n";
        // move momonton with new derivatives
        for (int uj = 0; uj < nuvar; uj++)
        {
          int j = ids_update[uj];
          for (int k = 0; k < K; k++)
          {
            momt(k, j) -= stepsizes(j) / 2 * DNlogpost(k, j);
          }
        }
        Rcpp::Rcerr << "here4.5\n";
      }

      bool isfault = false;
      for (int uj = 0; uj < nuvar; uj++)
      {
        int j = ids_update[uj];
        for (int k = 0; k < K; k++)
        {
          if (fabs(deltas(k, j)) > 20)
          {
            isfault = true;
            break;
          }
        }
        if (isfault)
          break;
      }

      // decide whether to accept it
      updateloglike();
      updatevardeltas();
      double nenergy = comp_nenergy();

      GetRNGstate();
      if (log(R::runif(0, 1)) > nenergy - nenergy_old || isfault)
      {
        restore_oldvalues();
        rej++;
      }
      PutRNGstate();

      /*********************** Sigmas Update  ***********************/
      double log_aw = logw + log(alpha);

      if (ptype.compare("t") == 0)
      {
        for (int j = 1; j < nvar; j++)
        {
          GetRNGstate();
          double alpha_post = (alpha + K) / 2;
          sigmasbt[j] =
              1.0 / R::rgamma(alpha_post, 1.0) * (alpha * exp(logw) + vardeltas(j)) / 2.0;
          PutRNGstate();
        }
        /********************** logw Update  **********************/
        if (eta > 1E-10)
        {
          if (eta < 0.01)
          {
            logw = s;
          }
          else
          {
            auto target = SamplerLogw(p, vardeltas, K, alpha, s, eta);
            auto spl = ARS(1, &target, logw);
            logw = spl.Sample()[0];
          }
        }
      }

      if (ptype.compare("ghs") == 0)
      {
        //GetRNGstate();
        auto target = SamplerSgmGhs(p, vardeltas, K, alpha, log_aw);
        for (int idx = 0; idx < p; idx++)
        {
          // performa ars on log(sigma_j), which is still saved in sigma_j
          //if (log(vardeltas[j]) > -20)
          target.set_idx(idx);
          auto spl = ARS(1, &target, log(vardeltas[idx] / K));
          //convert xi to sigma_j
          sigmasbt[idx] = exp(spl.Sample()[0]);
        }
        //PutRNGstate();
      }

      if (ptype.compare("neg") == 0)
      {
        //GetRNGstate();
        auto target = SamplerSgmNeg(p, vardeltas, K, alpha, log_aw);
        for (int idx = 0; idx < p; idx++)
        {
          // performa ars on log(sigma_j), which is still saved in sigma_j
          //if (log(vardeltas[j]) > -20)
          target.set_idx(idx);
          ARS spl = ARS(1, &target, log(vardeltas[idx] / K));
          //convert xi to sigma_j
          sigmasbt[idx] = exp(spl.Sample()[0]);
        }
        //PutRNGstate();
      }
    }

    no_uvar /= thin;
    rej /= thin;

    /****************** record the markov chain state ********************/
    if (i_rmc >= 0)
    {
      mcdeltas.slice(i_rmc + 1) = deltas;
      mcsigmasbt.col(i_rmc + 1) = sigmasbt;
      mcvardeltas.col(i_rmc + 1) = vardeltas;
      mclogw(i_rmc + 1) = logw;
      mcloglike(i_rmc + 1) = loglike;
      mcuvar(i_rmc + 1) = no_uvar;
      mchmcrej(i_rmc + 1) = rej;
    }

    // print some results on screen
    if (silence == 0)
    {
      Rprintf(
          "Iter%4d: deviance=%5.3f, logw=%6.2f, nuvar=%3.0f, hmcrej=%4.2f\n",
          i_rmc, -loglike / n, logw, no_uvar, rej);
    }
  }
}

Rcpp::List Fit::OutputR()
{
    auto prior_param = Rcpp::List::create
    (
      Rcpp::Named("ptype") = ptype,
      Rcpp::Named("alpha") = alpha,
      Rcpp::Named("s") = s,
      Rcpp::Named("eta") = eta,
      Rcpp::Named("sigmab0") = sigmab0
    ); 

    auto mc_param = Rcpp::List::create
    (
      Rcpp::Named("iters.rmc") = iters_rmc,
      Rcpp::Named("iters.h") = iters_h,
      Rcpp::Named("thin") = thin,
      Rcpp::Named("leap.L") = leap_L,
      Rcpp::Named("leap.L.h") = leap_L_h,
      Rcpp::Named("leap.step") = leap_step,
      Rcpp::Named("hmc.sgmcut") = hmc_sgmcut,
      Rcpp::Named("DDNloglike") = DDNloglike
    );

    return Rcpp::List::create
    (
      Rcpp::Named("p") = p, 
      Rcpp::Named("n") = n,
      Rcpp::Named("K") = K,
      Rcpp::Named("prior") = prior_param,
      Rcpp::Named("mc.param") = mc_param,
      Rcpp::Named("mcdeltas") = mcdeltas,
      Rcpp::Named("mclogw") = mclogw,
      Rcpp::Named("mcsigmasbt") = mcsigmasbt,
      Rcpp::Named("mcvardeltas") = mcvardeltas,
      Rcpp::Named("mcloglike") = mcloglike,
      Rcpp::Named("mcuvar") = mcuvar,
      Rcpp::Named("mchmcrej") = mchmcrej
    );
  }

// this function determines which features to be updated 
void Fit::whichupdate(double cut)
{
  nuvar = 0;
  nfvar = 0;
  for (int j = 0; j < nvar; j++) 
  {
    if (sigmasbt[j] > cut) 
      ids_update[nuvar++] = j;
    else 
      ids_fix[nfvar++] = j;
  }
}

// this function updates lv and predprob
void Fit::updatepredprob()
{
  arma::mat tmp1 = deltas.cols(ids_update);
  arma::mat tmp2 = X.rows(ids_update);
  arma::mat tmp3 = tmp1 * tmp2; // updated part
  lv.tail_cols(K) = lv_fix.tail_cols(K) + tmp3.t();
  normlv = find_normlv(lv);
  predprob = arma::exp(normlv);
}

void Fit::detach_fixlv()
{
  if (nuvar <= nvar / 2)
  {
    lv_fix = cpvec(lv);
    // remove updated part
    for (int uj = 0; uj < nuvar; uj++)
    {
      int j = ids_update[uj];
      for (int k = 0; k < K; k++)
      {
        for (int i = 0; i < n; i++)
        {
          lv_fix(i, k + 1) -= deltas(k, j) * X(j, i);
        }
      }
    }
  }
  else
  {
    lv_fix.tail_cols(K) = arma::mat(n, K, arma::fill::zeros);
    // add fixed part
    for (int fj = 0; fj < nfvar; fj++)
    {
      int j = ids_fix[fj];
      for (int k = 0; k < K; k++)
      {
        for (int i = 0; i < n; i++)
        {
          lv_fix(i, k + 1) += deltas(k, j) * X(j, i);
        }
      }
    }
  }
}

// this function computes derivatives of negative log like
void Fit::updateDNloglike()
{
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    for (int k = 0; k < K; k++)
    {
      DNloglike(k, j) = 0; 
      for (int i = 0; i < n; i++) 
      {
        DNloglike(k, j) += X(j, i) * (predprob(i, k + 1) - ymat(k, i));  
      }
    }
  }
}

// this function computes log likelihood given normlv
void Fit::updateloglike()
{
  loglike = 0;
  for (int i = 0; i < n; i++) 
  {
    loglike += normlv(i, ybase[i]);
  }
}
  
// this function update SUMdeltas and DNlogprior
void Fit::updateDNlogprior()
{
  //arma::mat tmp1 = deltas.cols(ids_update);
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    //sum deltas for each feature
    SUMdeltas[j] = 0;
    for (int k = 0; k < K; k++) 
    {
      SUMdeltas[j] += deltas(k, j);
    //}
    //for (int k = 0; k < K; k++)
    //{    
      DNlogprior(k, j) = deltas(k, j) - SUMdeltas[j]/C;
    }
  }
}

void Fit::updateDNlogpost()
{       
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    for (int k = 0; k < K; k++)
    {    
      DNlogpost(k, j) = DNloglike(k, j) + DNlogprior(k, j) / sigmasbt[j];
    }
  }
}

// this function updates SUMsqdeltas, vardeltas
void Fit::updatevardeltas()
{
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    // prior deltas
    SUMsqdeltas[j] = 0;
    for (int k = 0; k < K; k++) 
    {
      SUMsqdeltas[j] += R_pow_di(deltas(k, j),2);
    }
    vardeltas[j] = SUMsqdeltas[j] - R_pow_di(SUMdeltas[j],2) / C;
  }
}

// this function compute neg. energy
double Fit::comp_nenergy()
{
  double logprior = 0;
  double logprior_momt = 0;
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    logprior += vardeltas(j) / sigmasbt[j];
    for (int k = 0; k < K; k++)  
      logprior_momt += R_pow_di(momt(k, j),2);
  }
  return (loglike- logprior/2 - logprior_momt/2);
}
  
void Fit::gen_momt()
{
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    for (int k = 0; k < K; k++)
    {
      GetRNGstate();
      momt(k, j) = R::rnorm(0,1);
      PutRNGstate();
    }
  } 
}
  
void Fit::cache_oldvalues()
{
  lv_old = cpvec(lv);
  predprob_old = cpvec(predprob);
  loglike_old = loglike;

  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    for (int k = 0; k < K; k++)
    {
      deltas_old(k, j) = deltas(k, j);
      DNlogprior_old(k, j) = DNlogprior(k, j);
    }
    vardeltas_old(j) = vardeltas(j);
  }
}

void Fit::restore_oldvalues()
{
  lv = cpvec(lv_old);
  predprob = cpvec(predprob_old);
  loglike = loglike_old;

  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    for (int k = 0; k < K; k++)
    {
        deltas(k, j) = deltas_old(k, j);
        DNlogprior(k, j) = DNlogprior_old(k, j);
    }
    vardeltas(j) = vardeltas_old(j);
  }
}

void Fit::updatestepsizes()
{
  for (int uj = 0; uj < nuvar; uj++) 
  {
    int j = ids_update[uj];
    stepsizes(j) = leap_step;
    stepsizes(j) /= sqrt(DDNloglike(j) + K / sigmasbt[j] / C);
  }
}
