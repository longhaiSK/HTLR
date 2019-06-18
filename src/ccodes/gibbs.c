
# include <math.h>
# include <Rmath.h>
# include <stdio.h>
# include <R.h>
# include <string.h>

void testdiv (char *a[])
{
    printf ("%s prior\n", a[0]);
    printf ("%s prior\n", a[1]);
    printf ("%s prior\n", a[2]);
    
}

/********* sampling from posterior of sigma given var and prior ************/
void spl_sgm_neg (int p[1], double sigmasbt[p[0]], double vardeltas[p[0]], 
				  int K[1], double alpha[1], double log_aw [1])
{
  int j;
  
  //define the function computing logf and dlogf for ars
  void logpost_xi (double xi, double logf[1], double dlogf[1])
  {
	double vexi, eximinlogaw;
	
	logf[0]  = - (K[0] / 2.0 - 1.0) * xi;
	dlogf[0] = - (K[0] / 2.0 - 1.0);
	
	vexi = vardeltas [j] / 2.0 / exp (xi);
	logf[0]  -= vexi;
	dlogf[0] += vexi;
	
	eximinlogaw = exp(xi - log_aw[0]);
	logf[0] -= (alpha[0]/2.0 + 1.0) * log (1.0 + 2.0 * eximinlogaw);
	dlogf[0] -= (alpha[0] + 2.0) * eximinlogaw / (1.0 + 2.0 * eximinlogaw);
  }
  
  for (j = 0; j < p[0]; j++)
  {
	//performa ars on log(sigma_j), which is still saved in sigma_j	
// 	if (log(vardeltas[j]) > -20)	  
	  sample_ars (1, &sigmasbt[j], logpost_xi, -INFINITY, +INFINITY, 
				  log(vardeltas [j]/K[0]), 0);
	//convert xi to sigma_j
	sigmasbt [j] = exp (sigmasbt [j]); 
  }
}

/********* sampling from posterior of sigma given var and prior ************/
void spl_sgm_ghs (int p[1], double sigmasbt[p[0]], double vardeltas[p[0]], 
				  int K[1], double alpha[1], double log_aw [1])
{
  int j;
  
  //define the function computing logf and dlogf for ars
  void logpost_xi (double xi, double logf[1], double dlogf[1])
  {
	double vexi, eximinlogaw;
	
	logf[0]  = - (K[0]  - 1.0) / 2.0 * xi ;
	dlogf[0] = - (K[0]  - 1.0) / 2.0 ;
	
	vexi = vardeltas [j] / 2.0 / exp (xi);
	logf[0]  -= vexi;
	dlogf[0] += vexi;
	
	eximinlogaw = exp(xi - log_aw[0]);
	logf[0]  -= (alpha[0] + 1.0) / 2.0 * log (1.0 + eximinlogaw);
	dlogf[0] -= (alpha[0] + 1.0) / 2.0 * eximinlogaw / (1.0 + eximinlogaw);
  }
  
  for (j = 0; j < p[0]; j++)
  {
	//performa ars on log(sigma_j), which is still saved in sigma_j
	sample_ars (1, &sigmasbt[j], logpost_xi, -INFINITY, +INFINITY, 
				log(vardeltas [j]/K[0]), 0);
	//convert xi to sigma_j
	sigmasbt [j] = exp (sigmasbt [j]); 
  }
}


void samplew (double logw[1], int K[1], int p[1], double vardeltas[p[0]],
              double nu[1], double s[1], double eta[1])
{
    void logpost_logw (double u, double logf[1], double dlogf[1])
    {
        double w, sdu, wnu;
        int j;

        w = exp (u);
        sdu = (u - s[0])/eta[0];
        wnu = w * nu[0];
        // loglikelihood
        dlogf [0] = 0;
        logf[0] = 0;
        for (j = 0; j < p[0]; j++)
        {
          logf [0] += log (vardeltas[j] + wnu);
          dlogf[0] += wnu / (vardeltas[j] + wnu);
        }
        logf [0] *= (- (nu[0] + K[0]) / 2);
        dlogf[0] *= (- (nu[0] + K[0]) / 2);
        logf [0] += p[0] * nu[0] / 2 * u ;
        dlogf[0] += p[0] * nu[0] / 2;    
        logf [0] += - R_pow_di (sdu,2) / 2 - log (eta[0]);
        dlogf[0] += - sdu / eta[0];
    }
    
    if (eta[0] < 0.01)
    {
        logw[0] = s[0];
    }
    else
    {
        sample_ars (1, logw, logpost_logw, -INFINITY, +INFINITY, logw[0], 0);
    }
}

void find_normlv (int n[1], int C[1], 
                  double lv [n[0]][C[0]], double normlv[n[0]][C[0]] )
{
  int i,c;
  double lse[1];

  for (i = 0; i < n[0]; i++)
  {
    log_sum_exp (C, lv[i], lse);
    for (c = 0; c < C[0]; c++)
    {
      normlv[i][c] = lv[i][c] - lse[0];
    }
  }
}

void htlr_fit (        
        int p[1], int K[1], int n[1], 
        double X [p[0]+1][n[0]], double ymat[K[0]][n[0]], int ybase[n[0]],  
        // prior
        char* ptype[], double alpha [1], double s [1], double eta [1], 
        double sigmab0 [1], 
        // sampling
        int iters_rmc [1],int thin [1],int leap_L [1],   
        int iters_h [1], int leap_L_h [1], double leap_step [1],  
        double hmc_sgmcut [1], double DDNloglike [p[0] + 1],
        // fit result
        double mcdeltas  [iters_rmc[0] + 1] [K[0]] [p[0] + 1], 
        double mclogw [iters_rmc [0] + 1], 
        double mcsigmasbt [iters_rmc [0] + 1] [p[0]+1], 
        double mcvardeltas [iters_rmc [0] + 1] [p[0]+1], 
        double mcloglike [iters_rmc [0] + 1], 
        double mcuvar[iters_rmc[0] + 1],
        double mchmcrej [iters_rmc[0] + 1],
        // other control or result
        int silence [1],  int looklf [1] )
{
  int i_trj, k, c, i, j, C[1] = {K[0] + 1}, nvar[1] = {p[0] + 1};
  int ids_update [nvar[0]], ids_fix[nvar[0]], isfault = 0, 
      nuvar, nfvar, uj, fj, i_rmc, i_thin, L, i_mc;
      
  double lv [n[0]][C[0]], lv_old [n[0]][C[0]], lv_fix [n[0]][C[0]],
         normlv [n[0]][C[0]], loglike[1], loglike_old[1],
         predprob [n[0]][C[0]], predprob_old [n[0]][C[0]], 
         DNloglike [K[0]][nvar[0]], DNloglike_old [K[0]][nvar[0]], 
         deltas [K[0]][nvar[0]],deltas_old [K[0]][nvar[0]],
         momt[K[0]][nvar[0]],
         SUMsqdeltas[nvar[0]], SUMsqdeltas_old[nvar[0]], 
         SUMdeltas [nvar[0]], SUMdeltas_old [nvar[0]], 
         vardeltas [nvar[0]], vardeltas_old [nvar[0]], 
         DNlogprior [K[0]][nvar[0]], DNlogprior_old [K[0]][nvar[0]],
         DNlogpost [K[0]][nvar[0]], 
         sigmasbt [nvar[0]], logw[1],
         logprior, logprior_momt, nenergy, nenergy_old, aw, log_aw,
         stepsizes[nvar[0]], alpha_post = {(alpha[0] + K[0])/2},
         sgmsqcut, no_uvar, rej, log_alpha = {log (alpha[0])};

  if (hmc_sgmcut[0] > 0) sgmsqcut = R_pow_di(hmc_sgmcut [0],2);
  else sgmsqcut = hmc_sgmcut[0];

  // this function determines which features to be updated 
  void whichupdate (double cut)
  {
    nuvar = 0;
    nfvar = 0;
    for (j = 0; j < nvar[0]; j++) 
    {
        if (sigmasbt[j] > cut) ids_update [nuvar++] = j;
        else ids_fix [nfvar++] = j;
    }
  }
  
  // this function updates lv and predprob
  void updatepredprob (void)
  {
    // add fixed part
    for (i = 0; i < n[0]; i++) cpvec (K[0], &lv_fix[i][1], &lv[i][1]);
    // add updated part
    for (uj = 0; uj < nuvar; uj++) 
    {
         j = ids_update[uj];
         for (k = 0; k < K[0]; k++) 
         {
            for (i = 0; i < n[0]; i++)
            {
                lv [i][k+1] +=  deltas[k][j] * X [j][i];
            }
         }
    }
    
    find_normlv (n, C, lv, normlv);
    
    for (i = 0; i < n[0]; i++) 
    for (c = 0; c < C[0]; c++) 
    {
      predprob [i][c] =  exp (normlv[i][c]) ;
    }
  }
  
  void detach_fixlv ()
  {

    if (nuvar <= nvar[0]/2 )
    {
        for (i = 0; i < n[0]; i++)
        {
           for (c = 1; c < C[0]; c++)  lv_fix[i][c] = lv[i][c];
        } 

        // remove updated part
        for (uj = 0; uj < nuvar; uj++) 
        {
             j = ids_update[uj];
             for (k = 0; k < K[0]; k++) 
             {
                for (i = 0; i < n[0]; i++)
                {
                    lv_fix [i][k+1] -=  deltas[k][j] * X [j][i];
                }
             }
        }
    }
    else
    {
        for (i = 0; i < n[0]; i++) setvec (K[0], &lv_fix[i][1], 0);

        // add fixed part
        for (fj = 0; fj < nfvar; fj++) 
        {   
           j = ids_fix[fj];
           
           for (k = 0; k < K[0]; k++) 
           {
             for (i = 0; i < n[0]; i++)
             {
                lv_fix [i][k+1] +=  deltas[k][j] * X [j][i];
             }
           }
        }
    }
  }

  // this function computes derivatives of negative log like
  void updateDNloglike (void)
  {
     for (uj = 0; uj < nuvar; uj++) 
     {
        j = ids_update[uj];
        for (k = 0; k < K[0]; k++)
        {
            DNloglike[k][j] = 0; 
            for (i = 0; i < n[0]; i++) 
            {
              DNloglike[k][j] += X[j][i] * (predprob [i][k+1] - ymat [k][i]);  
            }
        }
     }
  }
  
  // this function computes log likelihood given normlv
  void updateloglike (void)
  {
    loglike[0] = 0;
    for (i = 0; i < n[0]; i++) loglike [0] += normlv [i][ybase[i]];
  }
  
  // this function update SUMdeltas and DNlogprior
  void updateDNlogprior ()
  {
    for (uj = 0; uj < nuvar; uj++) 
    {
      j = ids_update[uj];
      
      //sum deltas for each feature
      SUMdeltas [j] = 0;
      for (k = 0; k < K[0]; k++) SUMdeltas[j] += deltas[k][j];
      
      for (k = 0; k < K[0]; k++)
      {    
        DNlogprior [k][j] = deltas [k][j] - SUMdeltas[j]/C[0];
      }
    }
  }

  void updateDNlogpost ()
  {       
     for (uj = 0; uj < nuvar; uj++) 
     {
       j = ids_update [uj];
       for (k = 0; k < K[0]; k++)
       {    
         DNlogpost [k][j] = DNloglike [k][j] + DNlogprior [k][j] / sigmasbt[j];
       }
     }
  }
  // this function updates SUMsqdeltas, vardeltas
  void updatevardeltas (void)
  {
    for (uj = 0; uj < nuvar; uj++) 
    {
        j = ids_update[uj];
        // prior deltas
        SUMsqdeltas[j] = 0;
        
        for (k = 0; k < K[0]; k++) SUMsqdeltas[j] += R_pow_di (deltas[k][j],2);
        
        vardeltas[j] = SUMsqdeltas[j] - R_pow_di(SUMdeltas[j],2) / C[0];
    }
  }
  // this function compute neg. energy
  double comp_nenergy (void)
  {
    
    logprior = 0;
    logprior_momt = 0;
    for (uj = 0; uj < nuvar; uj++) 
    {
        j = ids_update[uj];
        logprior += vardeltas[j]/sigmasbt[j];
        
        for (k = 0; k < K[0]; k++)  
          logprior_momt += R_pow_di(momt[k][j],2);
    }
    return (loglike[0] - logprior/2 - logprior_momt/2);
  }
  
  void gen_momt (void)
  {
    for (uj = 0; uj < nuvar; uj++) 
    {
        j = ids_update[uj];
        for (k = 0; k < K[0]; k++)
        {
          GetRNGstate ();
          momt[k][j] = rnorm (0,1);
          PutRNGstate ();
        }
    } 
  }
  
  void  cache_oldvalues (void)
  {
    cpvec (n[0] * C[0], &lv[0][0], &lv_old[0][0]);
    cpvec (n[0] * C[0], &predprob[0][0], &predprob_old[0][0]);
    loglike_old [0] = loglike[0];

    for (uj = 0; uj < nuvar; uj++) 
    {
        j = ids_update[uj];
        for (k = 0; k < K[0]; k++)
        {
            deltas_old[k][j] = deltas[k][j];
            DNlogprior_old[k][j] = DNlogprior[k][j];
        }
        vardeltas_old [j] = vardeltas [j];
    }
  }

  void  restore_oldvalues (void)
  {
    cpvec (n[0] * C[0], &lv_old[0][0], &lv[0][0]);
    cpvec (n[0] * C[0], &predprob_old[0][0], &predprob[0][0]);
    loglike [0] = loglike_old[0];

    for (uj = 0; uj < nuvar; uj++) 
    {
        j = ids_update[uj];
        for (k = 0; k < K[0]; k++)
        {
            deltas [k][j] = deltas_old [k][j];
            DNlogprior[k][j] = DNlogprior_old[k][j];
        }
        vardeltas [j] = vardeltas_old [j];
    }
  }
  
  void updatestepsizes (void)
  {
    for (uj = 0; uj < nuvar; uj++) 
    {
        j = ids_update[uj];
        stepsizes [j] = leap_step [0];
        stepsizes [j] /=  sqrt(DDNloglike[j] + K[0]/sigmasbt[j]/C[0]);
    }
  }
  
  /*********************** getting initial values ************************/
  cpvec (K[0] * nvar[0], &mcdeltas[0][0][0], &deltas[0][0]); 
  logw[0] = mclogw[0];
  cpvec (nvar[0], &mcsigmasbt[0][0], sigmasbt);  
  setvec (n[0] * C[0], &lv_fix [0][0], 0);
  setvec (n[0] * C[0], &lv [0][0], 0);
  // compute some initial values
  whichupdate (-1);  // set to update all
  updatepredprob (); // lv is computed here
  updateloglike (); 
  updateDNlogprior ();
  updatevardeltas();  
  mcloglike [0] = loglike[0];
  cpvec (nvar[0], vardeltas, &mcvardeltas[0][0]);
  
  /************************ start gibbs sampling **************************/
  for (i_mc = 0; i_mc < iters_h[0] + iters_rmc[0]; i_mc ++)
  {
  
    i_rmc = i_mc - iters_h[0];

    if (i_mc < iters_h[0] / 2.0)
    {
        L = leap_L_h[0];
        logw[0] = -10;
    }  else if (i_mc < iters_h[0]) 
    {
        L = leap_L_h [0];
        logw[0] = s[0];
    } else 
    {
        L = leap_L [0];
        logw[0] = s[0];        
    }
    
    /***************** thin iterations of Gibbs sampling ******************/ 
    no_uvar = 0;
    rej = 0;
    for (i_thin = 0; i_thin < thin[0]; i_thin ++)
    {
        /*********************** HMC Metropolis Update ********************/    
        // initialize HMC
        whichupdate (sgmsqcut);
        no_uvar += nuvar;
        
        gen_momt ();
        updatestepsizes ();        
        detach_fixlv ();
        cache_oldvalues ();
        nenergy_old = comp_nenergy ();

        // start trajectory
        updateDNlogprior ();
        updateDNloglike (); // recompute derivatives of log likelihood
        updateDNlogpost (); // recompute derivatives of log prior       
        for (i_trj = 0; i_trj < L; i_trj ++)
        {
          for (uj = 0; uj < nuvar; uj++) 
          { 
            j = ids_update [uj];
            for (k = 0; k < K[0]; k++)
            {
                momt [k][j]   -= stepsizes [j] / 2 * DNlogpost [k][j];
                deltas [k][j] += stepsizes [j] * momt [k][j];
            }
          }
          // compute derivative of minus log joint distribution
          updatepredprob ();
          updateDNlogprior ();
          updateDNloglike ();
          updateDNlogpost ();
          // move momonton with new derivatives
          for (uj = 0; uj < nuvar; uj++) 
          {
            j = ids_update [uj];
            for (k = 0; k < K[0]; k++)
            {
                momt[k][j] -= stepsizes[j] / 2 * DNlogpost[k][j];
            }
          }
       }

       isfault = 0;
       for (uj = 0; uj < nuvar; uj ++)
       {
         j = ids_update [uj];
         for (k = 0; k < K[0]; k++)
         {
            if (fabs (deltas[k][j]) > 20 )
            {
                isfault = 1;
                break;
            }
         }
 
         if (isfault == 1) break;
 
       }       
       
       // decide whether to accept it
       updateloglike ();
       updatevardeltas ();
       nenergy = comp_nenergy ();

       GetRNGstate ();
       if (log (runif(0,1)) > nenergy - nenergy_old || isfault == 1) 
       {
           restore_oldvalues ();
           rej ++;
       }
       PutRNGstate ();
       
       /*********************** Sigmas Update  *****************************/
       aw = alpha[0] * exp (logw[0]);
       log_aw = logw[0] + log_alpha;
     
       if (strcmp (ptype[0], "t") == 0 )
       {
          for (j = 1; j < nvar[0]; j++)
          {       
            GetRNGstate ();
            sigmasbt [j]=1.0/rgamma (alpha_post,1.0) * (aw + vardeltas[j])/2.0; 
            PutRNGstate ();
          }
       }
       
       if (strcmp (ptype[0], "ghs") == 0 )
       {
		  GetRNGstate ();
		  spl_sgm_ghs (p, &sigmasbt[1],&vardeltas[1], K, alpha, &log_aw);
          PutRNGstate ();  
       }
       
       if (strcmp (ptype[0], "neg") == 0 )
       {
		  GetRNGstate ();
		  spl_sgm_neg (p, &sigmasbt[1],&vardeltas[1], K, alpha, &log_aw);
          PutRNGstate ();       
       }
       
       /********************** logw  Update  ******************************/
       if (strcmp (ptype[0], "t") == 0 & eta[0] > 1E-10)
            samplew (logw, K, p, &vardeltas[1], alpha, s, eta);
    }
    
    no_uvar /= thin [0];
    rej /= thin [0];
    
    /****************** record the markov chain state ********************/
    if (i_rmc >= 0 )
    {
        cpvec (K[0] * nvar[0], &deltas[0][0], &mcdeltas[i_rmc + 1][0][0]);
        cpvec (nvar[0], sigmasbt, &mcsigmasbt[i_rmc + 1][0]);
        cpvec (nvar[0], vardeltas, &mcvardeltas[i_rmc + 1][0]);
        mclogw[i_rmc + 1] = logw[0];
        mcloglike [i_rmc + 1] = loglike[0];
        mcuvar [i_rmc + 1] = no_uvar;
        mchmcrej [i_rmc + 1] = rej;
    }
    
    // print some results on screen
    if (silence[0] == 0)
    {
       Rprintf (
       "Iter%4d: deviance=%5.3f, logw=%6.2f, nuvar=%3.0f, hmcrej=%4.2f\n", 
       i_rmc, -loglike[0]/n[0],logw[0], no_uvar, rej);
    }
  }
}



