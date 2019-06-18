

void leapfrog (
  // integers contolling dimensions of matrices
  int looklf[1], int L[1], int K[1], int n[1], int nvar[1],
  // arugments along trajectory
  double deltas[K[0]][nvar[0]], double momt [K[0]][nvar[0]],
  double ddeltas [1], double loglike [1], double vardeltas [nvar[0]],
  double DNlogprior [K[0]][nvar[0]], double DNloglike [K[0]][nvar[0]],
  //arguments for leapfrog trajectory
  double logliketrj [L[0] + 1], double nenergytrj [L[0]+ 1],
  double ddeltastrj [L[0] + 1],
  // arguments for data and parameters for leapfrog
  double tnewX [nvar[0]][n[0]], double ymat [K[0]][n[0]], int y[n[0]],
  double stepsizes [nvar[0]], double sigmas [nvar[0]] )
{
  int i_trj, k, c, i, j, C[1] = {K[0] + 1};
  double lv [n[0]][C[0]], predprob [n[0]][C[0]], DNlogpost [K[0]][nvar[0]],
         SUMsqdeltas[nvar[0]], SUMdeltas [nvar[0]],
         svardeltas [1], smomtsq[1];

  for (k = 0; k < K[0]; k++)
  for (j = 0; j < nvar[0]; j++) {
    DNlogpost [k][j] = DNloglike[k][j] + DNlogprior[k][j]/sigmas[j];
  }

  // start leapfrog trajectory
  for (i_trj = 1; i_trj <= L[0]; i_trj++)
  {
    // move momonton and deltas with current derivatives
    for (k = 0; k < K[0]; k++)
    for (j = 0; j < nvar[0]; j++) {
      momt [k][j]   -= stepsizes [j] / 2 * DNlogpost [k][j];
      deltas [k][j] += stepsizes [j] * momt [k][j];
    }

    //compute derivative of minus log joint distribution
    colSums (K, nvar, deltas, SUMdeltas);
    for (k = 0; k < K[0]; k++)
    for (j = 0; j < nvar[0]; j++) {
      DNlogprior [k][j] = deltas [k][j] - SUMdeltas[j]/C[0];
    }

    // compute lv
    for (i = 0; i < n[0]; i++) lv[i][0] = 0;
    for (i = 0; i < n[0]; i++)
    for (c = 1; c < C[0]; c++) {
      lv[i][c] = 0; for (j = 0; j < nvar[0]; j++) {
          lv [i][c] +=  deltas[c-1][j] * tnewX[j][i];
      }
    }
    
    // normalize lv
    norm_lv (n, C, &lv[0][0]);

    //find predictive probabilities for training cases
    for (i = 0; i < n[0]; i++)
    for (c = 0; c < C[0]; c++) {
      predprob [i][c] = exp (lv[i][c]);
    }

    // compute derivatives of negative log like and log post
    // and move momonton with new derivatives
    for (k = 0; k < K[0]; k++)
    for (j = 0; j < nvar[0]; j++) {
      
      // compute derivative of log like
      DNloglike[k][j] = 0; for (i = 0; i < n[0]; i++) {
      DNloglike[k][j] += tnewX [j][i] * (predprob [i][k+1] - ymat [k][i]);
      }
      // compute derivative of log post
      DNlogpost [k][j] = DNloglike [k][j] + DNlogprior [k][j] / sigmas [j];

      // move momonton with new derivatives
      momt[k][j] -= stepsizes[j] / 2 * DNlogpost[k][j];
    }

    // computing log prior and log like, and energy
    if (i_trj == L[0] || looklf[0] == 1)
    {
        // compute vars and distance of deltas
        ddeltas [0] = 0;
        svardeltas [0] = 0;
        for (j = 0; j < nvar[0]; j++) {
          SUMsqdeltas[j] = 0;
          for (k = 0; k < K[0]; k++) SUMsqdeltas[j] += R_pow_di(deltas[k][j],2);

          vardeltas[j] = SUMsqdeltas[j] - R_pow_di(SUMdeltas[j],2) / C[0];
          svardeltas[0] += vardeltas[j]/sigmas[j];
          ddeltas[0] += SUMsqdeltas[j];
        }
        ddeltastrj[i_trj] = ddeltas[0]; 
        
        // compute distance of momonton
        smomtsq [0] = 0;
        for (j = 0; j < nvar[0]; j++)
        for (k = 0; k < K[0]; k++) smomtsq [0] +=  R_pow_di(momt[k][j],2);

        //compute log likelihood
        loglike [0] = 0;
        for (i = 0; i < n[0]; i++) loglike [0] += lv[i][y[i]];
        logliketrj [i_trj] = loglike[0];

        //compute negative energy
        nenergytrj [i_trj] = loglike[0] - svardeltas[0]/2 - smomtsq[0]/2;
    }
  }
}

