/***********************************************************************/
double find_double_max(int la,double a[]){
    int i;
    int m = a[0];
    for(i = 1; i < la; i++){
       if(a[i] > m) m = a[i];
    }
    return(m);
}

/***********************************************************************/
double log_sum_exp_c (int c, int nr, int nc, double M[nr][nc])
{
   int i;
   double m, out;
   
   //find maximum value
   m = M[0][c];
   for(i = 1; i < nr; i++) if(M[i][c] > m) m = M[i][c]; 
   
   out = 0;
   for(i = 0; i < nr; i++) {
      out += exp(M[i][c] - m); 
   }
   return(m + log(out));
}

/***********************************************************************/
double log_sum_exp1(int la, double a[])
{  int i;
   double m,out=0;
   m = find_double_max(la,a);
   for(i = 0; i < la; i++){
      out += exp(a[i] - m);
   }
   return (log(out) + m);
}

/***********************************************************************/
void pred_ht (
       int n [1], int k [1], int G [1], int spls [1], 
       double features [k[0]] [n[0]], double mu [spls[0]] [G[0]] [k[0]], 
       double  sd_x [spls[0]] [k[0]], double log_freq_y [G[0]], 
       double predprobs [G[0]] [n[0]] )
{
     double predprobs_spl [G[0]] [spls[0]], log_sum_probs;
     int i_sp, i_case, i_var, i_g;
     
     for(i_case = 0; i_case < n[0]; i_case++){
        for(i_g = 0; i_g < G[0]; i_g ++){
            //compute predictive probs given each sample
            for(i_sp = 0; i_sp < spls[0]; i_sp++){
                predprobs_spl [i_g] [i_sp] = 0;
                for(i_var = 0; i_var < k[0]; i_var++) {
                  predprobs_spl [i_g] [i_sp] += dnorm(features[i_var][i_case],
                     mu[i_sp][i_g][i_var],sd_x[i_sp][i_var],1);
                }
            }
            //sum predictive probs over samples, plus posterior of y
            predprobs [i_g] [i_case] = 
                log_sum_exp1 ( spls[0], &predprobs_spl[i_g][0]) + 
                log_freq_y [i_g];
        }
        //sum predictive probs over groups
        log_sum_probs = log_sum_exp_c (i_case, G[0], n[0], predprobs);
        //normalize predictive probs
        for(i_g = 0; i_g < G[0]; i_g ++) {
            predprobs [i_g] [i_case] = 
              exp( predprobs [i_g] [i_case] -log_sum_probs );
        }
     }
}

