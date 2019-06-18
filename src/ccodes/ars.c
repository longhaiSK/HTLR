#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdio.h>

// NOTE: in this program, I call a piece of linear function as 'hull', and
// the piecewise linear function above logf upperhulls, and
// the piecewise linear function below logf lowerhulls

/////////////////////////////////////////////////////////////////////////////
// header of this file
void sample_disc (int n[1], int rn[n[0]], int k[1], double lw[k[0]]);
void sample_elin
   (int n[1], double rn[n[0]], double lower[1], double upper[1],
    double dlogf[1], double tol_dlogf_is0[1]);
void logint_elin
  (double logf[1],double dlogf[1],double t[1],double lower[1],
   double upper[1], double lw[1], double tol_dlogf_is0[1]);
double interc (double t1[1], double t2[1],
               double logf1[1], double logf2[1],
               double dlogf1[1], double dlogf2[1], double tol_ddlogf_is0[1] );
/////////////////////////////////////////////////////////////////////////////

///////////////////////  the main ARS function //////////////////////////////
void sample_ars
    (int n, // number of samples
     double rn[n],  // sampling output
     void (*eval_logf) (double x, double logf[1], double dlogf[1]),
     // pointer to another function for evaluating logf and dlogf
     double lb, double ub, // lower and upper bounds of logf
     double ini_tpoint, //initial tangent point
     double dars // flag whether to print ARS information
     )
{
  /************************** Optional User Settings ***********************/
  //set maximum number of pieces of linear hulls
  int max_nhull = 1000;
  //set smallest difference of derivatives that can be thought of as the same
  //so that we don't need to insert one more hull
  double tol_dlogf_is0 = 1E-5, tol_ddlogf_is0 = 1E-5;
  double stepout = 10; // size of stepout in initializing linear hulls
  ////////////////////////////// end of user settings ///////////////////////

  // global working vectors used to represent linear hulls
  // an element corresponding to a tangent point
  double
  tpoints [max_nhull], // tangent points
  lws [max_nhull], // log integrals of exp hulls
  lowerbounds[max_nhull], upperbounds[max_nhull], //bounds of hulls
  //values of logf and dlogf at tangent points
  logfvs[max_nhull], dlogfvs[max_nhull],
  //slopes of left and right squeezings
  slopes_leftsq[max_nhull],slopes_ritesq[max_nhull];
  //indice of left and right hulls
  int lefthulls [max_nhull], ritehulls [max_nhull],
  no_hulls; // total number of hulls

  //construct the first hull
  tpoints [0] = ini_tpoint; //the tangent point
  eval_logf (tpoints[0], &logfvs[0], &dlogfvs[0]);
  if (!isfinite (logfvs[0]))
  {
    printf ("Error in adaptive rejection sampling:\n");
    printf ("the first tangent point doesn't have positive probability.\n");
    exit (1);
  }
  lowerbounds [0] = fmax2(lb, -INFINITY); // lower bound of the hull
  upperbounds [0] = fmin2(ub, +INFINITY); // upper bound of the hull
  lefthulls [0] = -1; // index of left hull
  ritehulls [0] = max_nhull; //index of right hull
  slopes_leftsq [0] = INFINITY; // slope of left squeezing arc
  slopes_ritesq [0] = -INFINITY; // slope of right sequeezing arc
  //compute log weights, updating lw[0]
  lws[0] = INFINITY;
  no_hulls = 1;

  // this function updates the envolop and squeezing functions.
  // newx --- new point to be inserted
  // h --- index of the hull where newx is from
  // logfv, dlogfv --- values of logf and dlogv at newx
  void update_hulls (int h, double newx, double logfv, double dlogfv)
  {
    int lh, rh, nh;

    if (no_hulls == max_nhull) return;// reaching the limit of working vector

    //specify left and right hulls of new hull
    if (newx > tpoints[h]) // to insert to the right of hull h
    {
      lh = h; rh = ritehulls[h];
      // if logfv is -infinity, only update the rightest hull rightbound and lw
      if (rh == max_nhull & logfv == - INFINITY)
      {
        if (upperbounds[h] != newx)
        {
          upperbounds [h] = newx;
          logint_elin (
            &logfvs[h], &dlogfvs[h],
            &tpoints[h], &lowerbounds[h],
            &upperbounds[h], &lws[h], &tol_dlogf_is0);
        }
        return;
      }
    }
    else // to insert to the left of hull h
    {
      lh = lefthulls[h]; rh = h;
      // if logfv is -infinity, only update the leftest hull leftbound and lw
      if (lh == -1 & logfv == - INFINITY)
      {
        if (lowerbounds [h] != newx)
        {
          lowerbounds [h] = newx;
          logint_elin (
            &logfvs[h], &dlogfvs[h],
            &tpoints[h], &lowerbounds[h],
            &upperbounds[h], &lws[h], &tol_dlogf_is0);
        }
        return;
      }
    }

    //insert a new hull
    nh = no_hulls;
    no_hulls ++;
    tpoints[nh] = newx;
    logfvs[nh] = logfv;
    dlogfvs[nh] = dlogfv;
    lefthulls[nh] = lh;
    ritehulls[nh] = rh;

    if (lh == -1) // nh will be the new leftest hull
    {
      lowerbounds [nh] = lowerbounds [h];
      slopes_leftsq [nh] = + INFINITY;
    }
    else
    {
      lowerbounds[nh] = interc (&tpoints[lh], &tpoints[nh],
        &logfvs[lh], &logfvs[nh], &dlogfvs[lh], &dlogfvs[nh],
        &tol_ddlogf_is0); // lowerbound
      slopes_leftsq [nh] = (logfvs[nh] - logfvs[lh]) /
                           (tpoints[nh] - tpoints[lh]);
    }
    if (rh == max_nhull)
    {
      upperbounds[nh] = upperbounds[h]; // upperbound
      slopes_ritesq[nh] = -INFINITY;
    }
    else
    {
      upperbounds[nh] = interc (&tpoints[nh], &tpoints[rh],
        &logfvs[nh], &logfvs[rh], &dlogfvs[nh], &dlogfvs[rh],
        &tol_ddlogf_is0); // upperbound
     slopes_ritesq [nh] = (logfvs[nh] - logfvs[rh]) /
                          (tpoints[nh] - tpoints[rh]);
    }

    logint_elin (
          &logfvs[nh], &dlogfvs[nh],
          &tpoints[nh], &lowerbounds[nh],
          &upperbounds[nh], &lws[nh], &tol_dlogf_is0);

    //update left hull of new null
    if (lh != -1)
    {
      upperbounds[lh] = lowerbounds[nh];
      ritehulls[lh] = nh;
      slopes_ritesq[lh] = slopes_leftsq [nh];
      logint_elin (
          &logfvs[lh], &dlogfvs[lh],
          &tpoints[lh], &lowerbounds[lh],
          &upperbounds[lh], &lws[lh], &tol_dlogf_is0);
    }

    //update right hull of newh if it exists
    if (rh != max_nhull)
    {
        lowerbounds[rh] = upperbounds[nh];
        lefthulls[rh] = nh;
        slopes_leftsq[rh] = slopes_ritesq [nh];

        logint_elin (
          &logfvs[rh], &dlogfvs[rh],
          &tpoints[rh], &lowerbounds[rh],
          &upperbounds[rh], &lws[rh], &tol_dlogf_is0);
    }
  }

  double newx, newlogf, newdlogf; // a new tangent point to be inserted
  int h; // index of the hull where newx is from
  int newh; // index of new inserted hull

  // if lb is finite, bound the first hull at left
  // or insert a hull tangent at lb if logf at lb is finite too
  if (isfinite (lb))
  {
    h = 0;
    newx = lb;
    eval_logf (newx, &newlogf, &newdlogf);
    update_hulls (h, newx, newlogf, newdlogf);
  }
  //expanding at the left until reaching a bound or integral to finite
  else
  {

    h = 0;
    newx = tpoints[0] - stepout;
    do
    {
      if (no_hulls == max_nhull)
      {
        printf ("Error in Rejection Sampling:\n");
        printf ("'max_nhull' is set too small, or your log-PDF NOT concave\n");
        exit (1);
      }
      eval_logf (newx, &newlogf, &newdlogf);
      update_hulls (h, newx, newlogf, newdlogf);
      // finding a new leftbound quite expanding
      if (newlogf == - INFINITY) break;
      newx -= stepout;
      h = no_hulls - 1;
    }
    while (newdlogf < tol_dlogf_is0);
  }

  // if ub is finite, bound the first hull at the right
  // or insert a hull tangent at ub if logf at ub is finite too
  if (isfinite (ub))
  {
    h = 0;
    newx = ub;
    eval_logf (newx, &newlogf, &newdlogf);
    update_hulls (h, newx, newlogf, newdlogf);
  }
  else // expanding at the right until reaching a bound or integral to finite
  {
    h = 0;
    newx = tpoints[0] + stepout;
    do
    {
      if (no_hulls == max_nhull)
      {
        printf ("Error in Rejection Sampling:\n");
        printf ("'max_nhull' is set too small, or your log-PDF NOT concave\n");
        exit (1);
      }

      eval_logf (newx, &newlogf, &newdlogf);
      update_hulls (h, newx, newlogf, newdlogf);
      if (!isfinite (newlogf)) break;
      newx += stepout;
      h = no_hulls - 1;
    }
    while (newdlogf > - tol_dlogf_is0);

  }

  double eval_upperhull (int h, double newx)
  {
    return ((newx - tpoints[h]) * dlogfvs[h] + logfvs[h]);
  }

  double eval_lowerhull (int h, double newx)
  {
     if (newx >= tpoints[h])
     {
       return ((newx - tpoints[h]) * slopes_ritesq[h] + logfvs[h]);
     }
     else
     {
       return ((newx - tpoints[h]) * slopes_leftsq[h] + logfvs[h]);
     }
  }

  //**************** Doing adaptive rejection sampling ********************//
  // define parameters used while sampling
  int one = 1, no_rejs = 0, i, rejected;
  double
  upperhullv, // value of upper hull at newx
  lowerhullv, // value of lower (squeezing) hull at newx
  u, // a random number used to determine whether to accept
  logacceptv; //if logacceptv is smaller than logf, newx  accepted ;

  for (i = 0; i < n; i++)
  {
    rejected = 1;
    while (rejected)
    {
      //draw a new point and a unif random number
      sample_disc (&one, &h, &no_hulls, lws);
      sample_elin (&one, &newx, &lowerbounds[h], &upperbounds[h],
                   &dlogfvs[h], &tol_dlogf_is0);
      upperhullv = eval_upperhull (h, newx);
      GetRNGstate (); u = unif_rand (); PutRNGstate ();
      logacceptv = upperhullv + log (u);
      lowerhullv = eval_lowerhull (h, newx);
      //check acceptance with squeezing function
      if (logacceptv <= lowerhullv)
      {
        rn [i] = newx;
        rejected = 0;
      }
      else
      {
        // check acceptance with logf
        // eval logf at newx and insert a new hull
        eval_logf (newx, &newlogf, &newdlogf);
        update_hulls (h, newx, newlogf, newdlogf);
        if (logacceptv <= newlogf)
        {
          rn [i] = newx;
          rejected = 0;
        }
        else no_rejs ++;
      }
    }
 }

 double rate_rej;
 if (dars == 1)
 {
  rate_rej = (no_rejs + 0.0) / (no_rejs + n + 0.0); // return rejection rate
  printf ("no of hulls = %d, rejection rate = %4.2f\n", no_hulls, rate_rej);
 }
}

// find maximum value in double vector a with length n
double fmaxm (int n, double a[n])
{
  double ma;
  int i;
  ma = a [0];
  if (n > 1)
  {
    for (i = 1; i < n; i++)
    {
      ma = fmax2 (a[i], ma);
    }

  }
  return (ma);
}

// n --- number of random numbers
// k --- number of discrete values
// lw --- log of probabilities
// rn --- vector of random numbers returned
void sample_disc (int n[1], int rn[n[0]], int k[1], double lw[k[0]])
{
  double cw[k[0]], u, max_lw;
  int i, j;

  // constructing probabilities from log probabilities
  max_lw = fmaxm (k[0], lw);
  cw[0] = exp(lw[0] - max_lw);
  for (i = 1; i < k[0]; i++) cw [i] = cw [i-1] + exp(lw [i] - max_lw);

  for (j = 0; j < n[0]; j++)
  {
    GetRNGstate ();
    u = unif_rand() * cw[k[0] - 1];
    PutRNGstate ();
    // convert u into a discrete value
    for (i = 0; i < k[0]; i++)
    {
      if (u <= cw[i])
      {
        rn[j] = i;
        break;
      }
    }
  }
}

// this function samples from: exp (a[0]*x) I (x in [lb, upper[0]])
// n is number of random numbers required, rn will store random numbers
void sample_elin
   (int n[1], double rn[n[0]], double lower[1], double upper[1],
    double dlogf[1], double tol_dlogf_is0[1])
{
  // set smallest value for derivative that can be thought of as 0
  double  y, dx;
  int j, type_lin, isfault = 0;

  // checking linear function type and fault
  if (fabs (dlogf[0]) <= tol_dlogf_is0 [0])
  {
     if (! ( isfinite (lower[0]) & isfinite (upper[0])) )
     {
       isfault = 1;
     }
     else
     {
       type_lin = 0;
     }
  }

  if (dlogf[0] >  tol_dlogf_is0 [0])
  {
    if (!isfinite (upper[0]))
    {
      isfault = 1;
    }
    else
    {
      type_lin = 1;
    }
  }

  if (dlogf[0] < -tol_dlogf_is0 [0])
  {
    if(!isfinite (lower[0]))
    {
      isfault = 1;
    }
    else
    {
      type_lin = -1;
    }
  }

  if (isfault)
  {
    printf ( "Error: in C function 'sample_elin':\n");
    printf ( "the exp linear function integrates to NAN/INFINITY\n");
    printf ( "(dlogf = %4.2f, lowerbound = %4.2f, upperbound = %4.2f)\n",
             dlogf[0], lower[0], upper[0]);
    exit (1);
  }

  dx = upper[0] - lower[0];

  for (j = 0; j < n[0]; j++)
  {
    GetRNGstate ();
    y = runif (0,1);
    PutRNGstate ();

    //converting uniform random number
    if (type_lin == 0 ) // slope is 0
      rn[j] = lower[0] + y * dx;

    if (type_lin == 1) // slope is postive
      rn[j] = upper[0] +
         log ( (1 - y) * exp (- dlogf[0] * dx) + y) / dlogf[0];

    if (type_lin == -1) //slope is negative
      rn[j] = lower[0] +
         log ( 1 - y + y * exp (dlogf[0] * dx) ) / dlogf[0];
  }
}

// this function evaluates the log of integral of exp linear hull
// logf --- value of linear hull at t
// dlogf --- value of derive of linear hull
// t --- tangent point where logf is calculated
// lower and upper --- lower and upper bounds of linear hull
// lw --- saved returned log integral
void logint_elin
  (double logf[1],double dlogf[1],double t[1],
   double lower[1],double upper[1], double lw[1],
   double tol_dlogf_is0[1])
{
  double  dx, abs_dlogf;
  int sgn_dlogf;

  dx = upper [0] - lower [0];
  abs_dlogf = fabs (dlogf[0]);


  if (abs_dlogf <= tol_dlogf_is0 [0]) // slope is 0
  {
    lw[0] = logf[0] + log (dx);
    return;
  }

  if (dlogf[0] > tol_dlogf_is0 [0]) sgn_dlogf = 1; else sgn_dlogf = -1;

  if (sgn_dlogf == 1) // slope is positive
  {
    lw[0] = logf[0] + dlogf[0] * (upper[0] - t[0]) - log (abs_dlogf)  +
            log (1 - exp (- abs_dlogf * dx) );
    return;
  }

  if (sgn_dlogf == -1) //slope is negative
  {
    lw[0] = logf[0] + dlogf[0] * (lower[0] - t[0]) - log (abs_dlogf)  +
            log (1 - exp (- abs_dlogf * dx) );
    return;
  }
}

// this function finds interception points between t1 and t2
double interc (
  double t1[1], double t2[1],
  double logf1[1], double logf2[1], double dlogf1[1], double dlogf2[1],
  double tol_ddlogf_is0[1] )
{
  if (fabs (dlogf1[0]-dlogf2[0]) > tol_ddlogf_is0 [0])
  {
    return (
      (logf2[0] - logf1[0] - dlogf2[0] * t2[0] + dlogf1[0] * t1[0]) /
      (dlogf1[0] - dlogf2[0])
    );
  }
  else
  {
    return ((t1[0] + t2[0])/2.0);
  }

}
