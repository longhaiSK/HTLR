#include <cmath>
#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include "ars.h"

// NOTE: in this program, I call a piece of linear function as 'hull', and
// the piecewise linear function above logf upperhulls, and
// the piecewise linear function below logf lowerhulls

// Header of helper functions
void sample_disc(const int n, int *rn, const int k, double *lw);
void sample_elin(const int n, double *rn, const double lower, const double upper,
                 double dlogf[1], double tol_dlogf_is0[1]);
void logint_elin(double logf[1], double dlogf[1], double t[1],
                 double lower[1], double upper[1], double lw[1],
                 double tol_dlogf_is0[1]);
double interc(double t1[1], double t2[1],
              double logf1[1], double logf2[1], double dlogf1[1], double dlogf2[1],
              double tol_ddlogf_is0[1]);

 
 
void SampleTarget::eval_logf(const double x, double &logf, double &dlogf) {}; 

// this function updates the envolop and squeezing functions.
// newx --- new point to be inserted
// h --- index of the hull where newx is from
// logfv, dlogfv --- values of logf and dlogv at newx
void ARS::update_hulls(const int h, const double newx, const double logfv, const double dlogfv)
{
  int lh, rh, nh;

  if (no_hulls == max_nhull) return;// reaching the limit of working vector

  //specify left and right hulls of new hull
  if (newx > tpoints[h]) // to insert to the right of hull h
  {
    lh = h; rh = ritehulls[h];
    // if logfv is -infinity, only update the rightest hull rightbound and lw
    if ((rh == max_nhull) & (logfv == -INFINITY))
    {
      if (upperbounds[h] != newx)
      {
        upperbounds[h] = newx;
        logint_elin(
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
    if ((lh == -1) & (logfv == -INFINITY))
    {
      if (lowerbounds[h] != newx)
      {
        lowerbounds[h] = newx;
        logint_elin(
          &logfvs[h], &dlogfvs[h],
          &tpoints[h], &lowerbounds[h],
          &upperbounds[h], &lws[h], &tol_dlogf_is0);
      }
      return;
    }
  }

  // insert a new hull
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
    lowerbounds[nh] = interc(&tpoints[lh], &tpoints[nh],
      &logfvs[lh], &logfvs[nh], &dlogfvs[lh], &dlogfvs[nh],
      &tol_ddlogf_is0); // lowerbound
    slopes_leftsq[nh] = (logfvs[nh] - logfvs[lh]) /
                        (tpoints[nh] - tpoints[lh]);
  }
  if (rh == max_nhull)
  {
    upperbounds[nh] = upperbounds[h]; // upperbound
    slopes_ritesq[nh] = -INFINITY;
  }
  else
  {
    upperbounds[nh] = interc(&tpoints[nh], &tpoints[rh],
      &logfvs[nh], &logfvs[rh], &dlogfvs[nh], &dlogfvs[rh],
      &tol_ddlogf_is0); // upperbound
    slopes_ritesq [nh] = (logfvs[nh] - logfvs[rh]) /
                        (tpoints[nh] - tpoints[rh]);
  }

  logint_elin(
        &logfvs[nh], &dlogfvs[nh],
        &tpoints[nh], &lowerbounds[nh],
        &upperbounds[nh], &lws[nh], &tol_dlogf_is0);

  // update left hull of new null
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

  // update right hull of newh if it exists
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

double ARS::eval_upperhull(const int h, const double newx)
{
  return ((newx - tpoints[h]) * dlogfvs[h] + logfvs[h]);
}

double ARS::eval_lowerhull(const int h, const double newx)
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

void ARS::Prepare()
{
  // if lb is finite, bound the first hull at left
  // or insert a hull tangent at lb if logf at lb is finite too
  if (std::isfinite (lb))
  {
    h = 0;
    newx = lb;
    target.eval_logf(newx, newlogf, newdlogf);
    update_hulls(h, newx, newlogf, newdlogf);
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
        printf("Error in Rejection Sampling:\n");
        printf("'max_nhull' is set too small, or your log-PDF NOT concave.\n");
        exit (1);
      }
      target.eval_logf(newx, newlogf, newdlogf);
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
  if (std::isfinite (ub))
  {
    h = 0;
    newx = ub;
    target.eval_logf(newx, newlogf, newdlogf);
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
        printf("Error in Rejection Sampling:\n");
        printf("'max_nhull' is set too small, or your log-PDF NOT concave\n");
        exit(1);
      }
      target.eval_logf(newx, newlogf, newdlogf);
      update_hulls (h, newx, newlogf, newdlogf);
      if (!std::isfinite(newlogf)) break;
      newx += stepout;
      h = no_hulls - 1;
    }
    while (newdlogf > - tol_dlogf_is0);
  }
}

ARS::ARS(int n, SampleTarget target, double ini_tpoint, 
    double lb/*= -INFINITY*/, double ub/*= +INFINITY*/,   
    bool verbose/*=false*/, int max_nhull/*=1000*/, double stepout/*=10*/,
    double tol_dlogf_is0/*= 1E-5*/, double tol_ddlogf_is0/*= 1E-5*/)
    : n(n), lb(lb), ub(ub), 
    verbose(verbose), max_nhull(max_nhull), stepout(stepout), 
    tol_dlogf_is0(tol_dlogf_is0), tol_ddlogf_is0(tol_ddlogf_is0), target(target)  
{
  // construct the first hull
  this->logfvs  = new double[max_nhull] {0};
  this->dlogfvs = new double[max_nhull] {0};
  this->tpoints = new double[max_nhull] {0};
  this->tpoints[0] = ini_tpoint; // the tangent point
  this->target.eval_logf(tpoints[0], logfvs[0], dlogfvs[0]);
  if (!std::isfinite(logfvs[0]))
  {
    printf("Error in Rejection Sampling:\n");
    printf("'max_nhull' is set too small, or your log-PDF NOT concave.\n");
    exit(1);
  }

  this->lowerbounds = new double[max_nhull]{0};
  this->upperbounds = new double[max_nhull]{0};
  lowerbounds[0] = fmax(lb, -INFINITY); // lower bound of the hull
  upperbounds[0] = fmin(ub, +INFINITY); // upper bound of the hull

  this->lefthulls = new int[max_nhull]{0};
  this->ritehulls = new int[max_nhull]{0};
  lefthulls[0] = -1;        // index of left hull
  ritehulls[0] = max_nhull; // index of right hull

  this->slopes_leftsq = new double[max_nhull]{0};
  this->slopes_ritesq = new double[max_nhull]{0};
  slopes_leftsq[0] = +INFINITY; // slope of left squeezing arc
  slopes_ritesq[0] = -INFINITY; // slope of right sequeezing arc

  this->lws = new double[max_nhull]{0};
  lws[0] = INFINITY; // compute log weights, updating lws[0]
  this->no_hulls = 1;
}

// Do adaptive rejection sampling
double *ARS::Sampling()
{
  Prepare();
  // define parameters used while sampling
  int one = 1, no_rejs = 0;
  double
      upperhullv, // value of upper hull at newx
      lowerhullv, // value of lower (squeezing) hull at newx
      u,          // a random number used to determine whether to accept
      logacceptv; // if logacceptv is smaller than logf, newx  accepted ;

  for (int i = 0; i < n; i++)
  {
    bool rejected = true;
    while (rejected)
    {
      // draw a new point and a unif random number
      sample_disc(one, &h, no_hulls, lws);
      sample_elin(one, &newx, lowerbounds[h], upperbounds[h],
                  &dlogfvs[h], &tol_dlogf_is0);
      upperhullv = eval_upperhull(h, newx);
      GetRNGstate();
      u = unif_rand();
      PutRNGstate();
      logacceptv = upperhullv + log(u);
      lowerhullv = eval_lowerhull(h, newx);
      //check acceptance with squeezing function
      if (logacceptv <= lowerhullv)
      {
        rn[i] = newx;
        rejected = false;
      }
      else
      {
        // check acceptance with logf
        // eval logf at newx and insert a new hull
        target.eval_logf(newx, newlogf, newdlogf);
        update_hulls(h, newx, newlogf, newdlogf);
        if (logacceptv <= newlogf)
        {
          rn[i] = newx;
          rejected = false;
        }
        else
          no_rejs++;
      }
    }
  }
  if (verbose)
  {
    double rate_rej = (no_rejs + 0.0) / (no_rejs + n + 0.0); // return rejection rate
    printf("no of hulls = %d, rejection rate = %4.2f\n", no_hulls, rate_rej);
  }
  return rn;
}

// find maximum value in double vector a with length n
double fmaxm(const int n, const double *a)
{
  double ma = a[0];
  if (n > 1)
  {
    for (int i = 1; i < n; i++)
    {
      ma = fmax(a[i], ma);
    }
  }
  return ma;
}

// n --- number of random numbers
// k --- number of discrete values
// lw --- log of probabilities
// rn --- vector of random numbers returned
void sample_disc (const int n, int *rn, const int k, double *lw)
{
  // constructing probabilities from log probabilities
  double max_lw = fmaxm(k, lw);
  double cw[k];
  cw[0] = exp(lw[0] - max_lw);
  for (int i = 1; i < k; i++) 
    cw[i] = cw [i-1] + exp(lw [i] - max_lw);

  for (int j = 0; j < n; j++)
  {
    GetRNGstate();
    double u = unif_rand() * cw[k - 1];
    PutRNGstate();
    // convert u into a discrete value
    for (int i = 0; i < k; i++)
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
   (const int n, double *rn, const double lower, const double upper,
    double dlogf[1], double tol_dlogf_is0[1])
{
  // set smallest value for derivative that can be thought of as 0
  double y;
  int type_lin; 
  bool isfault = false;

  // checking linear function type and fault
  if (fabs (dlogf[0]) <= tol_dlogf_is0 [0])
  {
     if (!(std::isfinite (lower) & std::isfinite (upper)) )
     {
       isfault = true;
     }
     else
     {
       type_lin = 0;
     }
  }

  if (dlogf[0] >  tol_dlogf_is0 [0])
  {
    if (!std::isfinite(upper))
    {
      isfault = true;
    }
    else
    {
      type_lin = 1;
    }
  }

  if (dlogf[0] < -tol_dlogf_is0 [0])
  {
    if(!std::isfinite(lower))
    {
      isfault = true;
    }
    else
    {
      type_lin = -1;
    }
  }

  if (isfault)
  {
    printf("Error: in C function 'sample_elin':\n");
    printf("the exp linear function integrates to NAN/INFINITY\n");
    printf("(dlogf = %4.2f, lowerbound = %4.2f, upperbound = %4.2f)\n",
           dlogf[0], lower, upper);
    exit(1);
  }

  double dx = upper - lower;

  for (int j = 0; j < n; j++)
  {
    GetRNGstate();
    y = runif(0, 1);
    PutRNGstate();

    //converting uniform random number
    if (type_lin == 0)  // slope is 0
      rn[j] = lower + y * dx;

    if (type_lin == 1)  // slope is postive
      rn[j] = upper + log((1 - y) * exp(-dlogf[0] * dx) + y) / dlogf[0];

    if (type_lin == -1) //slope is negative
      rn[j] = lower + log(1 - y + y * exp(dlogf[0] * dx)) / dlogf[0];
  }
}

// this function evaluates the log of integral of exp linear hull
// logf --- value of linear hull at t
// dlogf --- value of derive of linear hull
// t --- tangent point where logf is calculated
// lower and upper --- lower and upper bounds of linear hull
// lw --- saved returned log integral
void logint_elin(double logf[1], double dlogf[1], double t[1],
                 double lower[1], double upper[1], double lw[1],
                 double tol_dlogf_is0[1])
{
  double dx, abs_dlogf;
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
double interc(
    double t1[1], double t2[1],
    double logf1[1], double logf2[1], double dlogf1[1], double dlogf2[1],
    double tol_ddlogf_is0[1])
{
  if (fabs(dlogf1[0] - dlogf2[0]) > tol_ddlogf_is0[0])
  {
    return (
        (logf2[0] - logf1[0] - dlogf2[0] * t2[0] + dlogf1[0] * t1[0]) /
        (dlogf1[0] - dlogf2[0]));
  }
  else
  {
    return ((t1[0] + t2[0]) / 2.0);
  }
}
