#include "ars.h"

// NOTE: in this program, I call a piece of linear function as 'hull', and
// the piecewise linear function above logf upperhulls, and
// the piecewise linear function below logf lowerhulls

// Header of helper functions
int sample_disc(const int k, const double *lw);
double sample_elin(const double lower, const double upper,
                   const double dlogf, const double tol_dlogf_is0);
double logint_elin(const double logf, const double dlogf, const double t,
                   const double lower, const double upper, const double tol_dlogf_is0);
double interc(const double t1, const double t2,
              const double logf1, const double logf2,
              const double dlogf1, const double dlogf2,
              const double tol_ddlogf_is0);

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
    if ((rh == max_nhull) & (logfv == R_NegInf))
    {
      if (upperbounds[h] != newx)
      {
        upperbounds[h] = newx;
        lws[h] = logint_elin(logfvs[h], dlogfvs[h], tpoints[h],
                             lowerbounds[h], upperbounds[h], tol_dlogf_is0);
      }
      return;
    }
  }
  else // to insert to the left of hull h
  {
    lh = lefthulls[h]; rh = h;
    // if logfv is -infinity, only update the leftest hull leftbound and lw
    if ((lh == -1) & (logfv == R_NegInf))
    {
      if (lowerbounds[h] != newx)
      {
        lowerbounds[h] = newx;
        lws[h] = logint_elin(logfvs[h], dlogfvs[h], tpoints[h],
                             lowerbounds[h], upperbounds[h], tol_dlogf_is0);
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
    lowerbounds[nh] = lowerbounds[h];
    slopes_leftsq[nh] = R_PosInf;
  }
  else
  {
    lowerbounds[nh] = interc(
        tpoints[lh], tpoints[nh], logfvs[lh], logfvs[nh],
        dlogfvs[lh], dlogfvs[nh], tol_ddlogf_is0);
    slopes_leftsq[nh] = (logfvs[nh] - logfvs[lh]) /
                        (tpoints[nh] - tpoints[lh]);
  }
  if (rh == max_nhull)
  {
    upperbounds[nh] = upperbounds[h];
    slopes_ritesq[nh] = R_NegInf;
  }
  else
  {
    upperbounds[nh] =
        interc(tpoints[nh], tpoints[rh], logfvs[nh], logfvs[rh],
               dlogfvs[nh], dlogfvs[rh], tol_ddlogf_is0);
    slopes_ritesq[nh] = (logfvs[nh] - logfvs[rh]) /
                        (tpoints[nh] - tpoints[rh]);
  }

  lws[nh] = logint_elin(logfvs[nh], dlogfvs[nh], tpoints[nh],
                        lowerbounds[nh], upperbounds[nh], tol_dlogf_is0);

  // update left hull of new null
  if (lh != -1)
  {
    upperbounds[lh] = lowerbounds[nh];
    ritehulls[lh] = nh;
    slopes_ritesq[lh] = slopes_leftsq[nh];
    lws[lh] = logint_elin(logfvs[lh], dlogfvs[lh], tpoints[lh],
                          lowerbounds[lh], upperbounds[lh], tol_dlogf_is0);
  }

  // update right hull of newh if it exists
  if (rh != max_nhull)
  {
    lowerbounds[rh] = upperbounds[nh];
    lefthulls[rh] = nh;
    slopes_leftsq[rh] = slopes_ritesq[nh];

    lws[rh] = logint_elin(logfvs[rh], dlogfvs[rh], tpoints[rh],
                          lowerbounds[rh], upperbounds[rh], tol_dlogf_is0);
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
  if (R_FINITE(lb))
  {
    h = 0;
    newx = lb;
    target->eval_logf(newx, newlogf, newdlogf);
    update_hulls(h, newx, newlogf, newdlogf);
  }
  // expanding at the left until reaching a bound or integral to finite
  else
  {
    newx = tpoints[0] - stepout;
    do
    {
      if (no_hulls == max_nhull)
      {
        Rcpp::stop(
            "Error in Rejection Sampling: (finite lb)\n"
            "'max_nhull' is set too small, or your log-PDF is NOT concave.\n");
      }
      h = 0;
      target->eval_logf(newx, newlogf, newdlogf);
      update_hulls (h, newx, newlogf, newdlogf);
      // finding a new leftbound, quit expanding
      if (newlogf == R_NegInf) break;
      newx -= stepout;
      h = no_hulls - 1;
    }
    while (newdlogf < tol_dlogf_is0);
  }

  // if ub is finite, bound the first hull at the right
  // or insert a hull tangent at ub if logf at ub is finite too
  if (R_FINITE(ub))
  {
    h = 0;
    newx = ub;
    target->eval_logf(newx, newlogf, newdlogf);
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
        //Rcpp::Rcerr << no_hulls << " " << max_nhull << "\n";
        Rcpp::stop(
            "Error in Rejection Sampling: (finite ub)\n"
            "'max_nhull' is set too small, or your log-PDF is NOT concave.\n");
      }
      target->eval_logf(newx, newlogf, newdlogf);
      update_hulls (h, newx, newlogf, newdlogf);
      if (!R_FINITE(newlogf)) break;
      newx += stepout;
      h = no_hulls - 1;
    }
    while (newdlogf > - tol_dlogf_is0);
  }
}

ARS::ARS(int n, SampleTarget *target, double ini_tpoint, 
    double lb/*= -INFINITY*/, double ub/*= +INFINITY*/,   
    bool verbose/*=false*/, int max_nhull/*=1000*/, double stepout/*=10*/,
    double tol_dlogf_is0/*= 1E-5*/, double tol_ddlogf_is0/*= 1E-5*/)
    : n(n), lb(lb), ub(ub), verbose(verbose), max_nhull(max_nhull), stepout(stepout), 
    tol_dlogf_is0(tol_dlogf_is0), tol_ddlogf_is0(tol_ddlogf_is0), target(target)  
{
  // construct the first hull
  this->logfvs  = new double[max_nhull] {0};
  this->dlogfvs = new double[max_nhull] {0};
  this->tpoints = new double[max_nhull] {0};
  this->tpoints[0] = ini_tpoint; // the tangent point
  this->target->eval_logf(tpoints[0], logfvs[0], dlogfvs[0]);
  if (!R_FINITE(logfvs[0]))
  {
    Rcpp::stop(
        "Error in adaptive rejection sampling:\n"
        "the first tangent point doesn't have positive probability.\n");
  }

  this->lowerbounds = new double[max_nhull] {0};
  this->upperbounds = new double[max_nhull] {0};
  lowerbounds[0] = fmax(lb, R_NegInf); // lower bound of the hull
  upperbounds[0] = fmin(ub, R_PosInf); // upper bound of the hull

  this->lefthulls = new int[max_nhull]{0};
  this->ritehulls = new int[max_nhull]{0};
  lefthulls[0] = -1;        // index of left hull
  ritehulls[0] = max_nhull; // index of right hull

  this->slopes_leftsq = new double[max_nhull] {0};
  this->slopes_ritesq = new double[max_nhull] {0};
  slopes_leftsq[0] = R_PosInf; // slope of left squeezing arc
  slopes_ritesq[0] = R_NegInf; // slope of right sequeezing arc

  this->lws = new double[max_nhull] {0};
  lws[0] = R_PosInf; // compute log weights, updating lws[0]
  this->no_hulls = 1;
}

// Do adaptive rejection sampling
Rcpp::NumericVector ARS::Sample()
{
  Prepare();
  // define parameters used while sampling
  int one = 1, no_rejs = 0;
  double
      upperhullv, // value of upper hull at newx
      lowerhullv, // value of lower (squeezing) hull at newx
      u,          // a random number used to determine whether to accept
      logacceptv; // if logacceptv is smaller than logf, newx accepted
  Rcpp::NumericVector output (n); // sampling output

  for (int i = 0; i < n; i++)
  {
    bool rejected = true;
    while (rejected)
    {
      // draw a new point and a unif random number
      h = sample_disc(no_hulls, lws);
      newx = sample_elin(lowerbounds[h], upperbounds[h],
                         dlogfvs[h], tol_dlogf_is0);
      upperhullv = eval_upperhull(h, newx);
      GetRNGstate();
      u = unif_rand();
      PutRNGstate();
      logacceptv = upperhullv + log(u);
      lowerhullv = eval_lowerhull(h, newx);
      //check acceptance with squeezing function
      if (logacceptv <= lowerhullv)
      {
        output[i] = newx;
        rejected = false;
      }
      else
      {
        // check acceptance with logf
        // eval logf at newx and insert a new hull
        target->eval_logf(newx, newlogf, newdlogf);
        update_hulls(h, newx, newlogf, newdlogf);
        if (logacceptv <= newlogf)
        {
          output[i] = newx;
          rejected = false;
        }
        else
          no_rejs++;
      }
    }
  }
  if (verbose)
  {
    double rej_rate = (no_rejs + 0.0) / (no_rejs + n + 0.0);
    Rprintf("Sampling complete. Number of hulls: %d, Rejection rate: %4.2f\n",
            no_hulls, rej_rate);
  }
  return output;
}

// find maximum value in double array a with length n
double fmaxm(const int n, const double *a)
{
  double ma = a[0];
  if (n > 1)
  {
    for (int i = 1; i < n; i++)
      ma = fmax(a[i], ma);
  }
  return ma;
}

// k --- number of discrete values
// lw --- log of probabilities
int sample_disc(const int k, const double *lw)
{
  // constructing probabilities from log probabilities
  double max_lw = fmaxm(k, lw);
  double cw[k];
  cw[0] = exp(lw[0] - max_lw);
  for (int i = 1; i < k; i++)
    cw[i] = cw[i - 1] + exp(lw[i] - max_lw);

  GetRNGstate();
  double u = unif_rand() * cw[k - 1];
  PutRNGstate();
  // convert u into a discrete value
  int i = 0;
  while (i < k)
  {
    if (u <= cw[i])
      break;
    i++;
  }
  return i;
}

// this function samples one point from: exp (a[0]*x) I (x in [lb, upper[0]])
double sample_elin(const double lower, const double upper,
                   const double dlogf, const double tol_dlogf_is0)
{
  // set smallest value for derivative that can be thought of as 0
  int type_lin = -1; 
  bool isfault = false;

  // checking linear function type and fault
  if (fabs(dlogf) <= tol_dlogf_is0)
  {
    if (!(R_FINITE(lower) & R_FINITE(upper)))
      isfault = true;
    else
      type_lin = 0; // slope is zero
  }

  if (dlogf > tol_dlogf_is0)
  {
    if (!R_FINITE(upper))
      isfault = true;
    else
      type_lin = 1; // slope is postive
  }

  if (dlogf < -tol_dlogf_is0)
  {
    if (!R_FINITE(lower))
      isfault = true;
    else
      type_lin = 2; //slope is negative
  }

  if (isfault)
  {
    REprintf("(dlogf = %4.2f, lowerbound = %4.2f, upperbound = %4.2f)\n",
             dlogf, lower, upper);
    Rcpp::stop(
        "Error: in C function 'sample_elin':\n"
        "the exp linear function integrates to NAN/INFINITY\n");
  }

  double dx = upper - lower;
  GetRNGstate();
  double y = R::runif(0, 1);
  PutRNGstate();

  double output;

  if (type_lin == 0)       
    output = lower + y * dx;
  else if (type_lin == 1)  
    output = upper + log((1 - y) * exp(-dlogf * dx) + y) / dlogf;
  else if (type_lin == 2)  
    output = lower + log(1 - y + y * exp(dlogf * dx)) / dlogf;
  else
    Rcpp::stop("Error: in C function 'sample_elin': unexpected type_lin value\n");

  return output;
}

// this function evaluates the log of integral of exp linear hull
// logf --- value of linear hull at t
// dlogf --- value of derive of linear hull
// t --- tangent point where logf is calculated
// lower and upper --- lower and upper bounds of linear hull
double logint_elin(const double logf, const double dlogf, const double t,
                   const double lower, const double upper, const double tol_dlogf_is0)
{
  double output;

  double dx = upper - lower;
  double abs_dlogf = fabs(dlogf);

  if (abs_dlogf <= tol_dlogf_is0) // slope is 0
  {
    output = logf + log(dx);
  }
  else
  {
    if (dlogf > tol_dlogf_is0) // slope is positive
    {
      output = logf + dlogf * (upper - t) - log(abs_dlogf) +
               log(1 - exp(-abs_dlogf * dx));
    }
    else //slope is negative
    {
      output = logf + dlogf * (lower - t) - log(abs_dlogf) +
               log(1 - exp(-abs_dlogf * dx));
    }
  }
  return output;
}

// this function finds interception points between t1 and t2
double interc(const double t1, const double t2,
              const double logf1, const double logf2,
              const double dlogf1, const double dlogf2,
              const double tol_ddlogf_is0)
{
  if (fabs(dlogf1 - dlogf2) > tol_ddlogf_is0)
    return ((logf2 - logf1 - dlogf2 * t2 + dlogf1 * t1) / (dlogf1 - dlogf2));
  else
    return ((t1 + t2) / 2.0);
}
