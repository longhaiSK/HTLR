#include "ars.h"

// NOTE: in this program, I call a piece of linear function as 'hull', and
// the piecewise linear function above logf upperhulls, and
// the piecewise linear function below logf lowerhulls

// Header of helper functions
int sample_disc(const int k, const double *lw);
double sample_elin(const double lower, const double upper,
                   const double dlogf, const double tol_dlogf_is0_);
double logint_elin(const double logf, const double dlogf, const double t,
                   const double lower, const double upper, const double tol_dlogf_is0_);
double interc(const double t1, const double t2,
              const double logf1, const double logf2,
              const double dlogf1, const double dlogf2,
              const double tol_ddlogf_is0_);

// this function updates the envolop and squeezing functions.
// newx --- new point to be inserted
// h --- index of the hull where newx is from
// logfv, dlogfv --- values of logf and dlogv at newx
void ARS::update_hulls(const int h, const double newx, const double logfv, const double dlogfv)
{
  int lh, rh, nh;

  if (no_hulls_ == max_nhull_) return;// reaching the limit of working vector

  //specify left and right hulls of new hull
  if (newx > tpoints_[h]) // to insert to the right of hull h
  {
    lh = h; rh = ritehulls_[h];
    // if logfv is -infinity, only update the rightest hull rightbound and lw
    if ((rh == max_nhull_) & (logfv == R_NegInf))
    {
      if (upperbounds_[h] != newx)
      {
        upperbounds_[h] = newx;
        lws_[h] = logint_elin(logfvs_[h], dlogfvs_[h], tpoints_[h],
                             lowerbounds_[h], upperbounds_[h], tol_dlogf_is0_);
      }
      return;
    }
  }
  else // to insert to the left of hull h
  {
    lh = lefthulls_[h]; rh = h;
    // if logfv is -infinity, only update the leftest hull leftbound and lw
    if ((lh == -1) & (logfv == R_NegInf))
    {
      if (lowerbounds_[h] != newx)
      {
        lowerbounds_[h] = newx;
        lws_[h] = logint_elin(logfvs_[h], dlogfvs_[h], tpoints_[h],
                             lowerbounds_[h], upperbounds_[h], tol_dlogf_is0_);
      }
      return;
    }
  }

  // insert a new hull
  nh = no_hulls_;
  no_hulls_ ++;
  tpoints_[nh] = newx;
  logfvs_[nh] = logfv;
  dlogfvs_[nh] = dlogfv;
  lefthulls_[nh] = lh;
  ritehulls_[nh] = rh;

  if (lh == -1) // nh will be the new leftest hull
  {
    lowerbounds_[nh] = lowerbounds_[h];
    slopes_leftsq_[nh] = R_PosInf;
  }
  else
  {
    lowerbounds_[nh] = interc(
        tpoints_[lh], tpoints_[nh], logfvs_[lh], logfvs_[nh],
        dlogfvs_[lh], dlogfvs_[nh], tol_ddlogf_is0_);
    slopes_leftsq_[nh] = (logfvs_[nh] - logfvs_[lh]) /
                        (tpoints_[nh] - tpoints_[lh]);
  }
  if (rh == max_nhull_)
  {
    upperbounds_[nh] = upperbounds_[h];
    slopes_ritesq_[nh] = R_NegInf;
  }
  else
  {
    upperbounds_[nh] =
        interc(tpoints_[nh], tpoints_[rh], logfvs_[nh], logfvs_[rh],
               dlogfvs_[nh], dlogfvs_[rh], tol_ddlogf_is0_);
    slopes_ritesq_[nh] = (logfvs_[nh] - logfvs_[rh]) /
                        (tpoints_[nh] - tpoints_[rh]);
  }

  lws_[nh] = logint_elin(logfvs_[nh], dlogfvs_[nh], tpoints_[nh],
                        lowerbounds_[nh], upperbounds_[nh], tol_dlogf_is0_);

  // update left hull of new null
  if (lh != -1)
  {
    upperbounds_[lh] = lowerbounds_[nh];
    ritehulls_[lh] = nh;
    slopes_ritesq_[lh] = slopes_leftsq_[nh];
    lws_[lh] = logint_elin(logfvs_[lh], dlogfvs_[lh], tpoints_[lh],
                          lowerbounds_[lh], upperbounds_[lh], tol_dlogf_is0_);
  }

  // update right hull of newh if it exists
  if (rh != max_nhull_)
  {
    lowerbounds_[rh] = upperbounds_[nh];
    lefthulls_[rh] = nh;
    slopes_leftsq_[rh] = slopes_ritesq_[nh];

    lws_[rh] = logint_elin(logfvs_[rh], dlogfvs_[rh], tpoints_[rh],
                          lowerbounds_[rh], upperbounds_[rh], tol_dlogf_is0_);
  }
}

double ARS::eval_upperhull(const int h, const double newx)
{
  return ((newx - tpoints_[h]) * dlogfvs_[h] + logfvs_[h]);
}

double ARS::eval_lowerhull(const int h, const double newx)
{
    if (newx >= tpoints_[h])
    {
      return ((newx - tpoints_[h]) * slopes_ritesq_[h] + logfvs_[h]);
    }
    else
    {
      return ((newx - tpoints_[h]) * slopes_leftsq_[h] + logfvs_[h]);
    }
}

void ARS::Initialize()
{
  // if lb is finite, bound the first hull at left
  // or insert a hull tangent at lb if logf at lb is finite too
  if (R_FINITE(lb_))
  {
    h_ = 0;
    newx_ = lb_;
    target_->eval_logf(newx_, newlogf_, newdlogf_);
    update_hulls(h_, newx_, newlogf_, newdlogf_);
  }
  // expanding at the left until reaching a bound or integral to finite
  else
  {
    newx_ = tpoints_[0] - stepout_;
    do
    {
      if (no_hulls_ == max_nhull_)
      {
        Rcpp::stop(
            "Error in Rejection Sampling: (finite lb)\n"
            "'max_nhull_' is set too small, or your log-PDF is NOT concave.\n");
      }
      h_ = 0;
      target_->eval_logf(newx_, newlogf_, newdlogf_);
      update_hulls(h_, newx_, newlogf_, newdlogf_);
      // finding a new leftbound, quit expanding
      if (newlogf_ == R_NegInf) break;
      newx_ -= stepout_;
      h_ = no_hulls_ - 1;
    }
    while (newdlogf_ < tol_dlogf_is0_);
  }

  // if ub is finite, bound the first hull at the right
  // or insert a hull tangent at ub if logf at ub is finite too
  if (R_FINITE(ub_))
  {
    h_= 0;
    newx_ = ub_;
    target_->eval_logf(newx_, newlogf_, newdlogf_);
    update_hulls(h_, newx_, newlogf_, newdlogf_);
  }
  else // expanding at the right until reaching a bound or integral to finite
  {
    h_ = 0;
    newx_ = tpoints_[0] + stepout_;
    do
    {
      if (no_hulls_ == max_nhull_)
      {
        //Rcpp::Rcerr << no_hulls << " " << max_nhull << "\n";
        Rcpp::stop(
            "Error in Rejection Sampling: (finite ub)\n"
            "'max_nhull' is set too small, or your log-PDF is NOT concave.\n");
      }
      target_->eval_logf(newx_, newlogf_, newdlogf_);
      update_hulls(h_, newx_, newlogf_, newdlogf_);
      if (!R_FINITE(newlogf_)) break;
      newx_ += stepout_;
      h_ = no_hulls_ - 1;
    }
    while (newdlogf_ > - tol_dlogf_is0_);
  }
}

ARS::ARS(int n, SampleTarget *target_, double ini_tpoint, 
    double lb/*= -INFINITY*/, double ub/*= +INFINITY*/,   
    bool verbose/*=false*/, int max_nhull_/*=1000*/, double stepout_/*=10*/,
    double tol_dlogf_is0_/*= 1E-5*/, double tol_ddlogf_is0_/*= 1E-5*/)
    : n_(n), lb_(lb), ub_(ub), verbose_(verbose), max_nhull_(max_nhull_), stepout_(stepout_), 
    tol_dlogf_is0_(tol_dlogf_is0_), tol_ddlogf_is0_(tol_ddlogf_is0_), target_(target_)  
{
  // construct the first hull
  logfvs_  = new double[max_nhull_] {0};
  dlogfvs_ = new double[max_nhull_] {0};
  tpoints_ = new double[max_nhull_] {0};
  tpoints_[0] = ini_tpoint; // the tangent point
  target_->eval_logf(tpoints_[0], logfvs_[0], dlogfvs_[0]);
  if (!R_FINITE(logfvs_[0]))
  {
    Rcpp::stop(
        "Error in adaptive rejection sampling:\n"
        "the first tangent point doesn't have positive probability.\n");
  }

  lowerbounds_ = new double[max_nhull_] {0};
  upperbounds_ = new double[max_nhull_] {0};
  lowerbounds_[0] = fmax(lb, R_NegInf); // lower bound of the hull
  upperbounds_[0] = fmin(ub, R_PosInf); // upper bound of the hull

  lefthulls_ = new int[max_nhull_] {0};
  ritehulls_ = new int[max_nhull_] {0};
  lefthulls_[0] = -1;        // index of left hull
  ritehulls_[0] = max_nhull_; // index of right hull

  slopes_leftsq_ = new double[max_nhull_] {0};
  slopes_ritesq_ = new double[max_nhull_] {0};
  slopes_leftsq_[0] = R_PosInf; // slope of left squeezing arc
  slopes_ritesq_[0] = R_NegInf; // slope of right sequeezing arc

  lws_ = new double[max_nhull_] {0};
  lws_[0] = R_PosInf; // compute log weights, updating lws_[0]
  no_hulls_ = 1;
}

ARS::~ARS()
{
  delete[] logfvs_;
  delete[] dlogfvs_;
  delete[] tpoints_;
  delete[] lowerbounds_;
  delete[] upperbounds_;
  delete[] lefthulls_;
  delete[] ritehulls_;
  delete[] slopes_leftsq_;
  delete[] slopes_ritesq_;
  delete[] lws_;
}

// Do adaptive rejection sampling
Rcpp::NumericVector ARS::Sample()
{
  Initialize();

  /* define parameters used while sampling */ 
  double
      upperhullv, // value of upper hull at newx
      lowerhullv, // value of lower (squeezing) hull at newx
      u,          // a random number used to determine whether to accept
      logacceptv; // if logacceptv is smaller than logf, newx accepted
  int no_rejs = 0;
  Rcpp::NumericVector output (n_); // sampling output

  for (int i = 0; i < n_; i++)
  {
    bool rejected = true;
    while (rejected)
    {
      // draw a new point and a unif random number
      h_ = sample_disc(no_hulls_, lws_);
      newx_ = sample_elin(lowerbounds_[h_], upperbounds_[h_],
                         dlogfvs_[h_], tol_dlogf_is0_);
      upperhullv = eval_upperhull(h_, newx_);

      GetRNGstate();
      u = unif_rand();
      PutRNGstate();

      logacceptv = upperhullv + log(u);
      lowerhullv = eval_lowerhull(h_, newx_);
      //check acceptance with squeezing function
      if (logacceptv <= lowerhullv)
      {
        output[i] = newx_;
        rejected = false;
      }
      else
      {
        // check acceptance with logf
        // eval logf at newx and insert a new hull
        target_->eval_logf(newx_, newlogf_, newdlogf_);
        update_hulls(h_, newx_, newlogf_, newdlogf_);
        if (logacceptv <= newlogf_)
        {
          output[i] = newx_;
          rejected = false;
        }
        else
          no_rejs++;
      }
    }
  }
  if (verbose_)
  {
    double rej_rate = (no_rejs + 0.0) / (no_rejs + n_ + 0.0);
    Rprintf("Sampling complete. Number of hulls: %d, Rejection rate: %4.2f\n",
            no_hulls_, rej_rate);
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
  double *cw = new double[k];
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
  delete[] cw;
  return i;
}

// this function samples one point from: exp (a[0]*x) I (x in [lb, upper[0]])
double sample_elin(const double lower, const double upper,
                   const double dlogf, const double tol_dlogf_is0_)
{
  // set smallest value for derivative that can be thought of as 0
  int type_lin = -1; 
  bool isfault = false;

  // checking linear function type and fault
  if (fabs(dlogf) <= tol_dlogf_is0_)
  {
    if (!(R_FINITE(lower) & R_FINITE(upper)))
      isfault = true;
    else
      type_lin = 0; // slope is zero
  }

  if (dlogf > tol_dlogf_is0_)
  {
    if (!R_FINITE(upper))
      isfault = true;
    else
      type_lin = 1; // slope is postive
  }

  if (dlogf < -tol_dlogf_is0_)
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
                   const double lower, const double upper, const double tol_dlogf_is0_)
{
  double output;

  double dx = upper - lower;
  double abs_dlogf = fabs(dlogf);

  if (abs_dlogf <= tol_dlogf_is0_) // slope is 0
  {
    output = logf + log(dx);
  }
  else
  {
    if (dlogf > tol_dlogf_is0_) // slope is positive
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
              const double tol_ddlogf_is0_)
{
  if (fabs(dlogf1 - dlogf2) > tol_ddlogf_is0_)
    return ((logf2 - logf1 - dlogf2 * t2 + dlogf1 * t1) / (dlogf1 - dlogf2));
  else
    return ((t1 + t2) / 2.0);
}
