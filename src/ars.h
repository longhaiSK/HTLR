#ifndef ARS_H
#define ARS_H 

#include <Rcpp.h>

// An interface of sampling target, to be extended by users.
class SampleTarget 
{ 
  public: 
  virtual void eval_logf(const double x, double &logf, double &dlogf) = 0; 
};

// Main ARS class
class ARS
{
  private:

  const int n; // number of samples
  
  const double lb, ub;     // lower and upper bounds of logf
  //const double ini_tpoint; // initial tangent point
  const bool verbose;      // flag whether to print ARS information

  /* optional User Settings */

  const int max_nhull; // maximum number of pieces of linear hulls
  double stepout; // size of stepout in initializing linear hulls

  // the smallest difference of derivatives that can be thought of as the same
  // so that we don't need to insert one more hull
  double tol_dlogf_is0, tol_ddlogf_is0;

  /* global working vectors used to represent linear hulls */
  
  double 
    *tpoints,                       // tangent points 
    *lws,                           // log integrals of exp hulls
    *lowerbounds, *upperbounds,     // bounds of hulls 
    *logfvs, *dlogfvs,              // values of logf and dlogf at tangent points
    *slopes_leftsq, *slopes_ritesq; // slopes of left and right squeezings

  int *lefthulls, *ritehulls;     // indice of left and right hulls
  int no_hulls;                   // total number of hulls 

  double newx, newlogf, newdlogf; // a new tangent point to be inserted
  int h;                          // index of the hull where newx is from

  // the sampling target object which provides function for evaluating logf and dlogf  
  SampleTarget *target;

  void update_hulls(const int h, const double newx, const double logfv, const double dlogfv);
  double eval_upperhull(const int h, const double newx);
  double eval_lowerhull(const int h, const double newx);
  void Prepare();

  public:

  ARS(int n, SampleTarget *target, double ini_tpoint,
      double lb = R_NegInf, double ub = R_PosInf,
      bool verbose = false, int max_nhull = 1000,
      double tol_dlogf_is0 = 1E-5, double tol_ddlogf_is0 = 1E-5, double stepout = 10);

  Rcpp::NumericVector Sample(); 

};

#endif 
