#ifndef ARS_H
#define ARS_H 

class SampleTarget 
{ 
  public: 
  virtual void eval_logf(const double x, double &logf, double &dlogf); 
};

class ARS
{
  private:

  const int n; // number of samples
  double *rn;  // sampling output
  
  const double lb, ub;     // lower and upper bounds of logf
  const double ini_tpoint; // initial tangent point
  const bool verbose;      // flag whether to print ARS information

  const int max_nhull; // maximum number of pieces of linear hulls
  double stepout = 10; // size of stepout in initializing linear hulls

  // the smallest difference of derivatives that can be thought of as the same
  // so that we don't need to insert one more hull
  double tol_dlogf_is0, tol_ddlogf_is0;
  
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
  SampleTarget target;

  void update_hulls(const int h, const double newx, const double logfv, const double dlogfv);
  double eval_upperhull(const int h, const double newx);
  void Prepare();

  public:

  ARS(int n, SampleTarget target, double ini_tpoint,
      double lb = -INFINITY, double ub = +INFINITY,
      bool verbose = false, int max_nhull = 1000,
      double tol_dlogf_is0 = 1E-5, double tol_ddlogf_is0 = 1E-5, double stepout = 10);

  double* Sampling(); 

};

#endif 
