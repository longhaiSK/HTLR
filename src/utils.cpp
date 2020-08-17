#include "utils.h"

double log_sum_exp(const arma::vec &a)
{  
  double m = a.max();
  return log(accu(arma::exp(a - m))) + m;
}

// [[Rcpp::export]]
arma::vec log_sum_exp(const arma::mat &A)
{
  arma::colvec m = arma::max(A, 1);
  return 
    arma::log(
      row_sum(
        arma::exp(
          A.each_col() - m
        )
      )
    ) + m;
}

arma::mat find_normlv(const arma::mat &lv)
{  
  return lv.each_col() - log_sum_exp(lv);
}

// [[Rcpp::export]]
arma::vec spl_sgm_ig(double alpha, int K, double w, const arma::vec &vardeltas)
{
  arma::vec rn_gamma = Rcpp::rgamma(vardeltas.n_elem, (alpha + K) / 2); 
  return (1.0 / rn_gamma % (alpha * w + vardeltas) / 2.0);
}

// [[Rcpp::export]]
Rcpp::List std_helper(const arma::mat &A)
{
  arma::rowvec nuj = arma::median(A);
  arma::rowvec sdj = arma::stddev(A);
  arma::mat A_std = A.each_row() - nuj;
  A_std.each_row() /= sdj;

  return Rcpp::List::create(
    Rcpp::Named("median") = nuj,
    Rcpp::Named("sd") = sdj,
    Rcpp::Named("X") = A_std
  );
}

// [[Rcpp::export]]
arma::vec comp_vardeltas(const arma::mat &deltas)
{
  arma::vec sum_deltas = row_sum(deltas);
  arma::vec sum_sq_deltas = row_sum(arma::square(deltas));
  return sum_sq_deltas - (arma::square(sum_deltas) / (deltas.n_cols + 1));
}

// [[Rcpp::export]]
arma::vec comp_lsl(arma::mat &A)
{
  A.insert_cols(0, 1);
  return log_sum_exp(A);
}

// [[Rcpp::export]]
double log_normcons(arma::mat &A)
{
  return arma::accu(comp_lsl(A));
}

// [[Rcpp::export]]
Rcpp::List gendata_FAM_helper(int n, arma::mat &muj, const arma::mat &muj_rep, const arma::mat &A, double sd_g, bool stdx)
{
  int p = muj.n_rows;
  // int c = muj.n_cols;
  int k = A.n_cols;

  arma::vec rn_nk = Rcpp::rnorm(n * k);
  arma::mat KN = arma::reshape(rn_nk, k, n); 
  arma::mat X = A * KN + muj_rep; 

  arma::vec rn_np = Rcpp::rnorm(n * p); 
  arma::mat NP = arma::reshape(rn_np, arma::size(X));
  X += NP * sd_g; 

  arma::mat SGM = A * A.t() + arma::eye(p, p) * R_pow_di(sd_g, 2);
  
  if (stdx)
  {
    arma::vec mux = arma::mean(muj, 1);
    arma::vec sdx = arma::sqrt(SGM.diag() + arma::var(muj, 1, 1));
    muj.each_col() -= mux;
    muj.each_col() /= sdx;
    X.each_col() -= mux;
    X.each_col() /= sdx;
    SGM.each_col() /= sdx;
    SGM.each_row() /= sdx.t();
  }

  return Rcpp::List::create(
    Rcpp::Named("X") = X.t(),
    Rcpp::Named("muj") = muj,
    Rcpp::Named("SGM") = SGM
  );
}
