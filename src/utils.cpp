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
    arma::vec sum_sq_deltas = row_sum(arma::pow(deltas,2));
    return sum_sq_deltas - (arma::pow(sum_deltas, 2) / (deltas.n_cols + 1));
}