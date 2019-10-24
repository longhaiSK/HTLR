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
