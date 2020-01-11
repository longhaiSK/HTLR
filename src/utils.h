#ifndef UTILS_H
#define UTILS_H 

#include "RcppArmadillo.h"

double log_sum_exp(const arma::vec &a);
arma::vec log_sum_exp(const arma::mat &A);
arma::mat find_normlv(const arma::mat &lv);
arma::vec spl_sgm_ig(double alpha, int K, double w, const arma::vec &vardeltas);

inline arma::rowvec col_sum(const arma::mat &A)
{
  return arma::sum(A, 0);
}

inline arma::vec row_sum(const arma::mat &A)
{
  return arma::sum(A, 1);
}

#endif 
