#ifndef UTILS_H
#define UTILS_H 

#include "RcppArmadillo.h"

arma::mat cpvec(const arma::mat &A);
arma::vec col_sum(const arma::mat &A);
arma::vec row_sum(const arma::mat &A);
double log_sum_exp(const arma::vec &a);
arma::vec log_sum_exp(const arma::mat &A);
arma::mat find_normlv(const arma::mat &lv);

#endif 
