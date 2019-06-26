#include "RcppArmadillo.h"
//using namespace arma;

arma::mat cpvec(const arma::mat &A)
{
  return A(arma::span::all, arma::span::all);
}

arma::vec col_sum(const arma::mat &A)
{
  return arma::sum(A, 1);
}

double log_sum_exp(const arma::vec &a)
{  
  double m = a.max();
  return log(accu(arma::exp(a - m))) + m;
}

arma::vec log_sum_exp(const arma::mat &A)
{  
  arma::colvec m = arma::max(A, 1);
  return 
    arma::log(
      col_sum(
        arma::exp(
          A.each_row() - m
        )
      )
    ) + m;
}

arma::mat find_normlv(const arma::mat &lv)
{  
  return lv - log_sum_exp(lv);
}
