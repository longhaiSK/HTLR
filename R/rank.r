#' Rank features by Kruskal-Wallis test
#'
#' This function ranks all features in terms of Kruskal-Wallis test p-value.
#' 
#' @param X Input matrix, of dimension \code{nobs} by \code{nvars}; each row is an observation vector.
#' 
#' @param y Vector of response variables.
#' 
#' @return Rank of all features of length \code{nvars}.
#' 
#' @export
#'
#' @examples 
#' data("diabetes392")
#' rank_kruskal(diabetes392$X, diabetes392$y)
#' 
rank_kruskal <- function(X, y)
{
  pv_kt <- function(x, g) { kruskal.test(x, g)$p.value }
  order(apply(X, 2, pv_kt, g = y))
}

# Helper function for rank_ftest() 
panova <- function(X, y)
{
  ## This function computes the values of F-statistic of all the features.
  n <- length (y)
  nos_g <- as.vector (tapply (rep (1, n), INDEX = y, sum))
  G <- length (nos_g)
  
  gsum_X <- rowsum (X, y)
  sum_X2 <- colSums (X ^ 2)
  sum_gsum2 <- colSums (gsum_X ^ 2 / nos_g)
  
  pvars <- (sum_X2 - sum_gsum2) / (n - G) ## pooled variances
  
  sum_X <- colSums (X)
  gvars <-
    (sum_gsum2 - sum_X ^ 2 / n) / (G - 1) ## variances btw groups
  ## F-statistic
  fstats <- gvars / pvars
  
  pvalues <- 1 - pf (fstats,
                     df1 = G - 1,
                     df2 = n - G,
                     ncp = 0)
  pvalues
}

#' Rank features by F-statistic
#'
#' This function ranks all features in terms of ANOVA F-statistic.
#' 
#' @param X Input matrix, of dimension \code{nobs} by \code{nvars}; each row is an observation vector.
#' 
#' @param y Vector of response variables.
#' 
#' @return Rank of all features of length \code{nvars}.
#' 
#' @export
#'
#' @examples 
#' data("diabetes392")
#' rank_ftest(diabetes392$X, diabetes392$y)
#' 
rank_ftest <- function(X, y) { order(panova(X, y)) }

#' Plain rank function
#'
#' A placeholder rank function that returns the original order of given features.
#' 
#' @param X Input matrix, of dimension \code{nobs} by \code{nvars}; each row is an observation vector.
#' 
#' @param y Vector of response variables.
#' 
#' @return Sequence starting from 1 to \code{nvars}.
#' 
#' @export
#' 
#' @keywords internal
#' 
rank_plain <- function(X, y) { 1L:ncol(X) }

