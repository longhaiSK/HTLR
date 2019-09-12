#' Generate simulated data with multinomial logistic regression model
#' 
#' This function generates the response variables \code{y} given 
#' optional supplied \code{X} using a multinomial logistic regression model.
#' 
#' @param n Number of observations.
#' @param p Number of features.
#' @param NC Number of classes for response variables.
#' @param nu,w If \code{betas} is not supplied (default), the regression coefficients are generated with 
#' t prior with df = \code{nu}, scale = \code{sqrt(w)}; will be ignored if \code{betas} is supplied.
#' @param X The design matrix; will be generated from standard normal distribution if not supplied.
#' @param betas User supplied regression coefficients.  
#' 
#' @return A list contains input matrix \code{X}, response variables \code{y}, and regression coefficients \code{deltas}.      
#'
#' @export
#' 
#' @examples
#' set.seed(12345)
#' dat <- gendata_MLR(n = 100, p = 10)
#' ggplot2::qplot(dat$y, bins = 6)
#' corrplot::corrplot(cor(dat$X))
#' 
#' @seealso \code{\link{gendata_FAM}}
#' 
gendata_MLR <- function(n, p, NC = 3, nu = 2, w = 1, X = NULL, betas = NULL)
{
  if (is.null(X)) {
    X <- matrix(rnorm (n * p), n, p)
    colnames(X) <- paste0("V", 1:p)
  }

  if (is.null(betas))
  {
    sigmasbt <- 1 / rgamma(p, nu / 2, nu / 2) * w
    betas <- replicate(NC, rnorm(p + 1)) * c(0, sqrt(sigmasbt))
  }
  
  deltas <- betas[, -1, drop = FALSE] - betas[, 1]  
  lv <- cbind(1, X) %*% betas
  probs <- exp(lv - as.vector(log_sum_exp(lv)))

  y <- apply(probs, 1, function(prob){sample(1:NC, 1, TRUE, prob)})

  list("X" = X, "y" = y, "deltas" = deltas)
}

#' Generate simulated data with factor analysis model 
#' 
#' This function generates inputs \code{X} given by the response variable \code{y} 
#' using a multivariate normal model. 
#' 
#' The means of each covariate \eqn{x_j} depend on \code{y} specified by the 
#' matrix \code{muj}; the covariate matrix \eqn{\Sigma} of the multivariate normal 
#' is equal to \eqn{AA^t\delta^2I}, where \code{A} is the factor loading matrix
#' and \eqn{\delta} is the noise level.  
#' 
#' @param n Number of observations.
#' @param muj C by p matrix, with row c representing y = c, and column j representing \eqn{x_j}.
#' Used to specify \code{y}. 
#' @param A The factor loading matrix, see details.
#' @param sd_g Noise level \eqn{\delta}, see details.
#' @param stdx Logical; if \code{TRUE}, data \code{X} is standardized to have \code{mean = 0} and \code{sd = 1}.
#'
#' @return A list contains input matrix \code{X}, response variables \code{y},
#' covariate matrix \code{SGM} and \code{muj} (standardized if \code{stdx = TRUE}).
#' 
#' @export
#' 
#' @examples 
#' ## feature #1: marginally related feature
#' ## feature #2: marginally unrelated feature, but feature #2 is correlated with feature #1
#' ## feature #3-5: marginally related features and also internally correlated
#' ## feature #6-10: noise features without relationship with the y
#' 
#' set.seed(12345)
#' n <- 100
#' p <- 10
#' 
#' means <- rbind(
#'   c(0, 1, 0),
#'   c(0, 0, 0),
#'   c(0, 0, 1),
#'   c(0, 0, 1),
#'   c(0, 0, 1)
#' ) * 2
#' 
#' means <- rbind(means, matrix(0, p - 5, 3))
#' 
#' A <- diag(1, p)
#' A[1:5, 1:3] <- rbind(
#'   c(1, 0, 0),
#'   c(2, 1, 0),
#'   c(0, 0, 1),
#'   c(0, 0, 1),
#'   c(0, 0, 1)
#' )
#' 
#' dat <- gendata_FAM(n, means, A, sd_g = 0.5, stdx = TRUE)
#' ggplot2::qplot(dat$y, bins = 6)
#' corrplot::corrplot(cor(dat$X))
#' 
#' @seealso \code{\link{gendata_MLR}}
#'      
gendata_FAM <- function(n, muj, A, sd_g = 0, stdx = FALSE)
{
  p <- nrow(muj)
  C <- ncol(muj)
  k <- ncol(A)

  y <- rep(1:C, len = n)
  X <- A %*% matrix(rnorm(n * k), k, n) + muj[, y] + rnorm(n * p) * sd_g 
  SGM <- A %*% t(A) + diag(sd_g^2, p)
  
  if (stdx == TRUE)
  {
     mux <- rowMeans(muj)
     sdx <- sqrt(diag(SGM) + apply (muj, 1, var) * (C - 1) / C)
     muj <- (muj - mux) / sdx
     SGM <- sweep(SGM, 1, sdx, "/")
     SGM <- sweep(SGM, 2, sdx, "/")
     X <- (X - mux) / sdx
  }
  
  X <- t(X)
  colnames(X) <- paste0("V", 1:p)

  list("X" = X, "y" = y, "muj" = muj, "SGM" = SGM)
}

