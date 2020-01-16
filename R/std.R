#################################################### COPYLEFT NOTICE ####################################################
## This is a derivative work based on https://github.com/pbreheny/ncvreg/blob/master/R/std.R under GPLv3 licence.########
## Author of original version: Patrick Breheny, modified by Steven Liu. For more information see details section below.##
#########################################################################################################################

#' Standardizes a Design Matrix
#' 
#' This function accepts a design matrix and returns a standardized version of that matrix, 
#' the statistics of each column such as \code{median} and \code{sd} are also provided. 
#' 
#' @param X Design matrix, of dimension \code{nobs} by \code{nvars}; each row is an observation vector; 
#' can also be an object that can be coerced to a matrix, e.g. a data.frame.
#' 
#' @param tol The tolerance value; a column of \code{X} is considered as singular if the \code{sd}
#' of its entries (observations) is less than \code{tol}. Singular columns will be dropped by the end.
#' 
#' @return The standardized design matrix with the following attributes:
#' \describe{
#'   \item{nonsingular}{Indices of non-singular columns.}
#'   \item{center}{Median of each non-singular column which is used for standardization.}
#'   \item{scale}{Standard deviation of each non-singular column which is used for standardization.}
#' }
#' 
#' @details For each column of \code{X}, the standardization is done by first subtracting its median, 
#' then dividing by its sample standard deviation, while the original version in \code{ncvreg} uses
#' mean and population standard deviation. Its speed is slower than \code{ncvreg} because of the
#' complexity of median finding, but still substantially faster than \code{scale()} provided by R base.
#' 
#' @author Patrick Breheny (original) \cr Steven Liu (modification)
#' 
#' @seealso \url{http://pbreheny.github.io/ncvreg/reference/std.html} 
#' 
#' @export
#' 
#' @examples 
#' set.seed(123)
#' mat <- matrix(rnorm(n = 80 * 90, mean = 100, sd = 50), 80, 90)
#' mat %>% as.numeric() %>% ggplot2::qplot(bins = 30, xlab = '')
#' mat %>% std() %>% as.numeric() %>% ggplot2::qplot(bins = 30, xlab = '')
#'  
std <- function(X, tol = 1e-6)
{
  if (!is.matrix(X)) 
  {
    tmp <- try(X <- model.matrix(~0+., data = X), silent=TRUE)
    if (class(tmp)[1] == "try-error") 
      stop("X must be a matrix or able to be coerced to a matrix")
  }
  # checking NAs
  if (anyNA(X)) stop("data contains NA(s)")
  # calling std_helper() at src/utils.cpp
  std.info <- std_helper(X)
  # copying dimnames
  dimnames(std.info$X) <- dimnames(X)
  # checking singular columns
  vars.ns <- which(std.info$sd > tol)
  
  val <- std.info$X[, vars.ns, drop = FALSE]
  if (length(vars.ns) < ncol(X)) 
  {
    message(
      sprintf(
        "dropped %d singular column(s)",
        ncol(X) - length(vars.ns)
      )
    )
  }
  
  attr(val, "nonsingular") <- vars.ns
  attr(val, "center") <- std.info$median[, vars.ns]
  attr(val, "scale") <- std.info$sd[, vars.ns]
  return(val)
}
