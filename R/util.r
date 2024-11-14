#' Generate Prior Configuration
#' 
#' Configure prior hyper-parameters for HTLR model fitting
#' 
#' The output is a configuration list which is to be passed to \code{prior} argument of \code{htlr}.     
#' For naive users, you only need to specify the prior type and degree freedom, then the other hyper-parameters
#' will be chosen automatically. For advanced users, you can supply each prior hyper-parameters by yourself.
#' For suggestion of picking hyper-parameters, see \code{references}.  
#'
#' @param ptype The prior to be applied to the model. Either "t" (student-t, default), "ghs" (horseshoe), 
#' or "neg" (normal-exponential-gamma).
#' @param df The degree freedom (aka alpha) of t/ghs/neg prior for coefficients.
#' @param logw The log scale of priors for coefficients.
#' @param eta The \code{sd} of the normal prior for logw. When it is set to 0, logw is fixed. 
#' Otherwise, logw is assigned with a normal prior and it will be updated during sampling.
#' @param sigmab0 The \code{sd} of the normal prior for the intercept.
#' 
#' @return A configuration list containing \code{ptype}, \code{alpha}, \code{logw}, \code{eta}, and \code{sigmab0}.    
#' 
#' @references
#' Longhai Li and Weixin Yao. (2018). Fully Bayesian Logistic Regression 
#' with Hyper-Lasso Priors for High-dimensional Feature Selection.
#' \emph{Journal of Statistical Computation and Simulation} 2018, 88:14, 2827-2851.
#' 
#' @export
#' 
htlr_prior <- function(ptype = c("t", "ghs", "neg"), 
                       df = 1,
                       logw = -(1 / df) * 10,
                       eta = ifelse(df > 1, 3, 0),
                       sigmab0 = 2000)
{
  ptype <- match.arg(ptype)
  if (ptype != "t" & eta != 0)
    warning("random logw currently only supports t prior")
  list(
    "ptype" = ptype,
    "alpha" = df,
    "logw" = logw,
    "eta" = eta,
    "sigmab0" = sigmab0
  )
}

#' Split Data into Train and Test Partitions
#' 
#' This function splits the input data and response variables into training and testing parts.
#' 
#' @param X Input matrix, of dimension \code{nobs} by \code{nvars}; each row is an observation vector.
#' 
#' @param y Vector of response variables.
#' 
#' @param p.train Percentage of training set.
#' 
#' @param n.train Number of cases for training; will override \code{p.train} if specified.
#' 
#' @return List of training data \code{x.tr}, \code{y.tr} and testing data \code{x.te}, \code{y.te}.
#' 
#' @export
#' 
#' @examples 
#' dat <- gendata_MLR(n = 100, p = 10)
#' dat <- split_data(dat$X, dat$y, p.train = 0.7)
#' dim(dat$x.tr)
#' dim(dat$x.te)
#'    
split_data <- function(X,
                       y,
                       p.train = 0.7,
                       n.train = round(nrow(X) * p.train))
{
  stopifnot(nrow(X) == length(y), p.train > 0, p.train < 1)
  
  tr.row <- sample(1L:nrow(X), n.train, replace = FALSE)
  
  x.tr <- X[tr.row, , drop = FALSE]
  y.tr <- y[tr.row]
  x.te <- X[-tr.row, , drop = FALSE]
  y.te <- y[-tr.row]
  
  list(
    "x.tr" = x.tr, 
    "y.tr" = y.tr,
    "x.te" = x.te, 
    "y.te" = y.te
  )
}

# compute sd of betas
# @export
comp_sdb <- function (deltas, removeint = TRUE, normalize = FALSE)
{
    C <- ncol (deltas) + 1
    if (removeint)
    {
        deltas <- deltas[-1,,drop = F]
    }
    
    vardeltas <- comp_vardeltas (deltas)
    sdb <- sqrt (vardeltas/C)
    
    if (normalize) sdb <- sdb / max(sdb)
    
    sdb
}

#' @export
nobs.htlr.fit <- function(object, ...)
{
  object$n
}

#' @export
print.htlr.fit <- function(x, ...)
{
  info.data <- sprintf("Data:\n
  response:\t%d-class
  observations:\t%d
  predictors:\t%d (w/ intercept)
  standardised:\t%s",
  x$K + 1, x$n, x$p + 1, as.character(x$feature$stdx)) 
  
  info.model <- sprintf("Model:\n
  prior dist:\t%s (df = %d, log(w) = %.1f)
  init state:\t%s 
  burn-in:\t%d
  sample:\t%d (posterior sample size)",
  x$prior$ptype, x$prior$alpha, x$prior$logw, 
  x$mc.param$init, x$mc.param$iter.warm, x$mc.param$iter.rmc)
  
  info.est <- sprintf("Estimates:\n
  model size:\t%d (w/ intercept)
  coefficients: see help('summary.htlr.fit')",
  length(nzero_idx(x)) + 1)
  
  cat("Fitted HTLR model", "\n\n", info.data, "\n\n", info.model, "\n\n", info.est)
}

# try to install suggested packages when needed
# @author: Michael W. Kearney
# @source: https://github.com/ropensci/rtweet/blob/master/R/utils.R
try_require <- function(pkg, f = NULL) {
  if (is.null(f))
    f <- "this action"
  else
    f <- paste0("`", f, "`")
  
  if (requireNamespace(pkg, quietly = TRUE)) 
  {
    library(pkg, character.only = TRUE)
    return(invisible())
  }
  
  stop(paste0("Package `", pkg, "` required for ", f , ".\n",
              "Please install and try again."), call. = FALSE)
}

# @param rstudio launch in rstudio viewer instead of web browser? 
# @param ... passed to shiny::runApp
# launch_shiny <- function(launch.browser = TRUE, ...) {
#   try_require("shiny", f = "launch_shiny()")
#   shiny::runApp(system.file("app", package = "HTLR"), 
#                 launch.browser = launch.browser, ...)
# }

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#Plots feature importance scores
# 
#This function plots feature importance scores or coefficients using histogram line. 
# 
#param fscores Scores measuring feature importance, such as \code{wsdbs}, \code{msdbs}, or coefficients values.
# 
#export
# 
#seealso htlr_fss
# plot_fscore <- function (fscores, fsel=1:length (fscores), show_ix = 0.1, do.plot = TRUE, ...)
# {
# 	if (show_ix > 1) stop ("show_ix must be less than 1")
#     afscores <- abs (fscores)
#     mfscore <- max (afscores)
#     
#     plotargs <- list (...)
#     
#     p <- length (fscores) 
# 
#     if (is.null (plotargs$log))  plotargs$log <- ""
#     if (is.null (plotargs$type)) plotargs$type <- "h"
#     
#     if (is.null (plotargs$ylab) ) plotargs$ylab = "Feature Score"
#     if (is.null(plotargs$xlab))plotargs$xlab <-"Feature Index" 
#     if (is.null (plotargs$cex.axis))  plotargs$cex.axis <- 0.8
#     
#     # plot fscores   
#     if (do.plot) {
# 	    do.call (plot, c(list (x = fscores), plotargs) )    
#     		# show shresholds 0.1 and 0.01
#     		abline (h = mfscore * c(-0.01,0.01), lty = 2, col = "grey")
#     		abline (h = mfscore * c(-0.1,0.1), lty = 1, col = "grey")    
#     }
#     
#     # showtops       
#     itops <- which (afscores >= show_ix * mfscore)    
# 	if (do.plot)
# 	 text (itops, fscores [itops], fsel[itops], col = "red", srt = 0, adj = - 0.2, cex = 0.7)
#     
#     a <- fsel[itops]
#     
# }
