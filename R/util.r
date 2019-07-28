#' Generate prior configuration
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
#' @seealso htlr
htlr_prior <- function(ptype = c("t", "ghs", "neg"), 
                         df = 1,
                         logw = -(1 / df) * 10,
                         eta = (df > 2) * 100,
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

## a function for retrieve fithtlr objs saved in a RData file
# reload_fithtlr <- function (fithtlrfile)
# {
#     local
#     ({
#         fithtlr <- get (load (fithtlrfile))
#         return (fithtlr)
#     })
# }

## compute V (delta)
comp_vardeltas <- function (deltas)
{
    K <- ncol (deltas)
    SUMdeltas <- rowSums (deltas)
    SUMsqdeltas <- rowSums (deltas^2)
    SUMsqdeltas  - SUMdeltas^2 / (K + 1)
}

## compute sd of betas
#' @export
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


comp_lsl <- function(lv)
{
  log_sum_exp(cbind(0, lv))
}

log_normcons <- function(lv)
{
  sum(comp_lsl(lv))
}

get_ix <- function (sub, whole, digits= 0)
{
    p <- length (whole)
    wix <- 1:p
    names (wix) <- as.character (round(whole, digits))
    wix [as.character (round(sub, digits))]
}



#' Plots feature importance scores
#' 
#' This function plots feature importance scores or coefficients using histogram line. 
#' 
#' @param fscores Scores measuring feature importance, such as \code{wsdbs}, \code{msdbs}, or coefficients values.
#' 
#' @export
#' 
#' @seealso htlr_fss
plot_fscore <- function (fscores, fsel=1:length (fscores), show_ix = 0.1, do.plot = TRUE, ...)
{
	if (show_ix > 1) stop ("show_ix must be less than 1")
    afscores <- abs (fscores)
    mfscore <- max (afscores)
    
    plotargs <- list (...)
    
    p <- length (fscores) 

    if (is.null (plotargs$log))  plotargs$log <- ""
    if (is.null (plotargs$type)) plotargs$type <- "h"
    
    if (is.null (plotargs$ylab) ) plotargs$ylab = "Feature Score"
    if (is.null(plotargs$xlab))plotargs$xlab <-"Feature Index" 
    if (is.null (plotargs$cex.axis))  plotargs$cex.axis <- 0.8
    
    # plot fscores   
    if (do.plot) {
	    do.call (plot, c(list (x = fscores), plotargs) )    
    		# show shresholds 0.1 and 0.01
    		abline (h = mfscore * c(-0.01,0.01), lty = 2, col = "grey")
    		abline (h = mfscore * c(-0.1,0.1), lty = 1, col = "grey")    
    }
    
    # showtops       
    itops <- which (afscores >= show_ix * mfscore)    
	if (do.plot)
	 text (itops, fscores [itops], fsel[itops], col = "red", srt = 0, adj = - 0.2, cex = 0.7)
    
    a <- fsel[itops]
    
}

## try to install suggested packages when needed
try_require <- function(pkg, f = NULL) {
  if (is.null(f)) {
    f <- "this action"
  } else {
    f <- paste0("`", f, "`")
  }
  
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    return(invisible())
  }
  
  stop(paste0("Package `", pkg, "` required for ", f , ".\n",
              "Please install and try again."), call. = FALSE)
}