
########################### utility functions ###############################


## a function for retrieve fithtlr objs saved in a RData file
# reload_fithtlr <- function (fithtlrfile)
# {
#     local
#     ({
#         fithtlr <- get (load (fithtlrfile))
#         return (fithtlr)
#     })
# }


log_sum_exp <- function(lx)
{  
  mlx <- max(lx)
  log(sum(exp(lx - mlx))) + mlx
}

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


comp_lsl <- function (lv)
{
    apply (cbind (0,lv), 1, log_sum_exp)
}

# spl_sgm_ig <- function (alpha, K, w, vardeltas)
# {
#   1 / rgamma (length (vardeltas), (alpha + K)/2) * (alpha * w + vardeltas) / 2
# }

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