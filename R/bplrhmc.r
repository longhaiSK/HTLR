#' Fit a HTLR models and make predictions 
#'
#' This function trains linear logistic regression models with HMC in restricted Gibbs sampling. 
#' It also makes predictions for test cases if \code{X_ts} are provided.
#'
#' @param y_tr Vector of class labels in training or test data set. 
#' Must be coded as non-negative integers, e.g., 1,2,\ldots,C for C classes.
#' @param X_tr Design matrix of traning data; 
#' rows should be for the cases, and columns for different features.
#' @param fsel Subsets of features selected before fitting, such as by univariate screening.
#' @param stdzx Logical; if \code{TRUE}, the original feature values are standardized to have \code{mean} = 0 
#' and \code{sd} = 1.
#' 
#' @param iter_h A positive integer specifying the number of warmup (aka burnin).
#' @param iters_rmc A positive integer specifying the number of iterations after warmup.
#' @param thin A positive integer specifying the period for saving samples.
#' 
#' @param leap_L The length of leapfrog trajectory in sampling phase.
#' @param leap_L_h The length of leapfrog trajectory in burnin phase.
#' @param leap_step The stepsize adjustment multiplied to the second-order partial derivatives of log posterior.
#' 
#' @param initial_state The initial state of Markov Chain; can be NULL, 
#' or a gived parameter vector, or a previous markov chain results.
#' 
#' @param ptype The prior to be applied to the model. Either "t" (default), "ghs" (horseshoe), 
#' or "neg" (normal-exponential-gamma).
#' 
#' @param sigmab0 The \code{sd} of the normal prior for the intercept.
#' @param alpha The degree freedom of t/ghs/neg prior for coefficients.
#' @param s The log scale of priors (logw) for coefficients.
#' @param eta The \code{sd} of the normal prior for logw. When it is set to 0, logw is fixed. 
#' Otherwise, logw is assigned with a normal prior and it will be updated during sampling. 
#' 
#' @param hmc_sgmcut The coefficients smaller than this criteria will be fixed in each HMC updating step.
#' 
#' @param alpha.rda A user supplied alpha value for \code{bcbcsf_deltas}. The default is 0.2.
#' 
#' @param silence Setting it to \code{FALSE} for tracking MCMC sampling iterations. 
#' @param pre.legacy Logical; if \code{TRUE}, the output produced in \code{HTLR} versions up to 
#' legacy-3.1-1 is reproduced. 
#' 
#' @return A list of fitting results.  
#' 
#' @references
#' Longhai Li and Weixin Yao. (2018). Fully Bayesian Logistic Regression 
#' with Hyper-Lasso Priors for High-dimensional Feature Selection.
#' \emph{Journal of Statistical Computation and Simulation} 2018, 88:14, 2827-2851.
#' 
#' @useDynLib HTLR
#' 
#' @import Rcpp
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' data (gen_grpcor_data)
#' ## creating training and test data sets
#' tr <- 200
#' X_tr <- data$X[1:tr,]
#' y_tr <- data$y[1:tr]
#' X_ts <- data$X[-(1:tr), ]
#' y_ts <- data$y[-(1:tr)]
#' ## fit htlr models
#' fithtlr <- htlr_fitpred (
#'   y_tr = y_tr, X_tr = X_tr, stdzx = TRUE, #fsel = c (1,51,101),
#'   pty = "t", alpha = 1, s = -15, 
#'   iters_h = 1000, iters_rmc = 1000, thin = 50,
#'   leap_L_h = 5, leap_L = 50, leap_step = 0.5, hmc_sgmcut = 0.05, 
#'   initial_state = "bcbcsfrda", silence = F)
#' }
#'
htlr_fit <- function (
    y_tr, X_tr, X_ts = NULL, fsel = 1:ncol(X_tr), stdzx = TRUE, ## data
    sigmab0 = 2000, ptype = "t", alpha = 1, s = -10, eta = 0,  ## prior
    iters_h = 1000, iters_rmc = 1000, thin = 1,  ## mc iterations
    leap_L = 50, leap_L_h = 5, leap_step = 0.3,  hmc_sgmcut = 0.05, ## hmc
    initial_state = "lasso", alpha.rda = 0.2, silence = TRUE, pre.legacy = TRUE, ## initial state
    predburn = NULL, predthin = 1) ## prediction
{
  #.Deprecated("htlr")
  
  ## checking prior types
  if (!(ptype %in% c("t", "ghs", "neg"))) 
  {
    stop ("\"ptype\" NOT in (\"t\", \"ghs\", \"neg\")")
  }
  ## checking arguments
  if (iters_rmc <= 0||iters_h < 0||leap_L <=0 ||leap_L_h <= 0 || thin <= 0)
  {
    stop ("MC iterations and Leapfrog lengths must be nonnegative.")
  }
  
  if (length (y_tr) != nrow (X_tr) ) stop ("'y' and 'X' mismatch")

  ###################### Data preprocessing ######################
  if (min(y_tr) == 0)
    y_tr <- y_tr + 1
  
  ybase <- as.integer(y_tr - 1)
  ymat <- model.matrix( ~ factor(y_tr) - 1)[, -1]
  C <- length(unique(ybase))
  K <- C - 1
  ## feature selection
  X <- X_tr[, fsel, drop = FALSE]
  p <- length(fsel)
  n <- nrow(X)
  ## standardize selected features
  nuj <- rep(0, length(fsel))
  sdj <- rep(1, length(fsel))
  if (stdzx == TRUE & !is.numeric(initial_state))
  {
    nuj <- apply(X, 2, median)
    sdj <- apply(X, 2, sd)
    X <- sweep(X, 2, nuj, "-")
    X <- sweep(X, 2, sdj, "/")
  }
  ## add intercept
  X_addint <- cbind(1, X)
  ## stepsize for HMC from data
  DDNloglike <- 1 / 4 * colSums(X_addint ^ 2)

  #################### Markov chain state initialization ####################
  ## starting from a given deltas
  if (is.list(initial_state)) # use the last iteration of markov chain
  {
    no_mcspl <- length(initial_state$mclogw)
    deltas <- matrix(initial_state$mcdeltas[, , no_mcspl], nrow = p + 1)
    sigmasbt <- initial_state$mcsigmasbt[, no_mcspl]
    logw <- initial_state$mclogw[no_mcspl]
  }
  else if (is.numeric(initial_state))
  {
    deltas <- matrix(initial_state, nrow = p + 1)
    logw <- s
  }
  else if (initial_state == "lasso")
  {
    lambda <- ifelse(pre.legacy, NA, .01)
    deltas <- lasso_deltas(X = X, y = y_tr, lambda)
    logw <- s
  }
  else if (initial_state == "bcbcsfrda")
  {
    deltas <- bcbcsf_deltas(X, y_tr, alpha = alpha.rda)
    logw <- s
  }
  else if (initial_state == "random")
  {
    deltas <- matrix(rnorm((p + 1) * K) * 2, p + 1, K)
    logw <- s
  }
  if (nrow (deltas) != p + 1 || ncol (deltas) != K)
    stop("Initial `deltas' Mismatch Data")
  
  if (!exists("sigmasbt"))
  {
    vardeltas <- comp_vardeltas(deltas)[-1]
    sigmasbt <- c(sigmab0, spl_sgm_ig (alpha, K, exp(logw), vardeltas))
  }
    
  #################### Do Gibbs sampling ####################

  fit <- HtlrFit(
      ## data
      p = p, K = K, n = n,
      X = as.matrix(X_addint), 
      ymat = as.matrix(ymat), 
      ybase = as.vector(ybase),
      ## prior
      ptype = ptype, alpha = alpha, s = s, eta = eta, sigmab0 = sigmab0,
      ## sampling
      iters_rmc = iters_rmc, iters_h = iters_h, thin = thin, 
      leap_L = leap_L, leap_L_h = leap_L_h, leap_step = leap_step, 
      hmc_sgmcut = hmc_sgmcut, DDNloglike = as.vector(DDNloglike),
      ## fit result
      deltas = deltas, logw = logw, sigmasbt = sigmasbt,
      ## other control
      silence = as.integer(silence), legacy = pre.legacy)
  
  # adding data preprocessing information
  fit <- c(fit, list(fsel = fsel, nuj = nuj, sdj = sdj, y = y_tr))
        
  ################## Prediction for test cases #########################
  if (!is.null (X_ts))
  {
    fit$probs_pred <- htlr_predict(
      X_ts = X_ts,
      fithtlr = fit,
      burn = predburn,
      thin = predthin
    )
  }
  else
    fit$probs_pred <- NULL
  
  ############## Htlr fitting and prediction results #################
  return(fit)
}

#' Make prediction on new data
#' 
#' This function uses MCMC samples returned by \code{htlr_fit} to predict the class labels of test cases. 
#' 
#' @param X_ts Matrix of values at which predictions are to be made.
#' @param fithtlr Fitted HTLR model object.
#' @param deltas The values of deltas (for example true deltas) used to prediction. 
#' 
#' @return A matrix of predictive probabilities, with rows for cases, cols for classes.
#' 
#' @export
#' 
#' @seealso htlr_fit  
htlr_predict <- function(X_ts, fithtlr = NULL, deltas = NULL, burn = NULL, thin = NULL, usedmc = NULL)
{
  ## chaning X_ts as needed
  if (is.vector(X_ts))
    X_ts <- matrix(X_ts, 1)
  no_ts <- nrow(X_ts)
  
  if (is.null(deltas) & !is.null(fithtlr))
  {
    mcdims <- dim(fithtlr$mcdeltas)
    p <- mcdims[1] - 1
    K <- mcdims[2]
    no_mcspl <- mcdims[3]
    
    ## index of mc iters used for inference
    if (is.null(usedmc))
    {
      if (is.null(burn))
        burn <- floor(no_mcspl * 0.1)
      if (is.null(thin))
        thin <- 1
      usedmc <- seq(burn + 1, no_mcspl, by = thin)
    }
    
    no_used <- length(usedmc)
    
    ## read deltas for prediction
    longdeltas <- matrix(fithtlr$mcdeltas[, , usedmc], nrow = p + 1)
    
    ## selecting features and standardizing
    fsel <- fithtlr$fsel
    X_ts <- X_ts[, fsel, drop = FALSE]
    nuj <- fithtlr$nuj
    sdj <- fithtlr$sdj
    X_ts <- sweep(X_ts, 2, nuj, "-")
    X_ts <- sweep(X_ts, 2, sdj, "/")
  }
  else
  {
    if (is.vector(deltas) | is.matrix (deltas)) 
    {
      deltas <- matrix(deltas, nrow = ncol (X_ts) + 1)
      p <- nrow(deltas) - 1
      K <- 1
      longdeltas <- deltas
      no_used <- 1
    }
  }
  
  ## add intercept to all cases
  X_addint_ts <- cbind(1, X_ts)
  
  longlv <- X_addint_ts %*% longdeltas
  arraylv <- array(longlv, dim = c(no_ts, K, no_used))
  logsumlv <- apply(arraylv, 3, comp_lsl)
  array_normlv <- sweep(arraylv, c(1, 3), logsumlv)
  array_predprobs <- exp(array_normlv)
  probs_pred <- apply(array_predprobs, c(1, 2), mean)
  
  predprobs_c1 <- pmax(0, 1 - apply(probs_pred, 1, sum))
  probs_pred <- cbind(predprobs_c1, probs_pred)
  
  return(probs_pred)
}


######################## some functions not used currently ###################

if (F)
{

htlr_ci <- function (fithtlr, usedmc = NULL)
{
    mcdims <- dim (fithtlr$mcdeltas)
    p <- mcdims [1] - 1
    K <- mcdims [2]
    no_mcspl <- mcdims[3]

    ## index of mc iters used for inference

    mcdeltas <- fithtlr$mcdeltas[,,usedmc, drop = FALSE]
    
    cideltas <- array (0, dim = c(p+1, K, 3))
    for (j in 1:(p+1))
    {
        for (k in 1:K) {
          cideltas [j,k,] <- 
            quantile (mcdeltas[j,k,], probs = c(1-cp, 1, 1 + cp)/2)
        }
    }
    
    cideltas
}

## this function plots confidence intervals
htlr_plotci <- function (fithtlr, usedmc = NULL, 
                         cp = 0.95, truedeltas = NULL,   ...)
{
    
    cideltas <- htlr_coefs (fithtlr, usedmc = usedmc, showci = TRUE, cp = cp)
    K <- dim (cideltas)[2]
    
    for (k in 1:K)
    {
        plotmci (cideltas[,k,], truedeltas = truedeltas[,k], 
                 main = sprintf ("%d%% MC C.I. of Coefs (Class %d)", 
                                 cp * 100, k+1),
                ...)
        
    }
    
    return (cideltas)
}


htlr_outpred <- function (x,y,...)
{
  X_ts <- cbind (x, rep (y, each = length (x)))
  probs_pred <- htlr_predict (X_ts = X_ts, ...)$probs_pred[,2] 
  matrix (probs_pred, nrow = length (x) )
}


norm_coef <- function (deltas)
{
  slope <- sqrt (sum(deltas^2))
  deltas/slope
}

pie_coef <- function (deltas)
{
  slope <- sum(abs(deltas))
  deltas/slope
}

norm_mcdeltas <- function (mcdeltas)
{
  sqnorm <- function (a) sqrt(sum (a^2))
  dim_mcd <- dim (mcdeltas)
    
  slopes <- apply (mcdeltas[-1,,,drop=FALSE], MARGIN = c(2,3), sqnorm)
    
  mcthetas <- sweep (x = mcdeltas, MARGIN = c(2,3), STATS = slopes, FUN = "/")
  
  list (mcthetas = mcthetas, slopes = as.vector(slopes))
}

pie_mcdeltas <- function (mcdeltas)
{
  sumabs <- function (a) sum (abs(a))
  dim_mcd <- dim (mcdeltas)
    
  slopes <- apply (mcdeltas[-1,,,drop=FALSE], MARGIN = c(2,3), sumabs)
    
  mcthetas <- sweep (x = mcdeltas, MARGIN = c(2,3), STATS = slopes, FUN = "/")
  
  list (mcthetas = mcthetas, slopes = as.vector(slopes))
}

plotmci <- function (CI, truedeltas = NULL, ...)
{
    p <- nrow (CI) - 1

    plotargs <- list (...)
    
    if (is.null (plotargs$ylim)) plotargs$ylim <- range (CI)
    if (is.null (plotargs$pch))  plotargs$pch <- 4 
    if (is.null (plotargs$xlab)) 
       plotargs$xlab <- "Feature Index in Training Data"
    if (is.null (plotargs$ylab)) plotargs$ylab <- "Coefficient Value"
    
    do.call (plot, c (list(x= 0:p, y=CI[,2]), plotargs))
    
    abline (h = 0)
    
    for (j in 0:p)
    {
        
        points (c(j,j), CI[j+1,-2], type = "l", lwd = 2)
    }
    
    if (!is.null (truedeltas))
    {
        points (0:p, truedeltas, col = "red", cex = 1.2, pch = 20)
    }

}




htlr_plotleapfrog <- function ()
{
        if (looklf & i_mc %% iters_imc == 0 & i_mc >=0 )
        {
           if (!file.exists ("leapfrogplots")) dir.create ("leapfrogplots")

           postscript (file = sprintf ("leapfrogplots/ch%d.ps", i_sup),
           title = "leapfrogplots-ch", paper = "special",
           width = 8, height = 4, horiz = FALSE)
           par (mar = c(5,4,3,1))
           plot (-olp$nenergy_trj + olp$nenergy_trj[1],
                xlab = "Index of Trajectory", type = "l",
                ylab = "Hamiltonian Value",
                main =
                sprintf (paste( "Hamiltonian Values with the Starting Value",
                "Subtracted\n(P(acceptance)=%.2f)", sep = ""),
                min(1, exp(olp$nenergy_trj[L+1]-olp$nenergy_trj[1]) )
                )
           )
           abline (h = c (-1,1))
           dev.off()

           postscript (file = sprintf ("leapfrogplots/dd%d.ps", i_sup+1),
           title = sprintf("leapfrogplots-dd%d", i_sup + 1), 
           paper = "special",
           width = 8, height = 4, horiz = FALSE)
           par (mar = c(5,4,3,1))
           plot (olp$ddeltas_trj, xlab = "Index of Trajectory",type = "l",
                 ylab = "square distance of Deltas",
                 main = "Square Distance of `Deltas'")
           dev.off ()

           postscript (file = sprintf ("leapfrogplots/ll%d.ps", i_sup),
           title = "leapfrogplots-ll", paper = "special",
           width = 8, height = 4, horiz = FALSE)
           par (mar = c(5,4,3,1))
           plot (olp$loglike_trj, xlab = "Index of Trajectory", type = "l",
                 ylab = "log likelihood",
                 main = "Log likelihood of Training Cases")
           dev.off()
        }
}



}


