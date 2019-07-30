#' Bias-corrected Bayesian classification initial state
#' 
#' Generate initial Markov chain state with Bias-corrected Bayesian classification.
#' 
#' Caveat: This method can be used only for continuous predictors such as gene expression profiles, 
#' and it does not make sense for categorical predictors such as SNP profiles. 
#' 
#' @param X Design matrix of traning data; 
#' rows should be for the cases, and columns for different features.
#' 
#' @param y Vector of class labels in training or test data set. 
#' Must be coded as non-negative integers, e.g., 1,2,\ldots,C for C classes.
#' 
#' @param alpha The regularization proportion (between 0 and 1) for mixing the 
#' diagonal covariance estimates and the sample covariance estimated with the 
#' training samples. The default is 0, the covariance matrix is assumed to be diagonal, 
#' which is the most robust.
#' 
#' @return A matrix - the initial state of Markov Chain for HTLR model fitting.
#' 
#' @references Longhai Li (2012). Bias-corrected hierarchical Bayesian classification 
#' with a selected subset of high-dimensional features. 
#'\emph{Journal of the American Statistical Association}, 107(497), 120-134.
#' 
#' @importFrom BCBCSF bcbcsf_fitpred bcbcsf_sumfit  
#' 
#' @export
#' 
#' @keywords internal
#' 
#' @seealso \code{\link{lasso_deltas}} 
#' 
bcbcsf_deltas <- function(X, y, alpha = 0)
{
  n <- nrow(X)
  p <- ncol(X)
  G <- length(unique(y))
  
  if (alpha >= 1 & n - G <= p)
    stop ("'alpha' must be less than 1")
  
  muj <- matrix(0, p, G)
  ## estimate centroids
  fit_bcbcsf <- bcbcsf_fitpred(
    X_tr = X,
    y_tr = y,
    alpha1_mu = 1,
    standardize = FALSE,
    no_rmc = 500,
    no_imc = 1,
    no_mhwmux = 10,
    bcor = 0,
    fit_bcbcsf_filepre = NULL,
    monitor = F
  )$fit_bcbcsf
  
  muj[fit_bcbcsf$fsel, ] <- bcbcsf_sumfit(fit_bcbcsf = fit_bcbcsf)$muj
  
  ## transform X
  tX_tr <- t(X)
  
  ## proportions of groups
  gprior <- table(y) / n
  
  for (g in 1:G)
  {
    ## subtract muj[,g] from data in gth group
    tX_tr[, y == g] <- sweep(tX_tr[, y == g, drop = FALSE], 1, muj[, g])
  }
  X <- t(tX_tr)
  
  ## estimate inverse of covariance matrix
  lambda <- alpha / (1 - alpha) / n
  inv_SGM <- 1 / (1 - alpha) * 
    (diag(1, p) - lambda * tX_tr %*% solve(diag(1, n) + lambda * X %*% tX_tr) %*% X)
  
  deltas_bc_inv(muj = muj, inv_SGM = inv_SGM, gprior = gprior)
}

deltas_bc_inv <- function (muj, inv_SGM, gprior = rep(1, ncol (muj)))
{
  G <- ncol (muj)
  p <- nrow (muj)
  
  ## find coefficients of discriminative functions
  betas <- matrix (0, p + 1, G)
  for (g in 1:G)
  {
    ## coefficients for variables
    betas[-1, g] <- inv_SGM %*% muj[, g]
    ## intercept
    betas [1, g] <-
      log (gprior[g]) - 0.5 * sum (muj[, g] * betas[-1, g])
  }
  
  ## find 'deltas'
  deltas <- betas[, -1, drop = FALSE] - betas[, 1]
  
  deltas
}


deltas_bc <- function (muj, SGM, subset = NULL, gprior = NULL)
{
  G <- ncol (muj)
  p <- nrow (muj)
  
  if (is.null(gprior))
    gprior <- rep (1 / G, G)
  gprior <- gprior / sum(gprior)
  
  if (is.null (subset))
    subset <- 1:p
  
  muj <- muj [subset, , drop = FALSE]
  SGM <- SGM [subset, subset, drop = FALSE]
  
  ## find coefficients of discriminative functions
  betas <- matrix (0, length (subset) + 1, G)
  for (g in 1:G)
  {
    ## coefficients for variables
    betas[-1, g] <- solve(SGM) %*% muj[, g]
    ## intercept
    betas [1, g] <-
      log (gprior[g]) - 0.5 * sum (muj[, g] * betas[-1, g])
  }
  
  ## find 'deltas'
  betas[,-1, drop = FALSE] - betas[, 1]
  
}
#require (abind)

#############################################################################
########################## Utility Functions ################################
#############################################################################

### draw random numbers from inverse gamma distribution
#rinvgam <- function (n, alpha, lambda)
#{
#    1 / rgamma (n, alpha, 1) * lambda
#}

### define a binary operator doing row mulplication
#"%r%" <- function (A, a) sweep (A, 2, a, "*")

### this is a generic function for generating Markov chain samples
### from a given density with Metropolis method
#met_gauss <- function
#   ( iters = 100, log_f, ini_value, stepsize = 0.5, diag_mh = FALSE, ...)
#{
#    state <- ini_value
#    no_var <- length (state)
#    mchain <- matrix (state, no_var , iters)
#    nos_rej <- 0
#    logf <- log_f (state,...)

#    if (!is.finite (logf)) stop ("Initial value has 0 probability")

#    for (i in 1:iters)
#    {
#        new_state <- rnorm (no_var, state, stepsize)
#        new_logf <- log_f (new_state,...)

#        if (log (runif(1)) < new_logf - logf)
#        {
#            state <- new_state
#            logf <- new_logf
#        }
#        else nos_rej <- nos_rej + 1

#        ## save state in chain
#        mchain[,i] <- state
#    }

#    if (diag_mh)
#    {
#      cat ("Markov chain is saved in 'mchain' with columns for iterations\n")
#      cat ("Rejection rate = ", nos_rej / iters, "\n")
#      browser () ## pause to allow user to look at Markov chain
#    }
#    state
#}


#############################################################################
########################## parameter estimation #############################
#############################################################################

### This function estimates the parameters and hyerparameters based on
### the mle of means and variances of each feature.
#trpr_mle <- function (X_tr, y_tr, X_ts = NULL, nos_fsel = ncol (X_tr), 
#            rankf = rank_plain)
#{
#    ## read information about data
#    n <- nrow (X_tr) ## numbers of obs
#    p <- ncol (X_tr)
#    ## find number of observations in each group
#    nos_g <- as.vector (tapply (rep(1,n),INDEX = y_tr, sum))
#    G <- length (nos_g)
#    if (any(nos_g < 2)) stop ("Less than 2 cases in some group")


#    ## create result storage
#    nnfsel <- length (nos_fsel)

#    list_fit_mle <- rep (list (""), nnfsel)
#    if (!is.null (X_ts))
#        array_probs_pred <- array (0, dim = c (nrow (X_ts), G, nnfsel) )
#    else array_probs_pred <- NULL

#    ## feature selection
#    rankedf <- rankf (X_tr, y_tr)

#    for (i in 1:nnfsel)
#    {
#        k <- nos_fsel [i]
#        fsel <- rankedf [1:k]

#        X_tr_sel <- X_tr [, fsel, drop = FALSE]

#        ## sufficient statistic and pooled variances
#        gsum_X <- rowsum (X_tr_sel,y_tr)
#        sum_X2 <- colSums (X_tr_sel^2)
#        sum_gsum2 <- colSums (gsum_X^2 / nos_g )
#        pvars <- (sum_X2 - sum_gsum2) / (n-G) ## pooled variances

#        muj <- t(gsum_X / nos_g) ## group means
#        wxj <- pvars ## pooled variances
#        wx <- 1/mean (1/wxj)
#        nuj <- rowMeans (muj)
#        wnu <- mean (nuj^2)
#        cmuj <- muj - nuj
#        wmuj <- rowMeans (cmuj^2)
#        wmu <- 1/mean (1/wmuj)

#        fit_mle <- list (
#            muj = muj, wxj = wxj, wmuj = wmuj, nuj = nuj, cmuj = cmuj,
#            wx = wx, wmu = wmu, wnu = wnu, freqy = nos_g / sum (nos_g),
#            fsel = fsel)

#        list_fit_mle [[i]] <- fit_mle
#        ## note: in "muj", the row is for features, the column is for groups
#        if (!is.null (X_ts))
#        {
#            array_probs_pred [,,i] <- mlepred (X_ts = X_ts, fit_mle = fit_mle)
#        }
#    }

#    list (
#      nos_fsel = nos_fsel, array_probs_pred = array_probs_pred,
#      list_fit_mle = list_fit_mle)
#}


#trpr_hb <- function (
#    ## arguments specifying info of data sets
#    X_tr, y_tr, X_ts = NULL, nos_fsel = ncol (X_tr), standardize = TRUE,
#    ## arguments for Markov chain sampling
#    no_rmc = 1000, no_imc = 10, no_mhwmux = 10,
#    fit_bcbcsf_prefile = ".fit_bcbcsf_",  ini_mcstate = NULL,
#    ## arguments specifying priors for parameters and hyerparameters
#    w0_mu = 0.05, alpha0_mu = 0.5, alpha1_mu = 2,
#    w0_x  = 0.05, alpha0_x  = 0.5, alpha1_x  = 10,
#    w0_nu = 0.05, alpha0_nu = 0.5, prior_psi = NULL,
#    ## arguments for metropolis sampling for wmu, wx
#    stepadj_mhwmux = 1, diag_mhwmux = FALSE,
#    ## arguments for computing adjustment factor
#    bcor = 1, cut_qf = exp (-10), cut_dpoi = exp (-10), nos_sim = 1000,
#    ## arguments for prediction
#    burn = 200, thin = 10
#)

#{
#    ## read information about data
#    n <- nrow (X_tr) ## numbers of obs
#    p <- ncol (X_tr)
#    ## find number of observations in each group
#    nos_g <- as.vector (tapply (rep(1,n), INDEX = y_tr, sum))
#    G <- length (nos_g)
#    if (any(nos_g < 2)) stop ("less than 2 cases in some group in your data")

#    ## set prior for proportion of cases
#    if (is.null (prior_psi))
#    {
#      prior_psi <- rep (1, G)
#    }
#    else
#    {
#      if (length (prior_psi) != G) stop ("length of prior_psi is wrong")
#    }

#    ## starndardize features
#    nuj_std <- rep (0,p)
#    wxj_std <- rep (1,p)
#    if ( standardize )
#    {
#        mlepar <- trpr_mle (X_tr = X_tr, y_tr = y_tr)$list_fit_mle[[1]]
#        nuj_std <- mlepar $ nuj
#        wxj_std <- mlepar $ wxj
#        X_tr <- sweep (X_tr, 2, nuj_std)
#        X_tr <- sweep (X_tr, 2, sqrt (wxj_std), "/")
#    }

#    ## feature selection
#    info_sel <- list (vars = 1:p, fstats = rep (1, p))
#    ## creating result storage
#    nnfsel <- length (nos_fsel)
#    fit_bcbcsf_files <- rep ("", nnfsel)
#    if (!is.null (X_ts))
#        array_probs_pred <- array (0, dim = c(nrow (X_ts), G, nnfsel) )
#    else array_probs_pred <- NULL
#    time_train <- time_pred <- rep (0, nnfsel)

#    for (i  in seq (1, nnfsel) )
#    {
#        k <- nos_fsel [i]

#        fsel <- info_sel $ vars [1:k]
#        cut_F <- info_sel $ fstats [k]
#        X_tr_sel <- X_tr [, fsel, drop = FALSE]
#        nos_omit <- p - k

#        ## initialize Markov chain variables
#        ## note: in "muj", the rows for diff features, columns for diff groups
#        if (is.null (ini_mcstate) )
#        {
#            ## start from mle
#            mlepar <- trpr_mle (X_tr = X_tr_sel, y_tr = y_tr,
#                      nos_fsel = k)$list_fit_mle[[1]]
#            muj <- mlepar $ muj
#            wxj <- mlepar $ wxj
#            wmuj <- mlepar $ wmuj
#            nuj <- mlepar $ nuj
#            wx <- mlepar $ wx
#            logwx <- log (wx)
#            wmu <- mlepar $ wmu
#            logwmu <- log (wmu)
#            wnu <- mlepar $ wnu
#        }
#        else
#        {
#            ## start from a saved iteration
#            iters_inimcstate <- length (ini_mcstate $ wx)
#            muj <- ini_mcstate $ MUJ[ , , iters_inimcstate]
#            wxj <- ini_mcstate $ WXJ[, iters_inimcstate]
#            wmuj <- ini_mcstate $ WMUJ [ , iters_inimcstate]
#            nuj <- ini_mcstate $ NUJ [ , iters_inimcstate]
#            wx <- ini_mcstate $ WX [iters_inimcstate]
#            logwx <- log (wx)
#            wmu <- ini_mcstate $ WMU [iters_inimcstate]
#            logwmu <- log (wmu)
#            wnu <- ini_mcstate $ WNU [iters_inimcstate]
#        }

#        ## Marlov chain storage
#        MUJ <- array (muj, dim = c(k, G, no_rmc))
#        WXJ <- array (wxj, dim = c(k, no_rmc))
#        WMUJ <- array (wmuj, dim = c(k, no_rmc))
#        NUJ <- array (nuj, dim = c(k, no_rmc))
#        WX <- array (wx, dim = no_rmc) ## a vector
#        WMU <- array (wmu, dim = no_rmc) ## a vector
#        WNU <- array (wnu, dim = no_rmc) ## a vector

#        ## Gather sufficient statistic
#        tgsum_X <- t(rowsum (X_tr_sel,y_tr)) ## grouped sum of features
#        sum_X2 <- colSums (X_tr_sel^2)


#        ## static variables in Gibbs sampling
#        alpha_wmuj <- (alpha1_mu + G) / 2
#        alpha_wxj <- (alpha1_x + n) / 2
#        alpha_wnu <- (alpha0_nu + k) / 2
#        alpha_wmu <- alpha1_mu * k / 2  - alpha0_mu / 2
#        alpha_wx <- alpha1_x * k / 2  - alpha0_x / 2
#        lambda0_wmu <- alpha0_mu * w0_mu / 2 ##lambda for wmu from prior
#        lambda0_wx <- alpha0_x * w0_x / 2 ##lambda for wx from prior
#        stepsizes_mhwmux <-
#            1 / sqrt ( c (alpha_wmu, alpha_wx) ) / sqrt (8) * stepadj_mhwmux


#        ## Start Marlov chain super-transition
#        time_train [i] <- system.time (
#        for (i_record in 1 : no_rmc)
#        {

#            ## start Gibbs sampling
#            for (i_imc in 1 : no_imc)
#            {
#                ## update muj
#                vars_muj <- 1 / (1/wmuj + outer(1/wxj,nos_g) )
#                means_muj <- (nuj / wmuj + tgsum_X / wxj) * vars_muj
#                muj <- means_muj + replicate (G, rnorm (k)) * sqrt (vars_muj)

#                ## update wxj
#                lambda_wxj <-  ( alpha1_x * wx + sum_X2 -
#                    2 * rowSums (tgsum_X * muj) + rowSums (muj^2 %r% nos_g) )/2
#                wxj <- rinvgam (k, alpha_wxj, lambda_wxj)

#                ## update wmuj
#                lambda_wmuj <- (alpha1_mu * wmu + rowSums ((muj - nuj)^2) ) / 2
#                wmuj <- rinvgam (k, alpha_wmuj, lambda_wmuj)

#                ## update nuj
#                sum_muj <- rowSums (muj)
#                vars_nuj <- 1 / (1 / wnu + G / wmuj)
#                means_nuj <- sum_muj / wmuj * vars_nuj
#                nuj <- means_nuj + rnorm (k) * sqrt (vars_nuj)

#                ## update wnu
#                lambda_wnu <- (alpha0_nu * w0_nu + sum (nuj^2) ) / 2
#                wnu <- rinvgam (1, alpha_wnu, lambda_wnu)

#                ## update wx and wu together with M-H methods
#                ## log posterior of log (wmu, wx)
#                lambda_wmu <- alpha1_mu * sum(1/wmuj) / 2
#                lambda_wx <- alpha1_x * sum(1/wxj) / 2
#                logpost_logwmux <-  function (lw)
#                {
#                    w <- exp (lw)
#                    b4cor <-
#                      alpha_wmu * lw[1] - lambda_wmu * w[1] -
#                      lambda0_wmu / w[1] + alpha_wx * lw[2] -
#                      lambda_wx * w[2] - lambda0_wx / w[2]
#                    
#                    b4cor
#                }

#                log_wmu_wx <- met_gauss (
#                    iters = no_mhwmux, log_f = logpost_logwmux,
#                    ini_value = c(logwmu, logwx), stepsize = stepsizes_mhwmux,
#                    diag_mh = diag_mhwmux )
#                logwmu <- log_wmu_wx [1]
#                logwx <- log_wmu_wx [2]
#                wmu <- exp (logwmu)
#                wx <- exp (logwx)
#            }

#            ## write states into Marlov chain arrays
#            MUJ [,,i_record] <- muj
#            WXJ [,i_record] <- wxj
#            NUJ [,i_record] <- nuj
#            WMUJ [,i_record] <- wmuj
#            WX [i_record] <- wx
#            WMU [i_record] <- wmu
#            WNU [i_record] <- wnu
#        }
#        )[1]
#        ## posterior mean of y
#        post_y <- nos_g + prior_psi
#        freqy <- post_y / sum( post_y)

#        fit_bcbcsf <- list (
#            MUJ =   abind (ini_mcstate $ MUJ, MUJ),
#            WXJ =   abind (ini_mcstate $ WXJ, WXJ),
#            NUJ =   abind (ini_mcstate $ NUJ, NUJ),
#            WMUJ =  abind (ini_mcstate $ WMUJ, WMUJ),
#            WX =    abind (ini_mcstate $ WX, WX),
#            WMU =   abind (ini_mcstate $ WMU, WMU),
#            WNU =   abind (ini_mcstate $ WNU, WNU),
#            freqy = freqy,
#            fsel =  fsel,
#            nuj_std = nuj_std,
#            wxj_std = wxj_std,
#            no_imc =  no_imc,
#            no_rmc = no_rmc,
#            bias_corrected = bcor)

#        fit_bcbcsf_1file <- paste (fit_bcbcsf_prefile,
#        "w0_", w0_mu, "_alpha0_", alpha0_mu, "_alpha1_",alpha1_mu,
#        "_n_",n, "_p_", p, "_nfsel_", k, "_biascorected_", bcor,
#        ".RData", sep = "")

#        save (fit_bcbcsf, file = fit_bcbcsf_1file)

#        fit_bcbcsf_files[[i]] <- fit_bcbcsf_1file

#        ## prediction
#        if (!is.null (X_ts))
#        {
#            time_pred[i] <- system.time (
#            array_probs_pred[,,i] <-
#            mcmcpred (X_ts = X_ts, fit_bcbcsf = fit_bcbcsf,
#            burn = burn, thin = thin)
#            )[1]
#        }
#    }

#    list (
#     nos_fsel = nos_fsel, fit_bcbcsf_files = fit_bcbcsf_files,
#     array_probs_pred = array_probs_pred,
#     time_train = time_train, time_pred = time_pred)
#}

#reload_fit_bcbcsf <- function (fit_bcbcsf_1file)
#{
#    load (fit_bcbcsf_1file)
#    return (fit_bcbcsf)
#}

#avgmc_hb <- function (fit_bcbcsf_1file = NULL, fit_bcbcsf = NULL,
#            burn = NULL, thin = 1)
#{
#    if (is.null(fit_bcbcsf))
#        fit_bcbcsf <- reload_fit_bcbcsf (fit_bcbcsf_1file)

#    if (is.null (burn)) burn <- floor (fit_bcbcsf$no_rmc * 0.2)

#    muj <- apply (fit_bcbcsf $ MUJ[,, - (1:burn)], c(1,2), mean )
#    wxj <- apply (fit_bcbcsf $ WXJ[, - (1:burn)], 1, mean )
#    nuj <- apply (fit_bcbcsf $ NUJ[, - (1:burn)], 1, mean )
#    wmuj <- apply (fit_bcbcsf $ WMUJ[, - (1:burn)], 1, mean )
#    wx <-  mean (fit_bcbcsf $ WX[- (1:burn)] )
#    wmu <-  mean (fit_bcbcsf $ WMU[- (1:burn)] )
#    wnu <-  mean (fit_bcbcsf $ WNU[- (1:burn)] )
#    cmuj <- muj - nuj

#    list (muj = muj, wxj = wxj, wmuj = wmuj, nuj = nuj,
#    wx = wx, wmu = wmu, wnu = wnu, cmuj = cmuj)
#}



