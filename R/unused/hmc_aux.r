log_sum_exp <- function(lx)
{  mlx <- max(lx)
  log(sum(exp(lx - mlx))) + mlx
}

#     log_sum_exp <- function (la)
#     {   m <- max (la)
#         if (m > 0)
#         {
#             ## subtract the maximum value before evaluation
#             log (sum (exp(la - m)) + exp(- m) ) + m
#         }
#         else
#         {
#             ## directly evaluation
#             log (sum(exp(la)) + 1)
#         }
#     }
#


comp_lsl <- function (lv)
{
    apply (cbind (0,lv), 1, log_sum_exp)
}


ldt <- function (md, p, K, nu, lw)
{
    - (nu + K) / 2 * sum (log ( md +  exp (lw) * nu))  + lw * p * nu / 2
}

met_gauss2 <- function (iters = 100, stepsize = 0.5, log_f, ini_value,
             iters_imc = 1,  ...)
{
    state <- ini_value
    no_var <- length (state)
    logf <- log_f (ini_value,...)

    if (!is.finite (logf)) stop ("Initial value has 0 probability")

    one_mc <- function ()
    {
        new_state <- rnorm (no_var, state, stepsize)
        new_logf <- log_f (new_state,...)

        if (log (runif(1)) < new_logf - logf)
        {
            state <<- new_state
            logf <<- new_logf
        }
    }

    one_sup <-  function ()
    {
        replicate (iters_imc, one_mc())
        state
    }

    replicate (iters, one_sup () )
}

## y takes values 1,..., NC
bplr_hmc <- function (
    ## data
    y, X, NC = length (unique (y)),
    ## prior for hyperparameters
    width_b0 = 10, nu = 2, lw0 = 0, w_lw0 = 5,
    ## parameters for hybrid Monte Carlo sampling for deltas
    iters_sup = 200, iters_imc = 1, L_leapfrog = 50, hm_stepsize_adj = 0.3,
    ## parameters for M-H sampling for log w
    iters_lw = 50, lw_stepsize_adj = 1,
    ## parameters for visualizing leapfrog transition and initial_state
    lookleapfrog = FALSE, initial_state = NULL)
{
    if (max (y) > NC | min (y) < 1)  stop ("values of y is out of range")
    y <- y - 1 ## in R program, y takes 0,..., K
    K <- NC - 1

    p <- ncol (X)
    n <- nrow (X)
    no_rej <- 0
    mat_y <- matrix (0, n, K)
    sapply (1:n, function (i) {if (y[i]>0) mat_y[i,y[i]] <<- 1}  )
    newX <- cbind (1, X)
    tnewX <- t (newX)
    dbl_mdiff_loglike <- 1/4 * apply (newX^2, 2, sum)
    sigma_b0 <- width_b0 ^ 2
    ones <- matrix (1, K, 1)

    ## checking appropriateness of arguments
    if (iters_sup <= 0 || iters_imc <=0 ||  L_leapfrog <=0 || iters_lw <=0
        || lw_stepsize_adj <= 0 )
        stop ("Some parameter for Markov chain sampling is negative.")
    if (length (y) != n )
        stop ("Number of cases in 'y' and 'x' NOT match")

    ## Initialize Markov chain parameters
    if (is.null (initial_state))
    {
        deltas <- matrix (0, p + 1, K)
        lw <- lw0
    }
    else
    {
        deltas <- matrix (initial_state [1:((p + 1) * K)], p + 1, K)
        sigmas_bt <- initial_state[((p + 1) * K) + (1:p)]
        lw <- initial_state [ (K + 1) * (p + 1)]
    }

    ## compute derivative of log prior
    lv <- newX %*% deltas
    lsl <- comp_lsl (lv)
    predprob <- exp (lv - lsl)

    sumdeltas <- as.vector (deltas %*% ones)
    sumsqdeltas <- as.vector (deltas ^ 2 %*% ones)
    md <- sumsqdeltas  - sumdeltas^2 / (K + 1)
    mdiff_pp <- deltas - sumdeltas / (K + 1)
    mdiff_loglike <- tnewX %*% (predprob - mat_y)
    loglike <- sum (lv * mat_y) - sum (lsl)

    one_mc <- function()
    {
#         ## update lw with M-H sampling
#         log_post_lw <- function (lw)
#         {
#             ldt (md[-1], p, K, nu, lw) + dnorm (lw, lw0, w_lw0, log = TRUE )
#         }
#         lw <<- met_gauss2 (
#                iters = iters_lw, stepsize = lw_stepsize_adj/sqrt (p),
#                ini_value = lw, log_f = log_post_lw) [iters_lw]
#         w <- exp (lw)
#         ## update sigmas for betas
#         sigmas_bt <<- 1/rgamma (p, shape=(nu + K) / 2) * (nu * w + md [-1]) / 2

        # start leapfrog
        allsigmas <-  c (sigma_b0, sigmas_bt)

        ## set stepsizes for deltas
        stepsizes <- hm_stepsize_adj /
                     sqrt (dbl_mdiff_loglike + K/(K+1) / allsigmas)
        stepsizes <- matrix (stepsizes, p + 1, K)

        ## initialize trajectory state
        momt_trj <- matrix (rnorm ((p + 1) * K), p + 1, K)
        deltas_trj <- deltas
        md_trj <- md
        mdiff_pp_trj <- mdiff_pp
        mdiff_loglike_trj <- mdiff_loglike
        loglike_trj <- loglike

        ## compute initial derivatives
        mdiff_logpost_trj <- mdiff_loglike_trj + mdiff_pp_trj / allsigmas

        ## compute initial energy
        loglike_old <- loglike_trj
        logprior_old <- - sum (md_trj / allsigmas) / 2
        log_prob_momt_old <- - sum (momt_trj^2)/2

        if (lookleapfrog) ## recording initial values
        {
            magdeltas_lp <- sqrt (sum (deltas_trj ^ 2) )
            loglike_lp <- loglike_old
            logprior_lp <- logprior_old
            log_prob_momt_lp <- log_prob_momt_old
        }

        ## start leapfrog transition
        momt_trj <- momt_trj - stepsizes / 2 * mdiff_logpost_trj
        one_trj <- function ()
        {
            ## move deltas
            deltas_trj <<- deltas_trj + stepsizes * momt_trj

            ## compute derivative of minus log likelihood
            lv_trj <<- newX %*% deltas_trj
            lsl_trj <<- comp_lsl (lv_trj)
            predprob_trj <<- exp (lv_trj - lsl_trj)
            mdiff_loglike_trj <<- tnewX %*% (predprob_trj - mat_y)
            ## compute derivative of minus log prior
            sumdeltas_trj <<- as.vector (deltas_trj %*% ones)
			mdiff_pp_trj <<- deltas_trj - sumdeltas_trj/(K+1)
            ## compute derivative of minus log post
            mdiff_logpost_trj <<-  mdiff_loglike_trj + mdiff_pp_trj/allsigmas

            ## move momoton
            momt_trj <<- momt_trj - stepsizes * mdiff_logpost_trj

            if (lookleapfrog)
            {
                magdeltas_lp <<-
                   c(magdeltas_lp, sqrt (sum (deltas_trj^2)))
                loglike_lp <<-
                   c(loglike_lp, sum (lv_trj * mat_y) - sum (lsl_trj))
                sumsqdeltas_trj <- as.vector (deltas_trj ^ 2 %*% ones)
                md_lp <- sumsqdeltas_trj - sumdeltas_trj^2 / (K + 1)
                logprior_lp <<- c(logprior_lp, - sum(md_lp / allsigmas) / 2)
                momt_trj <- momt_trj + stepsizes/2 * mdiff_logpost_trj
                log_prob_momt_lp <<- c(log_prob_momt_lp, - sum (momt_trj^2)/2)
            }

        }
         replicate (L_leapfrog, one_trj())
        ## move back half step in changing momenton
        momt_trj <- momt_trj + stepsizes/2 * mdiff_logpost_trj

        ## decide whether to accept it
        loglike_trj <- sum (lv_trj * mat_y) - sum (lsl_trj)
        sumsqdeltas_trj <- as.vector (deltas_trj ^ 2 %*% ones)
        md_trj <- sumsqdeltas_trj - sumdeltas_trj^2 / (K + 1)
        logprior_trj <- - sum(md_trj/allsigmas)/2
        log_prob_momt_trj <- - sum (momt_trj^2) /2

        if ( log (runif(1)) <=
             loglike_trj + logprior_trj + log_prob_momt_trj -
            (loglike_old + logprior_old + log_prob_momt_old) )
        {
            deltas <<- deltas_trj
            md <<- md_trj
            mdiff_pp <<- mdiff_pp_trj
            loglike <<- loglike_trj
            mdiff_loglike <<- mdiff_loglike_trj
        }
        else no_rej <<- no_rej + 1

        if(lookleapfrog)
        {
            menergy_lp <- logprior_lp + loglike_lp + log_prob_momt_lp
            logpost_lp <- logprior_lp + loglike_lp
            par ( mfrow=c(2,2) )
            plot (magdeltas_lp, main = "Trajectory of Magnitude of Deltas",
                  ylab = "", xlab = "Trajectory Index",type = "l")
            plot (loglike_lp, main = "Trajectory of Log Likelihood",
                  ylab = "", xlab = "Trajectory Index", type = "l")
            plot (logpost_lp, main = "Trajectory of Log Posterior",
                  ylab = "", xlab = "Trajectory Index", type = "l")
            plot (menergy_lp, main = "Trajectory of Minus Total Energy",
                  ylab = "", xlab = "Trajectory Index", type = "l")

            browser()
        }

    }

    one_sup <- function ()
    {
       replicate (iters_imc, one_mc() )
       ## return mcstate after 'iters_imc' MC transitions
       c (deltas, sigmas_bt, lw, loglike)
    }

    mcchain <- replicate (iters_sup, one_sup() )

    attr (mcchain, "rej.rate") <- no_rej / iters_sup / iters_imc
    attr (mcchain, "K") <- K
    attr (mcchain, "p") <- p

    ## return Markov chain
    mcchain
}

read_deltas <- function (mcchain, burn = NULL, thin = 1, showmedian = FALSE)
{
  iters <- ncol (mcchain)
  if (is.null (burn)) burn <- ceiling(iters * 0.2)
  K <- attr (mcchain, "K")
  p <- attr (mcchain, "p")

  indice <- seq (burn + 1, iters, by = thin)
  lmc <- length (indice)

  deltas <- array (mcchain[1:((p+1)*K),  indice], dim = c(p + 1, K, lmc) )

  if (showmedian)
	deltas <- apply (deltas, MARGIN = c(1,2), median)

  deltas

}

read_sigmasbt <- function (mcchain, burn = NULL, thin = 1, log = TRUE,
				showmedian = FALSE)
{
  iters <- ncol (mcchain)
  if (is.null (burn)) burn <- ceiling(iters * 0.2)
  K <- attr (mcchain, "K")
  p <- attr (mcchain, "p")

  indice <- seq (burn + 1, iters, by = thin)

  sigmasbt <- mcchain [((p+1)*K + 1):((p+1)*K + p) ,  indice, drop = FALSE]

  if (log) sigmasbt <- log (sigmasbt)

  if (showmedian) sigmasbt <- apply (sigmasbt, 1, median)

  sigmasbt
}

read_lw <- function (mcchain, burn = NULL, thin = 1)
{
  iters <- ncol (mcchain)
  if (is.null (burn)) burn <- ceiling(iters * 0.2)
  K <- attr (mcchain, "K")
  p <- attr (mcchain, "p")

  indice <- seq (burn + 1, iters, by = thin)

  lw <- mcchain [ (K+1) * (p+1),  indice]
  lw
}

read_llik <- function (mcchain, burn = NULL, thin = 1)
{
  iters <- ncol (mcchain)
  K <- attr (mcchain, "K")
  p <- attr (mcchain, "p")

  if (is.null (burn)) burn <- ceiling(iters * 0.2)

  indice <- seq (burn + 1, iters, by = thin)

  llik <- mcchain [ (K+1) * (p+1) + 1,  indice]
  llik
}

summary_mcbplr <- function (mcchain)
{
  iters <- ncol (mcchain)
  K <- attr (mcchain, "K")
  p <- attr (mcchain, "p")
  rej.rate <- attr (mcchain, "rej.rate")
  c (iters = iters, K = K, p = p, rej.rate = rej.rate)

}

predict_bplr <- function (X_test, mcchain, burn = NULL, thin = NULL)
{
    if (is.vector (X_test)) X_test <- matrix (X_test, 1,)
    K <- attr (mcchain,"K")
    iters <- ncol (mcchain)
    p <- ncol (X_test)

    if (is.null (burn)) burn <- ceiling (iters * 0.2)

    if (is.null (thin)) thin <- floor(iters * 0.8 /100)

    X_test <- cbind (1, X_test)

    indice_pred <- seq (burn + 1, iters, by = thin)

    comp_prob <- function(deltas)
    {
       deltas <- matrix (deltas, p + 1, K)
       lv <- X_test %*% deltas
       as.vector (exp (cbind(0,lv) - comp_lsl(lv)) )
    }

    probs_pred <- matrix (
		apply (
		  apply (mcchain[1:((p+1)*K), indice_pred, drop = FALSE],
		  2, comp_prob),
		1, mean),
		nrow (X_test), K + 1)

    list (probs_pred = probs_pred,
		  values_pred = apply (probs_pred, 1, which.max))
}

## NC -- number of classes
## y_tr -- true values of y, taking values 1, ..., NC  (ie, 0,..., K)
trpr_bplr <- function (X_ts, X_tr, y_tr, NC, no_sel = ncol (X_tr))
{
  features_sel <- rank_k(X_tr,y_tr)[1:no_sel]
  X_tr_sel <- X_tr [ ,features_sel, drop = FALSE]
  X_ts_sel <- X_ts [ ,features_sel, drop = FALSE]
  mcchain <- bplr_hmc (y = y_tr, X = X_tr_sel, NC = NC)
  predict_bplr (X_test = X_ts_sel, NC = NC, mcchain = mcchain)
}

