
comp_lsl <- function (lv)
{
    log_sum_exp <- function (la)
    {   m <- max (la)
        if (m > 0) 
        {   
            ## subtract the maximum value before evaluation
            log (sum (exp(la - m)) + exp(- m) ) + m
        }
        else 
        { 
            ## directly evaluation
            log (sum(exp(la)) + 1)
        }
    }
    
    apply (lv, 1, log_sum_exp)
}
 
 
ldt <- function (md, p, K, nu, lw)
{
    (- (nu + K) * sum (log ( 1 + md / exp (2*lw) / nu)) ) / 2  - lw * p * K
} 

met_gauss <- function (iters = 100, stepsize = 0.5, log_f, ini_value, 
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
		
		accepted <- (log (runif(1)) < new_logf - logf)
		if (accepted)
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

bplr_hmc <- function (
    ## data 
    y, X, K = length (unique (y)) - 1,
    ## prior for hyperparameters
    width_b0 = 10, nu = 2, lw0 = 0, w_lw0 = 0.5, 
    ## parameters for hybrid Monte Carlo sampling for deltas
    iters_sup = 100, iters_imc = 1, L_leapfrog = 100, hm_stepsize_adj = 0.3,
    ## parameters for M-H sampling for log wsq
    iters_lw = 50, lw_stepsize_adj = 0.5, 
    ## parameters for visualizing leapfrog transition and initial_state
    lookleapfrog = FALSE, initial_state = NULL)
{
    p <- ncol (X)
    n <- nrow (X)
    no_rej <- 0 
    mat_y <- matrix (0, n, K)
    sapply (1:n, function (i) {if (y[i]>0) mat_y[i,y[i]] <<- 1}  )
    newX <- cbind (1, X)
    tnewX <- t (newX)
    dbl_mdiff_loglike <- 1/4 * apply (newX^2, 2, sum)
    sigma_b0 <- width_b0 ^ 2

    ## checking appropriateness of arguments
    if (iters_sup <= 0 || iters_imc <=0 ||  L_leapfrog <=0 || iters_lw <=0 
        || lw_stepsize_adj <= 0 )
        stop ("Some parameter for Markov chain sampling is negative.")
    if (length (y) != n )  
        stop ("Number of cases in 'y' and 'x' NOT match")
    
    ## Initialize Markov chain parameters
    if (is.null (initial_state))
    {
        deltas <- replicate (K, rnorm (p + 1)) 
        widths_bt <- rep (1,p)
        lw <- lw0
    }
    else
    {
        deltas <- matrix (initial_state [1:((p + 1) * K)], p + 1, K)
        widths_bt <- initial_state [(p + 1) * K + (1 : p)]
        lw <- initial_state [ (K + 1) * (p + 1)]
    }

    wsq <- exp (2 * lw)
    
    ## compute derivative of log prior
	sumdeltas <- apply (deltas, 1, sum)
    sumsqdeltas <- apply (deltas ^ 2, 1, sum)
    md <- sumsqdeltas  - sumdeltas^2 / (K + 1) 
    psuedo_sigmas <-  c (sigma_b0, (nu * wsq + md[-1]) / (nu + K) )
    mdiff_logprior <- (deltas - sumdeltas/(K + 1)) / psuedo_sigmas
    ## compute derivative of log likelihood
    lv <- newX %*% deltas
    lsl <- comp_lsl (lv)
    predprob <- exp (lv - lsl)
    mdiff_loglike <- tnewX %*% (predprob - mat_y) 

    mdiff_logpost <- mdiff_logprior + mdiff_loglike

    ## compute log prior
    logprior <- - md[1] / sigma_b0 / 2 + ldt (md[-1], p, K, nu, lw) 
    ## compute log likelihood
    loglike <- sum (lv * mat_y) - sum (lsl)
    logpost <- logprior + loglike

    ## update sigmas for betas
    widths_bt <-  sqrt ( 1 / rgamma (p, shape = (nu + K) / 2) *  
                         (nu * wsq + md [-1]) / 2 )

    one_mc <- function()
    {   
        ## set stepsizes for deltas 
#        stepsizes <- hm_stepsize_adj /  sqrt (dbl_mdiff_loglike + 
#                     K/(K+1)* c (1/sigma_b0, rep ((nu+K)/nu/wsq, p)))
		stepsizes <- hm_stepsize_adj / sqrt (dbl_mdiff_loglike + 
                     K/(K+1) / c (sigma_b0, widths_bt^2))

        stepsizes <- replicate (K, stepsizes)

        ## initialize position variables, and remember previous derivative
        deltas_trj <- deltas
        mdiff_logpost_trj <- mdiff_logpost

        ## initialize momenton variables
        momt_trj <- replicate (K, rnorm (p + 1))
        log_prob_momt <- sum (dnorm (momt_trj, log = TRUE) )

        if (lookleapfrog) ## recording original state
        {
            magdeltas_lp <- sqrt (sum (deltas ^ 2) )
            loglike_lp <- loglike
            menergy_lp <- loglike + logprior + log_prob_momt
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
        	sumdeltas_trj <<- apply (deltas_trj, 1, sum)
            sumsqdeltas_trj <<- apply (deltas_trj ^ 2, 1, sum)
            md_trj <<- sumsqdeltas_trj - sumdeltas_trj^2 / (K + 1) 
            psuedo_sigmas <- c (sigma_b0, (nu * wsq + md_trj[-1]) / (nu + K))
            mdiff_logprior_trj <<- 
                (deltas_trj - sumdeltas_trj/(K+1)) / psuedo_sigmas

            ## compute derivative of minus log post 
            mdiff_logpost_trj <<-  mdiff_loglike_trj + mdiff_logprior_trj
            
            momt_trj <<- momt_trj - stepsizes * mdiff_logpost_trj
            
            if (lookleapfrog)
            {
                    magdeltas_lp <<- 
                   c (magdeltas_lp, sqrt (sum (deltas_trj^2)))
                loglike_lp <<- 
                   c (loglike_lp, sum (lv_trj * mat_y) - sum (lsl_trj))
                logprior_lp <<- 
                   - md_trj[1] / sigma_b0 / 2 + ldt (md_trj[-1], p, K, nu, lw)
                log_prob_momt_lp <- sum (dnorm(momt_trj, log = TRUE) ) 
                menergy_lp <<- 
                   c (menergy_lp, loglike_lp + logprior_lp + log_prob_momt_lp)
            }
        }        
        replicate (L_leapfrog, one_trj())
        ## move back half step change of momenton
        momt_trj <- momt_trj + stepsizes/2 * mdiff_logpost_trj 
        
        ## decide whether to accept it
        loglike_trj <- sum (lv_trj * mat_y) - sum (lsl_trj)
        logprior_trj <- - md_trj[1]/sigma_b0/2 + ldt (md_trj[-1], p, K, nu, lw) 
        log_prob_momt_trj <- sum (dnorm(momt_trj, log = TRUE))

        if (log (runif(1)) <= loglike_trj + logprior_trj + log_prob_momt_trj - 
            (loglike + logprior + log_prob_momt ) )
        {
            deltas <<- deltas_trj
            loglike <<- loglike_trj
            logprior <<- logprior_trj
            mdiff_logpost <<- mdiff_logpost_trj
        }
        else no_rej <<- no_rej + 1

        if(lookleapfrog)
        {
            par (mfrow=c(1,3))
            plot (magdeltas_lp, main = "Magnitude of Deltas", xlab = "")
            plot (loglike_lp, main = "Log Likelihood", xlab = "")
            plot (menergy_lp, main = "Minus Total Energy", xlab = "")
            browser()
        }

        ## update lw with M-H sampling
        log_post_lw <- function (lw)
        {   
            ldt (md[-1], p, K, nu, lw) + dnorm (lw, lw0, w_lw0, log = TRUE )        
        }
        
        lw <<- met_gauss (
               iters = iters_lw, stepsize = lw_stepsize_adj/sqrt (p), 
        	   ini_value = lw, log_f = log_post_lw) [iters_lw]
        wsq <<- exp (2 * lw)	
        
        ## update sigmas for betas
        widths_bt <<-  sqrt ( 1 / rgamma (p, shape = (nu + K) / 2) *  
                              (nu * wsq + md [-1]) / 2 )
        
    }
    
    one_sup <- function () 
    {  
       replicate (iters_imc, one_mc() )
       ## return mcstate after 'iters_imc' MC transitions
       c (deltas, widths_bt, lw, loglike) 
    }

    mcchain <- replicate (iters_sup, one_sup() )
    attr(mcchain, "rej.rate") <- no_rej / iters_sup / iters_imc
    
    ## return Markov chain
    mcchain 
}


#        if (nu*p - nu0 > 0)
#        {
#            sum_tau <- sum (1/widths_bt)
#            accepted <- FALSE
#            no_rejs <- -1
#            while (!accepted)
#            {   no_rejs <- no_rejs + 1
#                if (no_rejs > 50)
#                   stop (
#  "More than 50 rejections in sampling 'sigma'. w0/nu0 may be too big.")
#                sigma <<- 
#                    rgamma (1, (nu * p  - nu0) / 2, nu * sum_tau / 2)
#                accepted <- 
#                    (log (runif(1)) <= - nu0 * w0/2/sigma) 
#            }
#        }

