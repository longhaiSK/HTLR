
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
    ## parameters for M-H sampling for log w2
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
        deltas <- replicate (K, rnorm (p + 1)) 
#        sigmas_bt <- rep (1,p)
        lw <- lw0
    }
    else
    {
        deltas <- matrix (initial_state [1:((p + 1) * K)], p + 1, K)
#        sigmas_bt <- initial_state [(p + 1) * K + (1 : p)]
        lw <- initial_state [ (K + 1) * (p + 1)]
    }

    ## compute derivative of log prior
	sumdeltas <- as.vector (deltas %*% ones)
    sumsqdeltas <- as.vector (deltas ^ 2 %*% ones)
    md <- sumsqdeltas  - sumdeltas^2 / (K + 1) 
    lv <- newX %*% deltas
    lsl <- comp_lsl (lv)
    predprob <- exp (lv - lsl)
    mdiff_loglike <- tnewX %*% (predprob - mat_y) 
    loglike <- sum (lv * mat_y) - sum (lsl) 

    one_mc <- function()
    {   
        ## update lw with M-H sampling
        log_post_lw <- function (lw)
        {   
            ldt (md[-1], p, K, nu, lw) + dnorm (lw, lw0, w_lw0, log = TRUE )        
        }
        
        lw <<- met_gauss (
               iters = iters_lw, stepsize = lw_stepsize_adj/sqrt (p), 
        	   ini_value = lw, log_f = log_post_lw) [iters_lw]
        w2 <- exp (2 * lw)	
        
        ## update sigmas for betas
        sigmas_bt <<- 1/rgamma (p, shape=(nu + K)/2) * (nu * w2 + md [-1]) / 2 
    
		## start leapfrog
        allsigmas <-  c (sigma_b0, sigmas_bt)
		
        ## set stepsizes for deltas 
		stepsizes <- hm_stepsize_adj /
					 sqrt (dbl_mdiff_loglike + K/(K+1) / allsigmas)
        stepsizes <- matrix (stepsizes, p + 1, K)

        ## initialize trajectory state
        momt_trj <- matrix (rnorm ((p + 1)*K), p + 1, K)
        deltas_trj <- deltas
        sumdeltas_trj <- sumdeltas
        md_trj <- md
        mdiff_loglike_trj <- mdiff_loglike

        mdiff_logprior_trj <- (deltas_trj - sumdeltas_trj/(K + 1)) / allsigmas
        mdiff_logpost_trj <- mdiff_loglike_trj + mdiff_logprior_trj 
        
        loglike_old <- loglike 
        logprior_old <- - sum (md_trj / allsigmas) / 2
        log_prob_momt_old <- - sum (momt_trj^2)/2

        if (lookleapfrog) ## recording original state
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
            mdiff_logprior_trj <<- (deltas_trj-sumdeltas_trj/(K+1))/allsigmas
            ## compute derivative of minus log post 
            mdiff_logpost_trj <<-  mdiff_loglike_trj + mdiff_logprior_trj
            
            momt_trj <<- momt_trj - stepsizes * mdiff_logpost_trj
            
            if (lookleapfrog)
            {
                magdeltas_lp <<- 
                   c (magdeltas_lp, sqrt (sum (deltas_trj^2)))
                loglike_lp <<- 
                   c (loglike_lp, sum (lv_trj * mat_y) - sum (lsl_trj) )
                logprior_lp <<- 
                   c (logprior_lp, - sum(md_trj / allsigmas) / 2)
                log_prob_momt_lp <- 
                   c (log_prob_momt_lp, - sum (momt_trj^2)/2) 
            }

        }        
        replicate (L_leapfrog, one_trj())
        ## move back half step in changing momenton
        momt_trj <- momt_trj + stepsizes/2 * mdiff_logpost_trj 

        if(lookleapfrog)
        {
            par (mfrow=c(1,3))
            plot (magdeltas_lp, main = "Magnitude of Deltas", xlab = "")
            plot (loglike_lp, main = "Log Likelihood", xlab = "")
            plot (menergy_lp, main = "Minus Total Energy", xlab = "")
            browser()
        }
        
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
            sumdeltas <<- sumdeltas_trj
            md <<- md_trj
            loglike <<- loglike_trj
            mdiff_loglike <<- mdiff_loglike_trj
        }
        else no_rej <<- no_rej + 1
    }
    
    one_sup <- function () 
    {  
       replicate (iters_imc, one_mc() )
       ## return mcstate after 'iters_imc' MC transitions
       c (deltas, sigmas_bt, lw, loglike) 
    }

    mcchain <- replicate (iters_sup, one_sup() )
    attr(mcchain, "rej.rate") <- no_rej / iters_sup / iters_imc
    
    ## return Markov chain
    mcchain 
}





