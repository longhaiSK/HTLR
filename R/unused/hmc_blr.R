blr_hmc <- function (
    ## data and prior for hyperparameters
    y, X, sigma_alpha = 100, df1 = 3, df0 = Inf, w0 = 0.1, 
    ## parameters for hybrid Monte Carlo
    iters_sup = 100, iters_imc = 1, L_leapfrog = 100, iters_sig, 
    stepsize_adj = 0.3,  lookleapfrog = FALSE, initial_state = NULL)
{
    if (iters_sup <= 0 || iters_imc <=0 ||  L_leapfrog <=0 )
       stop ("Some parameter for Markov chain sampling is negative.")
    
    col_X <- ncol (X)
    row_X <- nrow (X)
    if (length (y) != row_X ) 
        stop ("number of cases in y and X not match")
    no_rej <- 0 
   
    ## function for computing log priors
    comp_log_prior <- function (alpbeta)
    {
        sum (dnorm (alpbeta, 0, sqrt(c(sigma_alpha, sigmas_betas)), 
        log = TRUE))
    }

    ## add a column of 1 to X
    X <- cbind (1, X)
    X_t <- t (X)

    dbl_der_log_like <- 1/4 * apply (X^2, 2, sum)
    
    ## Initialize Markov chain state
    if (is.null(initial_state)){
        alpbeta <- rnorm(col_X+1)
        sigmas_betas <- rep(1,col_X)
        width_sigmas <- w0
    }
    else
    {
        alpbeta <- initial_state [1:(col_X+1)]
        sigmas_betas <- initial_state [col_X+1+(1:col_X)]
        width_sigmas <- initial_state [2*col_X+2]
    }
    
    lv <- X %*% alpbeta
    lsl <- comp_lsl (lv)
    prob1 <- exp (lv - lsl)
    der_log_like <- X_t %*% (prob1 - y) 
    log_like <- sum(y * lv - lsl)
 
    one_mc <- function()
    {   
        ## hybrid Monte Carlo for alpha and betas
        ## set stepsizes
        stepsizes <- stepsize_adj / sqrt (
          dbl_der_log_like + 1 / c(sigma_alpha,sigmas_betas) )

        ## set trajectory state
        alpbeta_trj <- alpbeta
        momvar_trj <- rnorm (col_X + 1)
        der_log_like_trj <- der_log_like
        grd_alpbeta_trj <- der_log_like_trj +
            alpbeta_trj / c(sigma_alpha,sigmas_betas)
        log_prob_momvar_old <- sum( dnorm(momvar_trj,log = TRUE) )
        
        ## start leapfrog
        momvar_trj <- momvar_trj - stepsizes/2 * grd_alpbeta_trj
        
        if(lookleapfrog)
        {
            magbetas_leap <- sqrt (sum (alpbeta_trj^2) )
            loglike_leap <- log_like
            minusenergy_leap <- 
                log_like + comp_log_prior (alpbeta_trj) + 
                sum (dnorm(momvar_trj, log = TRUE) )
        }
        
        one_trj <- function ()
        {   
            alpbeta_trj <<- alpbeta_trj + stepsizes * momvar_trj
            lv_trj <<- X %*% alpbeta_trj
            lsl_trj <<- comp_lsl (lv_trj)
            prob1_trj <<- exp (lv_trj - lsl_trj)
            der_log_like_trj <<- X_t %*% (prob1_trj - y) 
            grd_alpbeta_trj <<-  der_log_like_trj + 
                alpbeta_trj/c(sigma_alpha,sigmas_betas)
            
            momvar_trj <<- momvar_trj - stepsizes * grd_alpbeta_trj
            
            if (lookleapfrog){
                magbetas_leap <<- 
                   c (magbetas_leap,sqrt (sum (alpbeta_trj^2) ))
                loglike_leap <<- 
                   c (loglike_leap, sum (y * lv_trj - lsl_trj))
                minusenergy_leap <<- c(minusenergy_leap, 
                    sum (y * lv_trj - lsl_trj) + 
                    comp_log_prior (alpbeta_trj) + 
                    sum (dnorm(momvar_trj, log = TRUE) ) )
            }
        }
        
        replicate (L_leapfrog, one_trj())
        momvar_trj <- momvar_trj + stepsizes/2 * grd_alpbeta_trj
        
        ## decide whether to accept it
        log_like_new <- sum (y * lv_trj - lsl_trj)
        trj_accepted <- ( log (runif(1)) <= 
          log_like_new + comp_log_prior (alpbeta_trj) + 
          sum (dnorm(momvar_trj, log = TRUE)) - 
          (log_like + comp_log_prior (alpbeta) + log_prob_momvar_old ) )
        
        if(lookleapfrog)
        {
            par (mfrow=c(1,3))
            plot (magbetas_leap)
            plot (loglike_leap)
            plot (minusenergy_leap)
            browser()
        }
        
        if (trj_accepted) 
        {
            log_like <<- log_like_new
            alpbeta <<- alpbeta_trj
            der_log_like <<- der_log_like_trj
        }
        else no_rej <<- no_rej + 1
 
        update_hp <- function ()
        {       
            ## update sigmas for betas
            sigmas_betas <<- 
                1 / rgamma (col_X, shape = (df1 + 1) / 2, 1) * 
                (df1 * width_sigmas + alpbeta[-1]^2) / 2   
            
            ## update width for sigmas with rejection sampling
            if (df1*col_X - df0 > 0)
            {
                sum_tau <- sum (1/sigmas_betas)
                accepted <- FALSE
                no_rejs <- -1
                while (!accepted)
                {   no_rejs <- no_rejs + 1
                    if (no_rejs > 50)
                       stop ( paste("More than 50 rejections in", 
                              "sampling 'width_sigmas'. w0/df0 may be too big.") )
                    width_sigmas <<- 
                        rgamma (1, (df1 * col_X  - df0) / 2, df1 * sum_tau / 2)
                    accepted <- 
                        (log (runif(1)) <= - df0 * w0/2/width_sigmas) 
                }
            }
        }
        replicate (iters_sig, update_hp())
    }
    
    one_sup <- function () ## a super transition
    {  
       replicate (iters_imc, one_mc() )
       
       ## returen mcstate after 'iters_imc' MC transitions
       c (alpbeta, sigmas_betas, width_sigmas, log_like) 
    }

    mcchain <- replicate (iters_sup, one_sup() )
    attr(mcchain, "rej.rate") <- no_rej / iters_sup / iters_imc
    
    ## return Markov chain
    mcchain 
}

