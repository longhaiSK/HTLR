
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

p <- 100
K <- 2 
nu1 <- 2
width <- 0.1

sigmas <- 1 / rgamma (p, nu1/2, nu1 * width^2 / 2)
betas <- replicate (K + 1, rnorm(p)) * sqrt (sigmas)
deltas <- betas [, -1] - betas [, 1]
sumdeltas <- apply (deltas, 1, sum)
sqsumdeltas <- sumdeltas ^ 2
sumsqdeltas <- apply (deltas ^ 2, 1, sum)
md <- sumsqdeltas - sqsumdeltas / (K + 1) 

ldt <- function (md, p, K, nu, lw)
{
    (- (nu + K) * sum (log ( 1 + md / exp (2*lw) / nu)) ) / 2 - lw * p * K
} 

log_post_lw <- function (lw)
{   
    ldt (md, p, K, nu1, lw) +  dnorm (lw, 0, 5, log = TRUE )        
}

lw <-  met_gauss (iters = 1000, stepsize = 1/sqrt(p), 
                  ini_value = log (width), log_f = log_post_lw) 

plot (lw)
abline (h = log (width))



