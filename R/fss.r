# this function find mode for each markov chain iterations
# the results will be used to divide MCMC into subpools
# it will be very slow when p is very large
# it is recommended to apply it to fitting results with a pre-selected subset 
#htlr_mode <- function (fithtlr)
#{
#	if (fithtlr$K > 1) 
#		stop ("Currently this function is used for binary logistic regression")
#	fithtlr_mode <- fithtlr
#	iters <- fithtlr$iters_rmc
#	X <- fithtlr$X[, -1, drop = F]
#	y <- fithtlr$ybase
#	prior.scale.int <- sqrt(fithtlr$sigmab0)
#	prior.scale <- exp (fithtlr$s/2) * sqrt (2)
#	prior.df <- fithtlr$alpha
#	
#	for (i in 1:iters)
#	{
#		fithtlr_mode$mcdeltas[,1,i] <- 
#		bayesglm ( y ~ X, family = binomial (link = "logit"),
#		prior.scale = prior.scale,
#		prior.df = prior.df, prior.scale.for.intercept = prior.scale.int,
#		scaled = F, start = fithtlr$mcdeltas[,1,i])$coefficients
#	}
#	
#		
#}


htlr_map <- function (fithtlr, usedmc = NULL, tol = 1e-2, debug = FALSE)
{
	if (is.null (usedmc)) usedmc <- 2:ncol (fithtlr$mcvardeltas)
    alpha <- fithtlr$alpha 
    s <- fithtlr$s
    sigmab0 <- fithtlr$sigmab0
    X_addint <- fithtlr$X
    y_tr <- fithtlr$ybase
    mcvardeltas <- fithtlr$mcvardeltas 
    
    ddim <- dim (fithtlr$mcdeltas)[1:2]
    
   	fithtlr_mode <- fithtlr; fithtlr_mode$mode <- TRUE
   	
   	finish <- 0
   	pb <- utils::txtProgressBar(min = 0,  max = length (usedmc), style = 3)
   	
   	
   	for (i in 1:length (usedmc)){
	   	fithtlr_mode$mcdeltas[,,i] <- optim_logit(
    			array(fithtlr$mcdeltas[,,usedmc[i]], dim = ddim),
		    X_addint = X_addint, y_tr = y_tr, s = s, alpha = alpha, 
    			sigmab0 = sigmab0, tol = tol, debug = debug)
    		finish <- finish + 1
    		utils::setTxtProgressBar(pb, finish)
 	}
 	cat ("\n")
    fithtlr_mode
}

## y_tr is coded by 0, ..., K (=C-1)
## X_addint must be a matrix with rows for cases, columns for variables
optim_logit <- function (deltas, X_addint, y_tr, alpha, s, sigmab0, tol = 1e-2,
					     debug = FALSE)
{
    if (!is.matrix (deltas)) stop ("'deltas' in 'optim_logit' must be a matrix")
	# screening with relative sdbs first 
    mcsdb <- comp_sdb (deltas,removeint = T)
    update <- c(0, which (mcsdb/max(mcsdb) > 0.01)) + 1
	deltas [-update,] <- 0
    
    nvar <- ncol (X_addint)
    K <- ncol (deltas)
    log_aw <- log (alpha) + s
    
    post_logit_avar <- function (d_uvar, uvar, lv_fix)
    {
       d_uvar <- matrix (d_uvar, length (uvar), K)

       ## compute log likelihood
       loglike <- 0
       lv <- lv_fix + X_addint [, uvar, drop = FALSE] %*% d_uvar
       for (k in 1:K) loglike <- loglike + sum (lv[y_tr==k,k])
       loglike <- loglike - log_normcons (lv)
       
       ## compute log prior
       vardeltas <- comp_vardeltas (d_uvar)
       j_intc <- which (uvar == 1)
       j_vars <- which (uvar != 1)
       
       logprior <- 0
       if (length (j_intc) > 0) 
       {
         logprior <- logprior-vardeltas[j_intc]/2/sigmab0
       }
       if (length (j_vars) > 0) 
       {
           logprior <- logprior -  (alpha + K ) / 2 * 
             sum (log_sum_exp(cbind(log_aw, log(vardeltas[j_vars])))) 
       }
       ## return negative log posterior
       - (loglike + logprior)
    }

    lv <- X_addint %*% deltas
    
    ## co-ordinate optimization if number of variables > 10
	while (length(update) > 5)
	{   
		if (debug) cat ("uo:", update, "\n")

		deltas_old <- deltas [update,,drop = F]

		for (uvar in update)
		{
		    d_uvar <- deltas [uvar,,drop = FALSE]
		    lv_fix <- lv - X_addint [,uvar,drop = FALSE] %*% d_uvar
		    
		    out_nlm <- nlm (f = post_logit_avar, p = d_uvar, gradtol = tol,
		                    uvar = uvar, lv_fix = lv_fix)
		    deltas[uvar, ] <- out_nlm$estimate                       
		    d_uvar <- deltas [uvar,,drop = FALSE]
		    lv <- lv_fix + X_addint [, uvar,drop = FALSE] %*% d_uvar
		}
		
		deltas_new <- deltas[update,,drop = F]
		diff_deltas <- apply(abs(deltas_new-deltas_old)/abs(deltas_old),1,max)
		abs_deltas <- apply (abs (deltas_new), 1, max)
		update <- update [diff_deltas > tol & abs_deltas > tol]

	}
	
	## multivariate optimization
	if (length (update) > 0) 
	{
		if (debug) cat ("mo:", update, "\n")
		uvar <- update
		d_uvar <- deltas [uvar,,drop = FALSE]
		lv_fix <- lv - X_addint [, uvar, drop = FALSE] %*% d_uvar            
		out_nlm <- nlm (f = post_logit_avar, p = d_uvar, gradtol = tol,
				    uvar = uvar, lv_fix = lv_fix)
		deltas[uvar, ] <- out_nlm$estimate
    }
    
    deltas
}


# Feature subsets selection with MCMC Samples
# 
# This function divides all of the Markov chain samples into a certain number of subpools, 
# each representing a feature subset.
# 
# The fitting results of \code{htlr_fitpred} are a mixture of subpools (posterior modes) of 
# Markov chain samples for different feature subsets. \code{htlr_fss} divides the Markov chain samples into 
# different subpools associated with different feature subsets (posterior modes).
# 
# @param fithtlr A list containing fitting results by \code{htlr_fitpred}.
# 
# @return 
# \itemize{
#   \item ftab - The rows of \code{ftab} are different feature subsets, with the 1st column showing 
#   the indices of feature subset, the 2nd column showing the relative frequency of Markov chain 
#   iterations in the subpool, and the 3rd column showing the coefficients estimated by median with the subpool.
#   \item mcids - A list of vectors containing the indices of Markov chain iterations associated with 
#   different feature subsets.
#   \item wsdbs - Weighted average of sdbs in all feature subsets found. 
#   The weights are the relative frequencies of Markov chain samples for different feature subsets.
#   \item msdbs - The maximum of sdbs in all feature subsets, for discovering features used in feature subsets with small frequencies.
# }
# 
# @export
# 
# @seealso htlr_fitpred 
htlr_fss <- function  (fithtlr, threshold = NULL, mfreq = NULL, print = TRUE)
{
    if (is.null (threshold) | is.null (mfreq)) {
        if (!is.null (fithtlr$mode) & fithtlr$mode == TRUE) {
            threshold <- 0.01
            mfreq <- 0.01
        }
        else {
            threshold <- 0.05
            mfreq <- 0.05
        }
    }
    
    mcsdb <- apply (fithtlr$mcdeltas,3,comp_sdb)
    mcsdb_max <- apply (mcsdb, 2, max)
    mcsdb_norm <- sweep (mcsdb, 2, mcsdb_max, "/")
    
    is_used <- 1 * (mcsdb_norm >= threshold)
    
    ## selecting features using marginal probabilities
    mprobs <- apply (is_used, 1, mean)
    fset_mf <- which (mprobs >= mfreq)
    
    if (length (fset_mf) == 0) stop ("mfreq may be too big; no features used")
    
    is_used_mf <- is_used [fset_mf,,drop = F]
    mciters <- ncol (mcsdb)
    
    ## find feature subsets
    iclust <- rep (1, mciters)
    clusts <- is_used_mf[,1,drop = FALSE]
    nclust <- 1
    nos_clust <- 1
    for (imc in 2:mciters)
    {
        new_cls <- TRUE
        i <- 2
        while (i <= nclust)
        {
            if ( all (is_used_mf[,imc] == clusts[,i]) )
            {
                new_cls <- FALSE
                iclust [imc] <- i
                nos_clust [i] <- nos_clust [i] + 1
                break
            }
            i <- i + 1
        }
        
        if (new_cls == TRUE) ## found a new cluster
        {
            nclust <- nclust + 1
            iclust [imc] <- nclust
            nos_clust[nclust] <- 1
            clusts <- cbind (clusts, is_used_mf[,imc])
        }
        
    } 
    
    ## find ranks of clusters by subset frequencies
    rank_clust <- rank (-nos_clust, ties.method = "first" )
    iclust  <- rank_clust[iclust]
    clust_order <- order (rank_clust)
    
    ## re-order clusters with ranks
    clusts <- clusts[, clust_order, drop = FALSE]
    nos_clust <- nos_clust [clust_order]
    freqs <- nos_clust/mciters
    cfreqs <- cumsum (freqs)
    
    
    ## summarize feature subsets
    nsubsets <- nclust
    fsubsets <- rep (list (""), nsubsets)
    mcids <- rep (list (""), nsubsets)
    coefs <- rep (list (""), nsubsets)
    coefs_sel <- rep (list (""), nsubsets)   
    sdbs <- rep (list (""), nsubsets)
    ftab <- data.frame (matrix (0, nsubsets, 3))
    colnames (ftab) <-  c("fsubsets", "freqs", "int.&coefs")
    
    for (i in 1:nsubsets)
    {
        ## indices of variables selected in this subset
        ix_fsubsets <- fset_mf [which (clusts[,i] == 1)]
        
        fsubsets[[i]] <- fithtlr$fsel[ix_fsubsets] ## original indices
        mcids [[i]] <- which (iclust == i) ## markov chain indices
        ## coefs of all features
        coefs[[i]] <- htlr_mdcoef (fithtlr, usedmc = mcids[[i]])
        ## coefs for selected features        
        coefs_sel [[i]] <- coefs [[i]][c(0,ix_fsubsets) + 1,,drop = FALSE] 
        ## features importance indices 
        sdbs [[i]] <- comp_sdb (coefs [[i]])        
        ## reporting tables
        ftab[i, 1] <- paste(fsubsets[[i]], collapse = ",")
        ftab[i, 2] <- round(freqs[i],3)
        ftab[i, 3] <- paste (round(coefs_sel[[i]],2), collapse = ",")
    }
    
    ## find weighted sdbs
    mat_sdbs <- matrix(unlist (sdbs), ncol = nsubsets)
    wsdbs <- rowSums(sweep (mat_sdbs, 2, freqs, "*") )
    ## find max sdbs         
    msdbs <- apply(mat_sdbs, 1, max)
    ## find averaged sdbs
    asdbs <- apply (sqrt(fithtlr$mcvardeltas[-1,-1]),1,mean)/sqrt(fithtlr$K + 1)

    fssout <- list (
                ftab = ftab, nsubsets = nsubsets, mcids = mcids, 
                fsubsets = fsubsets, coefs = coefs, iclust = iclust,
                coefs_sel = coefs_sel, freqs = freqs,cfreqs = cfreqs, 
                sdbs = sdbs, wsdbs = wsdbs, msdbs = msdbs,  asdbs = asdbs )
            
    if (print)  print_fss (fssout, sfreq = 0.01)
    
    fssout <- fssout
}

print_fss <- function (fss, sfreq = 0.01)
{
	nsubsets_sel <- max(length(which(fss$freqs >= sfreq)),1)
    freqs <- fss$freqs [1:nsubsets_sel]
    cfreq <- fss$cfreqs [nsubsets_sel]

    print (fss$ftab[1:nsubsets_sel,])
    cat (sprintf(
    "\n%d subsets with freq >= %.2f found from %.0f%% MC Iters\n",
     nsubsets_sel, sfreq, fss$cfreqs [nsubsets_sel] * 100))
     
    bfsel <- sort(unique (unlist (fss$fsubsets[1:nsubsets_sel]))) 
     
    cat (sprintf("%d features used in the above %d subsets are:\n",
         length (bfsel), nsubsets_sel) )
    
    cat(paste (bfsel, collapse = ","), "\n")
    
}


htlr_mapfss <- function (fithtlr, threshold = 0.01, mfreq = 0.01, 
						 tol = 0.01, print = TRUE, ...)
{
	fithtlr_mode <- htlr_map (fithtlr, tol = tol, ...)
	fssout <- htlr_fss (fithtlr_mode,  threshold = threshold, mfreq = mfreq, 
						print = print)

}

# loocv_fss <- function (fss, X, y, sfreq = 0.01)
# {
# 	nsubsets_sel <- max(length(which(fss$freqs >= sfreq)),1)
#     freqs <- fss$freqs [1:nsubsets_sel]
#     cfreq <- fss$cfreqs [nsubsets_sel]
#     cveval <- cv_table (fss$fsubsets[1:nsubsets_sel], x_total = X, y_total = y - 1)
#     cbind (fss$ftab[1:nsubsets_sel,], cveval[,-1])
# 
# }

