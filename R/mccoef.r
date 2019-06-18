

## this function find summary of deltas over markov chain
htlr_mdcoef <- function (fithtlr, usedmc = "all", features = "all", method = median)
{
    mcdims <- dim (fithtlr$mcdeltas)
    p <- mcdims [1] - 1
    K <- mcdims [2]
    no_mcspl <- mcdims[3]
    
    if (usedmc [1] == "all")  usedmc <- 2:no_mcspl
    
    
    if (is.null(features) || length (features) == 0) ix.f <- c()
    else if (features [1] == "all")  ix.f <- 1:p
    else        
    {
        ix.f <- get_ix (features, fithtlr$fsel)
    }
    
    mcdeltas <- fithtlr$mcdeltas[c(1,ix.f + 1),, usedmc, drop = FALSE]
    
    mddeltas <- apply (mcdeltas, MARGIN = c(1,2), FUN = method)
    
    mddeltas
    
}


## function for plotting coefficients
htlr_mccoef <-  function (
    fithtlr, features = 1, class = 2,  usedmc = "all", 
    symlim = FALSE, drawq = c(0,1), truedeltas = NULL)
{
    mcdims <- dim (fithtlr$mcdeltas)
    p <- mcdims [1] - 1
    K <- mcdims [2]
    no_mcspl <- mcdims[3]
    features <- features [!is.na (features)]
    if (usedmc [1] == "all") usedmc <- 2:no_mcspl
    
    if (length (features) == 1)
    {
        if (features == 0) j <- 1
        else  j <- which(fithtlr$fsel == features) + 1
        k <- class - 1
        deltas <- fithtlr$mcdeltas[j, k,usedmc, drop = FALSE]
        
        plot (deltas, pch = 20, 
              xlab = "Markov Chain Index (after burning and thinning)", 
              ylab = sprintf ("Coef. Value of Feature %d",  features),
              main = sprintf("MC Coefficients for Feature %d (Class %d)",
                             features, class)
        )
        qdeltas <- quantile (deltas, probs = c(0.025,0.5,0.975))
        abline (h = qdeltas[2], lwd = 2)
        abline (h = qdeltas[1], lty = 2, lwd = 2)
        abline (h = qdeltas[3], lty = 2, lwd = 2)
        
        if (!is.null (truedeltas))
        {   
            abline (h = truedeltas [j,k], lwd = 2, col = "red")
        }
    }
    else
    {
        j <- 1:2
        if (features[1] == 0) j[1] <- 1 else
            j[1] <- which (fithtlr$fsel == features[1]) + 1
        if (features[2] == 0) j[2] <- 1 else        
            j[2] <- which (fithtlr$fsel == features[2]) + 1
        k <- class - 1
        deltas <- fithtlr$mcdeltas[j, k,usedmc, drop = TRUE]
        
        if (symlim)
        {
            lim <- quantile (deltas, probs = drawq)
            xlim <- lim 
            ylim <- lim
        }
        else
        {
            xlim <- quantile (deltas[1,], probs = drawq)
            ylim <- quantile (deltas[2,], probs = drawq)
        }
        
        plot (t(deltas), 
              xlab = sprintf ("feature %d",  features[1]),
              ylab = sprintf ("feature %d",  features[2]),
              xlim = xlim,
              ylim = ylim,
              type = "l", col = "grey", 
              main = sprintf("MC Coef. for Features %d and %d (Class %d)",
                             features[1], features[2], class) )
        points (t(deltas), pch = 20) 
    }
    out <- deltas
}



htlr_sdb <- function (fitbplr,usedmc = NULL, burn = NULL, thin = NULL)
{
    
    mcdims <- dim (fitbplr$mcdeltas)
    p <- mcdims [1] - 1
    K <- mcdims [2]
    C <- K + 1
    no_mcspl <- mcdims[3]
    
    ## index of mc iters used for inference
    if (is.null (burn)) burn <- floor(no_mcspl * 0.2)
    if (is.null (thin)) thin <- 1
    if (is.null (usedmc)) usedmc <- seq (burn + 1, no_mcspl, by = thin)
    
    
    deltas <- htlr_mdcoef (fitbplr, usedmc = usedmc, method = mean)
    sdbs <- comp_sdb (deltas, removeint = T, normalize = F)
        
}

