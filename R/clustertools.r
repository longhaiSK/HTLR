
iaxp <- function (xlim)
{
    ilim1 <- floor(min (xlim)) - 1
    ilim2 <- ceiling (max (xlim)) + 1
    c(ilim1, ilim2, ilim2 - ilim1)
}

pathplot <- function (scales, coefs, subset = NULL, ...)
{
    ncoef <- ncol (coefs)
    npath <- length (scales)
    nozeroindex <- seq (1, npath, by = floor(npath/15))
    if (is.null (subset)) subset <- 2:ncoef

    plotargs <- list (...)
    if (is.null (plotargs$ylim)) plotargs$ylim <- range (coefs) 
    if (is.null (plotargs$xlim)) plotargs$xlim <- range (scales) + c(-0.1, 0.1)
    if (is.null (plotargs$xlab)) plotargs$xlab <- "log Scale"
    if (is.null (plotargs$ylab)) plotargs$ylab <- "Coefficient Estimates"
    if (is.null (plotargs$yaxp)) plotargs$yaxp <- iaxp (plotargs$ylim)    
    if (is.null (plotargs$xaxp)) plotargs$xaxp <- iaxp (plotargs$xlim)    
    
    
    do.call (plot, c(plotargs, list (x=scales,y=rep(0,npath),type = "n")) )
    
    for (j in subset)
    {
      points (scales, coefs[,j], type = "l", lty = j, col = j)
      text (scales[npath], coefs[npath, j], j-1, col = j, adj = 1, cex=0.8 )
    }

    abs_coefs <- abs(coefs) [nozeroindex,subset]
    max_abs_coefs <- apply (abs_coefs, 1, max)
    abs_coefs <- abs_coefs/max_abs_coefs        
    no_nozeros <- apply( 1*(abs_coefs > 0.1), 1, sum) 
       
    for (i in 1:length (nozeroindex))
    {
        k <- nozeroindex [i]
        text (scales[k], plotargs$ylim[2], no_nozeros[i], cex=0.8, srt = -90,
              adj = 0)
    }
}


