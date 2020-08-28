#' Create a Matrix of Markov Chain Samples
#' 
#' The Markov chain samples (without warmup) included in a \code{htlr.fit} object will be coerced to a matrix.
#' 
#' @param x An object of S3 class \code{htlr.fit}.
#' 
#' @param k Coefficients associated with class \code{k} will be drawn. Must be a positive integer in 
#' 1,2,\ldots,C-1 for C-class traning labels (base class 0 can not be chosen). By default the last class
#' is selected. For binary logistic model this argument can be ignored.
#' 
#' @param include.warmup Whether or not to include warmup samples 
#' 
#' @param ... Not used.
#' 
#' @return A matrix with \code{(p + 1)} columns and \code{i} rows, where \code{p} is the number of features 
#' excluding intercept, and \code{i} is the number of iterations after burnin.
#' 
#' @export
#' 
#' @examples 
#' ## No. of features used: 100; No. of iterations after burnin: 15 
#' fit <- htlr(X = colon$X, y = colon$y, fsel = 1:100, iter = 20, warmup = 5)
#' 
#' dim(as.matrix(fit))
#'   
as.matrix.htlr.fit <- function(x, k = NULL, include.warmup = FALSE, ...)
{
  if (is.null(k))
  {
    k <- x$K
    if (k > 1)
      message(
        "'k' was not specified, coefficients associated with the last class will be drawn"
      )
  }
  if (include.warmup)
    mcdeltas <- t(x$mcdeltas[ , k, -1])
  else
    mcdeltas <- t(x$mcdeltas[ , k, get_sample_indice(dim(x$mcdeltas)[3], x$mc.param$iter.rmc)]) 
  colnames(mcdeltas) <- colnames(x$feature$X)
  mcdeltas
}

#' Posterior Summaries  
#' 
#' This function gives a summary of posterior of parameters.
#' 
#' @param object An object of S3 class \code{htlr.fit}.
#' 
#' @param usedmc Indices of Markov chain iterations used for inference. By default all iterations are used.
#' 
#' @param method A function that is used to aggregate the MCMC samples. The default is \code{median}, 
#' other built-in/customized statistical functions such as \code{mean}, \code{sd}, and \code{mad}
#' can also be used. 
#' 
#' @param features A vector of indices (int) or names (char) that specify the parameters we will look at.
#' By default all parameters are selected. 
#' 
#' @param ... Not used.  
#' 
#' @return A point summary of MCMC samples. 
#' 
#' @importFrom magrittr extract2
#' 
#' @export
#' 
#' @examples 
#' set.seed(12345)
#' data("colon")
#' 
#' fit <- htlr(X = colon$X, y = colon$y, fsel = 1:100, iter = 20)
#' summary(fit, features = 1:16)
#'   
summary.htlr.fit <-
  function (object,
            features = 1L:object$p,
            method = median,
            usedmc = get_sample_indice(dim(object$mcdeltas)[3], object$mc.param$iter.rmc),
            ...)
{
  ix.f <- object$feature$fsel[features] %>% na.omit()
  
  mddeltas <- object$mcdeltas[c(1, ix.f + 1), , usedmc, drop = FALSE] %>%
    apply(MARGIN = c(1, 2), FUN = method)
  
  attr(mddeltas, "stats") <- 
    match.call() %>% 
    as.list() %>% 
    extract2("method") %>%
    as.character()
  if (is.null(colnames(mddeltas)))
    attr(mddeltas, "stats") <- "median"
  rownames(mddeltas) <- c("Intercept", names(ix.f))
  colnames(mddeltas) <- paste("class", levels(factor(object$feature$y)))[-1]
  
  return(mddeltas)
}

# @param p.burn.extra Percentage of iterations to be burned out after known warmup iterations
# @param n.burn.extra Number of iterations to be burned out after known warmup iterations,
# will overwrite p.burn.extra if specified
# @param ignore.first First entry should be ignored because it is supposed to be the initial state,
# in rare cases you may want it be FALSE
get_sample_indice <- function(n.sample,
                              iter.rmc,
                              thin = 1,
                              p.burn.extra = 0,
                              n.burn.extra = floor(iter.rmc * p.burn.extra),
                              ignore.first = TRUE)
{
  stopifnot(p.burn.extra >= 0, p.burn.extra < 1)
  n.warmup <- n.sample - iter.rmc
  seq((n.warmup + n.burn.extra + ignore.first), n.sample, thin)
}

# @export
htlr_sdb <- function(fit,
                     burn = NULL,
                     usedmc = NULL,
                     thin = 1)
{
  if (is.null(usedmc))
  {
    if (is.null(burn))
      usedmc <- get_sample_indice(dim(fit$mcdeltas)[3], fit$mc.param$iter.rmc, p.burn.extra = 0.2, thin = thin)
    else
      usedmc <- get_sample_indice(dim(fit$mcdeltas)[3], fit$mc.param$iter.rmc, n.burn.extra = burn, thin = thin)
  }
  deltas <- summary(fit, usedmc = usedmc, method = mean)
  comp_sdb(deltas, removeint = TRUE, normalize = FALSE)
}

#' Get Indices of Non-Zero Coefficients
#' 
#' Get the indices of non-zero coefficients from fitted HTLR model objects.
#' 
#' @param fit An object of S3 class \code{htlr.fit}.
#' 
#' @param cut Threshold on relative SDB to distinguish zero coefficients. 
#' 
#' @return Indices vector of non-zero coefficients in the model.
#' 
#' @export
#' 
#' @examples 
#' set.seed(12345)
#' data("colon")
#' 
#' fit <- htlr(X = colon$X, y = colon$y, fsel = 1:100, iter = 20)
#' nzero_idx(fit)
#' 
nzero_idx <- function(fit, cut = 0.1)
{
  sdb <- htlr_sdb(fit)
  which(sdb > cut * max(sdb))
}

# Plots Markov chain trace or scatterplot
# 
# This function plots Markov chain samples of 1 or 2 features. In plotting for 2 features, 
# gray lines show Markov chain transitions.
# 
# @param fithtlr A list containing fitting results by \code{\link{htlr_fit}}.
# @param features A vector of 1 or 2 numbers representing 1 or 2 features one wishes to look.
# @param class Coefficients associated with \code{class} will be drawn. Must be a positive integer in 
# 1,2,\ldots,C for C-class traning labels.
# @param usedmc Indices of Markov chain iterations used in plottings; one can set it to the 
# indices of Markov chain iterations belonging to the ith feature subset, \code{mcids[[i]]}, 
# found by \code{\link{htlr_fss}}.
# 
# @return A vector of Markov chain sample of 1 coefficient, 
# or an array of Markov chain samples of 2 coefficients.
# 
# @export
# 
htlr_mccoef <-  function (fithtlr,
                          features = 1,
                          class = 2,
                          usedmc = get_sample_indice(dim(fithtlr$mcdeltas)[3], fithtlr$mc.param$iter.rmc),
                          symlim = FALSE,
                          drawq = c(0, 1),
                          truedeltas = NULL
)
{
  mcdims <- dim (fithtlr$mcdeltas)
  p <- mcdims [1] - 1
  K <- mcdims [2]
  features <- features [!is.na (features)]
  
  if (length (features) == 1)
  {
    if (features == 0) j <- 1
    else  j <- which(fithtlr$feature$fsel == features) + 1
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
      j[1] <- which (fithtlr$feature$fsel == features[1]) + 1
    if (features[2] == 0) j[2] <- 1 else        
      j[2] <- which (fithtlr$feature$fsel == features[2]) + 1
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
