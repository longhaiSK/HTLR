#' Make Prediction on New Data (Advanced)
#' 
#' This function uses MCMC samples from fitted \code{htlrfit} object OR user supplied 
#' regression coefficient to predict the class labels of test cases. 
#' 
#' @param X_ts Matrix of values at which predictions are to be made.
#' 
#' @param fithtlr Fitted HTLR model object.
#' 
#' @param deltas The values of deltas (for example true deltas) used to make prediction; 
#' will override \code{fithtlr} if provided. 
#'  
#' @param burn,thin \code{burn} of Markov chain (super)iterations will be discarded for prediction,
#' and only every \code{thin} are used.
#' 
#' @param usedmc Indices of Markov chain iterations used for inference. 
#' If supplied, \code{burn} and \code{thin} will be ignored.
#' 
#' @param rep.legacy To reproduce (actually incorrect) results in legacy version.
#' See \url{https://github.com/longhaiSK/HTLR/issues/7}. 
#' 
#' @return A matrix of predictive probabilities, with rows for cases, cols for classes.
#' 
#' @export
#' 
#' @keywords internal
#' 
htlr_predict <- function(X_ts, fithtlr = NULL, deltas = NULL, burn = NULL, thin = 1, usedmc = NULL, rep.legacy = TRUE)
{
  if (is.vector(X_ts))
    X_ts <- matrix(X_ts, 1)
  else if (!is.matrix(X_ts))
    X_ts <- as.matrix(X_ts)
  no_ts <- nrow(X_ts)
  
  if (is.null(deltas) & !is.null(fithtlr))
  {
    mcdims <- dim(fithtlr$mcdeltas)
    p <- mcdims[1] - 1
    K <- mcdims[2]
    no_mcspl <- mcdims[3]
    
    ## index of mc iters used for inference
    if (is.null(usedmc))
    {
      if (is.null(burn))
        usedmc <- get_sample_indice(no_mcspl, fithtlr$mc.param$iter.rmc, p.burn.extra = 0.1, thin = thin, ignore.first = !rep.legacy)
      else
        usedmc <- get_sample_indice(no_mcspl, fithtlr$mc.param$iter.rmc, n.burn.extra = burn, thin = thin, ignore.first = !rep.legacy)
    }
    
    no_used <- length(usedmc)
    
    ## read deltas for prediction
    longdeltas <- matrix(fithtlr$mcdeltas[, , usedmc], nrow = p + 1)

    ## selecting features and standardizing
    fsel <- fithtlr$feature$fsel
    X_ts <- X_ts[, fsel, drop = FALSE]
    nuj <- fithtlr$feature$nuj
    sdj <- fithtlr$feature$sdj
    X_ts <- sweep(X_ts, 2, nuj, "-")
    X_ts <- sweep(X_ts, 2, sdj, "/")
  }
  else
  {
    if (is.vector(deltas) | is.matrix(deltas)) 
    {
      deltas <- matrix(deltas, nrow = ncol(X_ts) + 1)
      p <- nrow(deltas) - 1
      K <- 1
      longdeltas <- deltas
      no_used <- 1
    }
  }
  
  ## add intercept to all cases
  X_addint_ts <- cbind(1, X_ts)
  
  longlv <- X_addint_ts %*% longdeltas
  arraylv <- array(longlv, dim = c(no_ts, K, no_used))
  logsumlv <- apply(arraylv, 3, comp_lsl)
  array_normlv <- sweep(arraylv, c(1, 3), logsumlv)
  array_predprobs <- exp(array_normlv)
  probs_pred <- apply(array_predprobs, c(1, 2), mean)
  
  predprobs_c1 <- pmax(0, 1 - apply(probs_pred, 1, sum))
  probs_pred <- cbind(predprobs_c1, probs_pred)
  colnames(probs_pred) <- paste("class", levels(factor(fithtlr$feature$y)))
  
  return(probs_pred)
}

#' Make Prediction on New Data
#' 
#' Similar to other predict methods, this function returns predictions from a fitted \code{htlrfit} object.
#' 
#' @param object A fitted model object with S3 class \code{htlrfit}.
#' 
#' @param newx A Matrix of values at which predictions are to be made.
#' 
#' @param type Type of prediction required. Type "response" gives the fitted probabilities.
#' Type "class" produces the class label corresponding to the maximum probability.
#' 
#' @param ... Advanced options to specify the Markov chain iterations used for inference. 
#' See \code{\link{htlr_predict}}.   
#' 
#' @return The object returned depends on type.
#' 
#' @export
#' 
predict.htlr.fit <- function(object, newx, type = c("response", "class"), ...)
{
  if (!exists("burn")) burn <- NULL
  if (!exists("usedmc")) usedmc <- NULL
  if (!exists("thin")) thin <- 1
  
  pred.prob <- htlr_predict(X_ts = newx, fithtlr = object, burn = burn, thin = thin, usedmc = usedmc, rep.legacy = FALSE)
  
  type <- match.arg(type)
  if (type == "response")
    return(pred.prob)
  if (type == "class") {
    classes <- object$feature$y %>%
      factor() %>%
      levels() %>%
      as.numeric()
    newy <- classes[apply(pred.prob, 1, which.max)] %>% as.matrix()
    colnames(newy) <- "y.pred"
    return(newy)
  }
}
