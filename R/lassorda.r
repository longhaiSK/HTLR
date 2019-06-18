if (T) ## not needed in building package
{

library ("glmnet")
library ("rda")


lasso_fitpred <- function (X_tr, y_tr, X_ts = NULL, 
                          rank_fn = rank_plain, k = ncol (X_tr) )
{
    ## read information about data
    n <- nrow (X_tr) ## numbers of obs
    p <- ncol (X_tr)
    ## find number of observations in each group
    nos_g <- as.vector (tapply (rep(1,n),INDEX = y_tr, sum))
    G <- length (nos_g)
    if (any(nos_g < 2)) stop ("Less than 2 cases in some group")

    ## feature selection
    if (k < p)
    {
      fsel <- rank_fn (X_tr, y_tr) [1:k]
      X_tr <- X_tr [, fsel, drop = FALSE]
      X_ts <- X_ts[,fsel, drop = FALSE]
    }

    ## choosing the best lambda
    cvfit <- cv.glmnet (x = X_tr, y = y_tr,nlambda = 500, 
                    family = "multinomial")
    lambda <- cvfit$lambda[which.min(cvfit$cvm)]
    cat ("The best lambda chosen by CV:",lambda,"\n")

    ## fit lasso with the best lambda
    lassofit <- glmnet (x = X_tr, y= y_tr, nlambda = 500,  
                        family = "multinomial")

    betas <- coef (lassofit, s = lambda)

    mbetas <- matrix (0, p + 1, G)
    for (g in 1:G)
    {
        mbetas[,g] <- as.numeric (betas[[g]])
    }
    deltas <- mbetas[, -1, drop = FALSE] - mbetas[,1]

    ## predicting for new cases
    if (is.null (X_ts))
    {
        return (deltas)
    }
    else
    {
        probs_pred <- matrix(
                  predict (lassofit, newx = X_ts,
                  s = lambda, type = "response")[,,1], 
                  nrow = nrow (X_ts))
        return (
            list (probs_pred = probs_pred,
                  values_pred = apply (probs_pred, 1, which.max),
                  deltas = deltas)
        )
    }
}


lassocv_fsel_trpr <- function (
  y_tr, X_tr, X_ts, nos_fsel = ncol (X_tr), rankf = rank_k)
{
  no_ts <- nrow (X_ts)
  rankedf <- rankf (X = X_tr, y = y_tr)

  nfsel <- length (nos_fsel)
  NC <- length (unique (y_tr))
  array_probs_pred <- array (0, dim = c(no_ts, NC, nfsel))

  for (i in 1:nfsel)
  {
    fsel <- rankedf [1:nos_fsel[i]]
    array_probs_pred[,,i] <- trpr_lassocv (
      y_tr = y_tr, X_tr = X_tr[, fsel, drop = FALSE],
      X_ts = X_ts[, fsel, drop = FALSE])$probs_pred
  }

  list (nos_fsel = nos_fsel, array_probs_pred = array_probs_pred)
}




convert_ix <- function (n, nrow, ncol)
{
  m <- floor ((n - 1) %/% nrow)
  r <- n - 1 - m * nrow
  c (r + 1, m + 1)
}

trpr_rdacv <- function (X_tr, y_tr, X_ts, nos_fsel = ncol (X_tr),
              rankf = rank_plain)
{
  topgenes <- rankf (X_tr, y_tr)
  x_tr <- t (X_tr)
  x_ts <- t (X_ts)

  n <- ncol (x_ts)
  C <- length (unique (y_tr))

  K <- length (nos_fsel)

  array_probs_pred <- array (0, dim = c(n, C, K))

  for (k in 1:K )
  {
    nofsel <- nos_fsel [k]
    fsel <- topgenes [1:nofsel]
    x_tr_sel <- x_tr [fsel,, drop = FALSE]
    x_ts_sel <- x_ts [fsel,, drop = FALSE]

    fit <- rda (x = x_tr_sel, y = y_tr)
    fitcv <- rda.cv (fit, x = x_tr_sel, y = y_tr)
    no_delta <- length (fitcv$delta)
    no_alpha <- length (fitcv$alpha)
    ixmin <- which.min(fitcv$cv.err)
    ixad <- convert_ix (ixmin, no_alpha, no_delta)
    opt_alpha <- fitcv$alpha[ixad[1]]
    opt_delta <- fitcv$delta[ixad[2]]

    array_probs_pred[,,k] <- rda (x = x_tr_sel, y = y_tr, xnew = x_ts_sel,
    alpha = opt_alpha, delta = opt_delta )$posterior[1,1,,]

    array_probs_pred[,,k] <- exp (array_probs_pred[,,k])
    sumprobs <- apply (array_probs_pred[,,k], 1, sum)
    array_probs_pred[,,k] <- array_probs_pred[,,k] / sumprobs
  }

  array_probs_pred

}


rdacv_fsel_trpr <- function (X_tr, y_tr, X_ts, nos_fsel = ncol (X_tr),
    rankf = rank_plain)
{
  list (array_probs_pred = trpr_rdacv (X_tr = X_tr, y_tr = y_tr, X_ts = X_ts, nos_fsel = nos_fsel, rankf = rankf), nos_fsel = nos_fsel )
}

}
