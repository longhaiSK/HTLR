
htlr_gendata <- function (n, p, NC = 3, nu = 2, w = 1, X = NULL, betas = NULL )
{
  if (is.null (X))
  {
	X <- matrix (rnorm (n * p), n, p)
  }

  if (is.null (betas))
  {
    sigmasbt <- 1/rgamma (p, nu/2, nu/2) * w
    betas <-  replicate (NC, rnorm (p+1)) * c(0,sqrt (sigmasbt))
  }
  deltas <- betas[,-1, drop = FALSE] - betas[,1]  
  lv <- cbind(1,X) %*% betas
  probs <- exp (lv - as.vector(log_sum_exp(lv)))

  y <- apply (probs, 1, function (prob) sample (1:NC, 1, TRUE, prob))

  list (X = X, y = y, deltas = deltas)

}

gendata_fam <- function (n, muj, A, sd_g = 0, stdx = FALSE)
{
  p <- nrow (muj)
  C <- ncol (muj)
  k <- ncol (A)

#  y <- sample (1:C, n, replace = T)
  y <- rep (1:C, len = n)
  X <- A %*% matrix (rnorm (n*k), k, n) + muj[,y] + rnorm (n * p) * sd_g 
  SGM <- A %*% t (A) + diag (sd_g^2,p)
  
  if (stdx == TRUE)
  {
     mux <- rowMeans (muj)
     sdx <- sqrt ( diag (SGM) + apply (muj,1, var) * (C-1)/C )
     muj <- (muj - mux) / sdx
     SGM <- sweep (SGM, 1, sdx, "/")
     SGM <- sweep (SGM, 2, sdx, "/")
     X <- (X - mux) / sdx
     
  }

  list (X = t(X), y = y, muj = muj, SGM = SGM)

}

