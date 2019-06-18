rank_k <- function (X,y)
{
    pv_kt <- function (x, g) kruskal.test (x,g)$p.value
    order (apply (X, 2, pv_kt, g = y))
}

## This function ranks all features in terms of F-statistic.
panova <- function(X, y)
{
    ## This function computes the values of F-statistic of all the features.
    n <- length (y)
    nos_g <- as.vector (tapply (rep (1,n), INDEX = y, sum))
    G <- length (nos_g)

    gsum_X <- rowsum (X,y)
    sum_X2 <- colSums (X^2)
    sum_gsum2 <- colSums (gsum_X^2 / nos_g )

    pvars <- (sum_X2 - sum_gsum2) / (n-G) ## pooled variances

    sum_X <- colSums (X)
    gvars <- (sum_gsum2 - sum_X^2 / n) / (G-1) ## variances btw groups
    ## F-statistic
    fstats <- gvars / pvars

    pvalues <- 1 - pf (fstats,  df1 = G - 1, df2 = n - G, ncp = 0)
    pvalues
}

rank_f <- function (X,y)
{
  order (panova (X, y))
}

rank_plain <- function (X, y)
{
  seq(1,ncol (X), by = 1)
}

