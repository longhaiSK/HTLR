rank_k <- function (X,y)
{
	pv_kt <- function (x, g) kruskal.test (x,g)$p.value
	order (apply (X, 2, pv_kt, g = y))
}

## This function ranks all features in terms of F-statistic.
rank_F <- function(X, y)
{   
    ## This function computes the values of F-statistic of all the features.
    comp_fstats <- function (X,y)
    {
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
	gvars / pvars
    }

    fstats <- comp_fstats (X, y)
    vars <- order (fstats, decreasing = TRUE)
    list (vars=vars, fstats = fstats [vars])
}

rank_plain <- function (X, y)
{	
  seq(1,ncol (X), by = 1)
}
