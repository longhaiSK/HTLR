predict_bplr <- function (X_test, K, mcchain, start, thin)
{
    p <- ncol (X_test)

    X_test <- cbind (1, X_test)
    
    indice_pred <- 
       start + seq (0, floor((ncol (mcchain) - start)/thin)) * thin
    
    comp_prob <- function(deltas)
    {  
       deltas <- matrix (deltas, p + 1, K)
       lv <- X_test %*% deltas
       as.vector (exp (lv - comp_lsl(lv)) )
    }

    predprob <- matrix (
                   apply (apply (mcchain[1:((p+1)*K), indice_pred], 2, comp_prob),
                   1, mean), nrow (X_test), K)
    predprob <- cbind (1 - apply (predprob,1, sum), predprob)
    
    list (predprob = predprob, predvalue = apply (predprob, 1, which.max) - 1)
}

# cv <- function (y, X, no_fold, )

