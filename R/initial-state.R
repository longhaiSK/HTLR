as.htlr.init <- function(x) UseMethod("as.htlr.init")

as.htlr.init.cv.glmnet <- function(x) 
{
  coefs <- coef(x, s = "lambda.min") %>% Reduce(f = cbind) %>% as.matrix()
  deltas <- coefs[, -1, drop = FALSE] - coefs[, 1]
  return (deltas)
}




