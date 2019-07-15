htlr <- function (
  x, y, fsel = 1:ncol(X_tr), prior = c("t", "neg", "ghs"), 
  stdzx = TRUE,
  sigmab0 = 2000,  alpha = 1, s = -10, eta = 0,
  iters_h = 1000, iters_rmc = 1000, thin = 100,
  leap_L = 50, leap_L_h = 5, leap_step = 0.3,  hmc_sgmcut = 0.05,
  initial_state = "lasso", alpha.rda = 0.2, verbose = FALSE, .legacy = FALSE)
{
  prior <- match.arg(prior)
  
  stopifnot(length (y_tr) == nrow (X_tr),
            iters_rmc > 0, iters_h > 0, leap_L > 0, leap_L_h > 0, thin > 0)
  
  if (any(table(y) < 2)) 
    stop("Less than 2 cases in some group.")
  
  #attr(fit, "class") <- "htlr"
  
}