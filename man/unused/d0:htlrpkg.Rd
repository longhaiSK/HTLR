\name{d0:htlrpkg}

\alias{htlrpkg}

\title{Examples of Using HTLR}

\description{
This help file gives a step-by-step example of using HTLR package how to fit HTLR model, select feature subsets, look at selected feature subsets, plot feature importance score, and make prediction for test cases with fitting results.
}

\examples{

library (HTLR)
data (gen_grpcor_data)

## creating training and test data sets
tr <- 200
X_tr <- data$X[1:tr,]
y_tr <- data$y[1:tr]
X_ts <- data$X[-(1:tr), ]
y_ts <- data$y[-(1:tr)]

## fit htlr models
fithtlr <- htlr_fitpred (
    y_tr = y_tr, X_tr = X_tr, stdzx = TRUE, ## fsel = c (1,51,101),##  data
    pty = "t", alpha = 1, s = -15,  ## alpha = df and s= log (w) 
    iters_h = 1000, iters_rmc = 1000, thin = 50, ## mcmc iterations
    ## hmc and restricted gibbs sampling settings
    leap_L_h = 5, leap_L = 50, leap_step = 0.5, hmc_sgmcut = 0.05, 
    silence = F)  ## silence = F for printing at mcmc iterations info

## select feature subsets by dividing mcmc samples into subpools
fsshtlr <- htlr_fss (fithtlr)

## draw averaged sdbs in all feature subsets
plot_fscore (fsshtlr$wsdbs)

## draw sdbs (feature importance) for top 1 feature subset
plot_fscore (fsshtlr$sdbs[[1]])

## make prediction on test cases
probs_pred_ts <- htlr_predict (X_ts = X_ts, fithtlr = fithtlr)

## report detailed prediction results with true response
pred_report <- evaluate_pred (probs_pred_ts, y_ts)

## look at MCMC sample if you like
# 1) in all mcmc iterations (maybe see multimodes)
htlr_mccoef (fithtlr, features = c(1,51))
# 2) in only mcmc iterations for feature subset 1
htlr_mccoef (fithtlr, features = c(1,51), usedmc = fsshtlr$mcids[[1]])

}

