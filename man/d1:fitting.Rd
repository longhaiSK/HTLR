\name{d1:fitting}

\alias{htlr_fit}
\alias{htlr_predict}

\title{Fitting HTLR Models and Making Predictions}

\description{
\code{htlr_fit} trains linear logistic regression models with HMC in restricted Gibbs sampling. This function also makes predictions for test cases if \code{X_ts} are provided.

\code{htlr_predict} uses MCMC samples returned by \code{htlr_fit} to predict the class labels of test cases. Prediction results are a matrix of predictive probabilities and a vector of predicted class labels }

\usage{

htlr_fit (
    y_tr, X_tr, X_ts = NULL, fsel = 1:ncol(X_tr), stdzx = TRUE,
    sigmab0 = 2000, ptype = "t", alpha = 1, s = -15, eta = 0,  
    iters_h = 2000, iters_rmc = 2000, thin = 10,  
    leap_L = 50, leap_L_h = 5, leap_step = 0.5,  hmc_sgmcut = 0.05,
    looklf = 0,  initial_state = NULL, silence = TRUE,
    predburn = NULL, predthin = 1)

htlr_predict (X_ts, fithtlr = NULL, deltas = NULL,
              burn = NULL, thin = NULL, usedmc = NULL)

}

\arguments{

\item{X_tr, X_ts}{matrices containing data; rows should be for the cases, and columns for different features; \code{X_tr} are training data, \code{X_ts} are test data or future data for which prediction are needed.}

\item{fsel}{subsets of features selected before fitting, such as by univariate screening.}

\item{stdzx}{if it is set to TRUE, the original features values are standardized to have mean 0 and sd 1 for each gene; by default, it is TRUE.}

\item{y_tr,y}{a vector of class labels in training or test data set. Must be coded as positive integers 1,2,...,C for C classes.}

\item{iters_h, iters_rmc, thin}{\code{iters_h} and \code{iters_rmc}  of super Markov chain transitions, each with \code{thin} Markov chain iterations, are run for burning (aka head/initial) and sampling phases; only the last state of each super transition in the sampling phase is saved.}

\item{alpha, s, sigmab0, ptype, eta}{Prior settings for coefficients. \code{alpha} is the degree freedom, and s is log (w), equal to twice of log scale of priors for coefficients.} 

\item{silence}{logical. Setting it to TRUE for tracking MCMC sampling iterations.}

\item{fithtlr}{a list containing fitting results by \code{htlr_fit}}

\item{predburn,predthin}{\code{predburn} of Markov chain (super)iterations will be discarded for prediction, and only every \code{predthin}th are used; by default, 20\% of (super)iterations are burned, and \code{thin}=1.}

\item{burn, thin}{the meanings are the same as \code{predburn} and \code{predthin}}

}

\value{

\code{htlr_fit} returns a list of fitting results. \code{htlr_predict} returns a matrix of predictive probabilities, with rows for cases, cols for classes. }

\keyword{classif}

