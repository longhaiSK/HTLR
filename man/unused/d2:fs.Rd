\name{d2:fs}

\alias{htlr_fss}
\alias{plot_fscore}

\title{Feature Subsets Selection with MCMC Samples}

\description{\code{htlr_hss} divides all of the Markov chain samples into a certain number of subpools, each representing a feature subset. \code{plot_fscore} is a function for plotting feature importance scores or coefficients using histogram lines.}

\usage{
htlr_fss (fithtlr, threshold = 0.1, mfreq = 0.05, sfreq = 0.01, 
          print = TRUE)

plot_fscore (fscores, fsel=1:length (fscores),show_ix = 0.01, ...)
}

\arguments{

\item{fithtlr}{a list containing fitting results by \code{htlr_fitpred}}

\item{fscore}{scores measuring feature importance, such as \code{wsdbs}, \code{msdbs}, or coefficients values}

}

\value{

\item{mcids}{a list of vectors containing the indices of Markov chain iterations associated with different feature subsets.}

\item{ftab}{the fitting results of \code{htlr_fitpred} are a mixture of subpools (posterior modes) of Markov chain samples for different feature subsets. \code{htlr_fss} divides the Markov chain samples into different subpools associated with different feature subsets (posterior modes). The rows of \code{ftab} are different feature subsets, with the 1st column showing the indices of feature subset, the 2nd column showing the relative frequency of Markov chain iterations in the subpool, and the 3rd column showing the coefficients estimated by median with the subpool.}

\item{sdbs}{a list of sdbs (feature importance scores) for all feature subsets.}

\item{wsdbs}{weighted average of sdbs in all feature subsets found. the weights are the relative frequencies of Markov chain samples for different feature subsets.}

\item{msdbs}{the maximum of sdbs in all feature subsets, for discovering features used in feature subsets with small frequencies.}

}


