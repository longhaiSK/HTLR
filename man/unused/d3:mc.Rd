\name{d3:mc}
\alias{htlr_mccoef}

\title{Plotting Markov chain Trace or Scatterplot}

\description{This function plots Markov chain samples of 1 or 2 features. In plotting for 2 features, gray lines show Markov chain transitions.}

\usage{
htlr_mccoef (
        fithtlr, features = 1, class = 2,  usedmc = "all", 
        symlim = FALSE, drawq = c(0,1), truedeltas = NULL)
}

\arguments{
\item{fithtlr}{a list containing fitting results by \code{htlr_fitpred}}

\item{features}{a vector of 1 or 2 numbers representing 1 or 2 features one wishes to look}

\item{class}{coefficients associated with \code{class} will be drawn.}

\item{usedmc}{indices of Markov chain iterations used in plottings; one can set it to the indices of Markov chain iterations belonging to the ith feature subset, \code{mcids[[i]]}, found by \code{htlr_fss}.}

}

\value{
a vector of Markov chain sample of 1 coefficient, or an array of Markov chain samples of 2 coefficients
}
