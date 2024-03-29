% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betafunctions.R
\name{rBetacBinom}
\alias{rBetacBinom}
\title{Random Number Generation under Lord's Beta Compound-Binomial Distribution.}
\usage{
rBetacBinom(x, N, k, l, u, alpha, beta)
}
\arguments{
\item{x}{Number of draws.}

\item{N}{Number of trials.}

\item{k}{Lord's k (see documentation for the \code{Lords.k()} function).}

\item{l}{The lower-bound location parameter of the four-parameter Beta distribution.}

\item{u}{The upper-bound location parameter of the four-parameter Beta distribution.}

\item{alpha}{The first shape-parameter of the four-parameter Beta distribution.}

\item{beta}{The second shape-parameter of the four-parameter Beta distribution.}
}
\description{
Random number generation under Lord's Beta Compound-Binomial distribution, where the Compound-Binomial distribution is Lord's two-term approximation.
}
\note{
For larger values of \code{k}, the distribution can yield negative probabilities which returns an error.
}
\examples{
# To draw a sample of 50 values from a Beta Compound-Binomial distribution
# with number of trials = 100, Lord's k = 1, and probabilities of successful
# trials drawn from a four-parameter Beta distribution with location-
# parameters l = .15 and u = .85, and shape parameters alpha = 6 and
# beta = 4:
rBetacBinom(x = 50, N = 100, k = 1, l = .15, u = .85, alpha = 6, beta = 4)
}
