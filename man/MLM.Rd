% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/betafunctions.R
\name{MLM}
\alias{MLM}
\title{Most Likely Mean of the Standard Beta Probability Density Distribution, Given that the Observation is Considered the Most Likely Observation of the Standard Beta Probability Density Distribution (i.e., the mode).}
\usage{
MLM(alpha, beta, x = NULL, n = NULL)
}
\arguments{
\item{alpha}{Observed alpha value for fitted Standard Beta Probability Density Distribution.}

\item{beta}{Observed beta value for fitted Standard Beta Probability Density Distribution.}

\item{x}{Observed proportion-correct outcome.}

\item{n}{Test-length.}
}
\value{
The expected mean of the Standard Beta probability density distribution, for which the observed mean is the most likely value.
}
\description{
Assuming a prior Standard (two-parameter) Beta Distribution, returns the expected mean of the distribution under the assumption that the observed value is the most likely value of the distribution.
}
\examples{
# Assuming a prior Standard (two-parameter) Beta distribution is fit, which
# yield an alpha parameter of 10 and a beta parameter of 8, calculate the
# true-mean most likely to have produced the observations:
MLM(a = 10, b = 8)
}
