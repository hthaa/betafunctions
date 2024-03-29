% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classification.R
\name{caStats}
\alias{caStats}
\title{Classification Accuracy Statistics.}
\usage{
caStats(tp, tn, fp, fn)
}
\arguments{
\item{tp}{The frequency or rate of true-positive classifications.}

\item{tn}{The frequency or rate of true-negative classifications.}

\item{fp}{The frequency or rate of false-positive classifications.}

\item{fn}{The frequency or rate of false-negative classifications.}
}
\value{
A list of diagnostic performance statistics based on true/false positive/negative statistics. Specifically, the sensitivity, specificity, positive predictive value (PPV), negative predictive value (NPV), Youden's J. (Youden.J), and Accuracy.
}
\description{
Provides a set of statistics often used for conveying information regarding the certainty of classifications based on tests.
}
\examples{
# Generate some fictional data. Say, 1000 individuals take a test with a
# maximum score of 100 and a minimum score of 0.
set.seed(1234)
testdata <- rbinom(1000, 100, rBeta.4P(1000, 0.25, 0.75, 5, 3))
hist(testdata, xlim = c(0, 100))

# Suppose the cutoff value for attaining a pass is 50 items correct, and
# that the reliability of this test was estimated to 0.7. First, compute the
# estimated confusion matrix using LL.CA():
cmat <- LL.CA(x = testdata, reliability = 0.7, cut = 50, min = 0,
max = 100)$confusionmatrix

# To estimate and retrieve diagnostic performance statistics using caStats(),
# feed it the appropriate entries of the confusion matrix.
caStats(tp = cmat["True", "Positive"], tn = cmat["True", "Negative"],
fp = cmat["False", "Positive"], fn = cmat["False", "Negative"])
}
