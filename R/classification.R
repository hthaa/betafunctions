#' Livingston and Lewis' "Effective Test Length".
#'
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param min The lower-bound (minimum possible value) of the observed-score distribution. Default is 0 (assuming observed scores represent proportions).
#' @param max The upper-bound (maximum possible value) of the observed-score distribution. Default is 1 (assuming observed scores represent proportions).
#' @param reliability The reliability of the observed scores (proportion of observed-score distribution variance shared with true-score distribution).
#' @return An estimate of the effective length of a test, given the stability of the observations it produces.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the effective test length using ETL():
#' ETL(mean = mean(testdata), variance = var(testdata), min = 0, max = 100,
#' reliability = .7)
#' @export
ETL <- function(mean, variance, min = 0, max = 1, reliability) {
  ((mean - min) * (max - mean) - (reliability * variance)) / (variance * (1 - reliability))
}

#' An Implementation of the Livingston and Lewis (1995) Approach to Estimate Classification Consistency and Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores for which a Beta true-score distribution is to be estimated, or a list of pre-defined true-score distribution parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location parameters, \code{alpha} and \code{beta} for the shape parameters, and \code{etl} for the effective test length (see documentation for the \code{ETL} function).
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0 (assuming \code{x} represent proportions).
#' @param max The maximum value possible to attain on the test. Default is 1 (assuming \code{x} represent proportions).
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param error.model The probability distribution to be used for producing the sampling distributions at different points of the true-score scale. Options are \code{binomial} and \code{beta}. The binomial distribution is discrete, and is the distribution used originally by Livingston and Lewis. Use of the binomial distribution involves a rounding of the effective test length to the nearest integer value. The Beta distribution is continuous, and does not involve rounding of the effective test length.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves (see documentation for the \code{LL.ROC()} function).
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param failsafe Logical value indicating whether to engage the automatic failsafe defaulting to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{FALSE} (i.e., the function will not engage failsafe, and will likely produce an error if impermissible parameter estimates were produced.
#' @param l If \code{true.model = "2P"} or \code{failsafe = TRUE}, the lower-bound location parameter to be used in the two-parameter fitting procedure. Default is 0 (i.e., the lower-bound of the Standard Beta distribution).
#' @param u If \code{true.model = "2P"} or \code{failsafe = TRUE}, the upper-bound location parameter to be used in the two-parameter fitting procedure. Default is 1 (i.e., the upper-bound of the Standard Beta distribution).
#' @param override Inert artifact from betafunctions version 1.3.1 (replaced by the \code{failsafe} argument). Will be removed completely in a later update.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note It should be noted that this implementation differs from the original articulation of Livingston and Lewis (1995) in some respects. First, the procedure includes a number of diagnostic performance (accuracy) indices which the original procedure enables but that were not included. Second, the possibility of employing a two-parameter Beta error distribution in place of the binomial error distribution is not part of the original procedure. Third, the way consistency is calculated differs substantially from the original articulation of the procedure, which made use of a split-half approach. Rather, this implementation uses the approach to calculating classification consistency outlined by Hanson (1991).
#' @note A shiny application providing a GUI for this method is available at https://hthaa.shinyapps.io/shinybeta/ .
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, 0.25, 0.75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the estimated parameters, confusion matrix, consistency and
#' # accuracy statistics using LL.CA():
#' LL.CA(x = testdata, reliability = .7, cut = 50, min = 0, max = 100)
#'
#' # Suppose the true-score parameter estimation procedure arrived at
#' # impermissible parameter estimates (i.e., l < 0, u > 1, alpha < 0, or
#' # beta < 0). For example:
#' set.seed(9)
#' testdata <- rbinom(100, 25, rBeta.4P(100, 0.25, 1, 5, 3))
#' Beta.tp.fit(testdata, 0, 25, 25, failsafe = TRUE)
#'
#' # Suppose further that you have good grounds for assuming that the lower-
#' # bound parameter is equal to 0.25 (e.g., the test consists of multiple-
#' # choice questions with four response options, leading to a 25% probability
#' # of guessing the correct answer per question), and good reason to believe
#' # that the upper-bound parameter is equal to 1 (i.e., there is no reason to
#' # believe that there are no members of the population who will attain a
#' # perfect score across all possible test-forms.) To set these lower and
#' # upper bounds for the fitting procedure in the LL.CA() function, set
#' # the argument true.model = "2p", and specify the location parameters
#' # l = 0.25 and u = 1:
#' LL.CA(testdata, 0.6287713, 12, 0, 25, true.model = "2p", l = 0.25, u = 1)
#'
#' # Alternatively to supplying scores to which a true-score distribution is
#' # to be fit, a list with true-score distribution parameter values can be
#' # supplied manually along with the effective test length (see documentation
#' # for the ETL() function), foregoing the need for actual data. The list
#' # entries must be named. "l" is the lower-bound and "u" the upper-bound
#' # location parameters of the true-score distribution, "alpha" and "beta" for
#' # the shape parameters, and "etl" for the effective test-length..
#' trueparams <- list("l" = 0.25, "u" = 0.75, "alpha" = 5, "beta" = 3, "etl" = 50)
#' LL.CA(x = trueparams, cut = 50, min = 0, max = 100)
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @export
LL.CA <- function(x = NULL, reliability, cut, min = 0, max = 1, true.model = "4P", error.model = "binomial", truecut = NULL, output = c("accuracy", "consistency"), failsafe = FALSE, l = 0, u = 1, override = NULL) {
  if (!is.null(override)) {
    warning("The override argument is rendered inert as of betafunctions v. 1.4.0, replaced by the failsafe argument.
            The argument will be removed entirely in a future update.")
  }
  out <- base::list()
  if (class(x) != "list") {
    if ((base::min(x) < min) | (base::max(x) > max)) {
      warning(paste("Observed values not within the specified [", min, ", ", max, "] bounds (observed min = ",
                    min(x), ", observed max = ", max(x), ").", sep = ""))
    }
    N <- ETL(base::mean(x), stats::var(x), min = min, max = max, reliability = reliability)
    params <- Beta.tp.fit(x, min = min, max = max, etl = N, true.model = true.model, failsafe = failsafe, l = l, u = u)
    if (params$l < 0 | params$u > 1 | params$alpha < 0 | params$beta < 0) {
      warning(paste("Parameter out of bounds: l = ", round(params$l, 4), ", u = ", round(params$u, 4), ", alpha = ", round(params$alpha, 4), ", beta = ", round(params$beta, 4),
                    ". Consider constraining the fitting procedure further (e.g., set the location-parameters).", sep = ""))
    }
    x <- (x - min) / (max - min)
  } else {
    params <- x
    N <- params$etl
  }
  if (base::is.null(truecut)) {
    truecut <- cut
  }
  cut <- (cut - min) / (max - min)
  truecut <- (truecut - min) / (max - min)
  if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
    N <- base::round(N)
  }
  out[["effectivetestlength"]] <- N
  out[["parameters"]] <- params
  if (any(output == "accuracy") | any(output == "Accuracy") | any(output == "ca") |
      any(output == "CA") | any(output == "a") | any(output == "A")) {
    if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
      p.tp <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut)
        }, lower = truecut, upper = 1)$value
      p.fp <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut)
        }, lower = 0, upper = truecut)$value
      p.ff <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE)
        }, lower = truecut, upper = 1)$value
      p.tf <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE)
        }, lower = 0, upper = truecut)$value
    }
    if (error.model == "beta" | error.model == "Beta") {
      p.tp <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut)
        }, lower = truecut, upper = 1)$value
      p.fp <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut)
        }, lower = 0, upper = truecut)$value
      p.ff <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE)
        }, lower = truecut, upper = 1)$value
      p.tf <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE)
        }, lower = 0, upper = truecut)$value
    }
    camat <- confmat(p.tf, p.tp, p.ff, p.fp, "prop")
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(output == "consistency") | any(output == "Consistency" ) | any(output == "cc") |
      any(output == "CC") | any(output == "c") | any(output == "C")) {
    if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
      p.ii <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE) *
          stats::pbinom(floor(cut * N), N, x, lower.tail = TRUE )
        }, lower = 0, upper = 1)$value
      p.ij <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE) *
          (1 - stats::pbinom(floor(cut * N), N, x, lower.tail = TRUE))
        }, lower = 0, upper = 1)$value
      p.jj <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = FALSE) *
          (1 - stats::pbinom(floor(cut * N), N, x, lower.tail = TRUE))
        }, lower = 0, upper = 1)$value
      p.ji <- stats::integrate(function(x) {
        dBeta.pBinom(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = FALSE) *
          stats::pbinom(floor(cut * N), N, x, lower.tail = TRUE )
        }, lower = 0, upper = 1)$value
    }
    if (error.model == "beta" | error.model == "Beta") {
      p.ii <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE) *
          stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = TRUE)
        }, lower = 0, upper = 1)$value
      p.ij <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = TRUE) *
          (1 - stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = TRUE))
        }, lower = 0, upper = 1)$value
      p.jj <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = FALSE) *
          (1 - stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = TRUE))
        }, lower = 0, upper = 1)$value
      p.ji <- stats::integrate(function(x) {
        dBeta.pBeta(x, params$l, params$u, params$alpha, params$beta, N, cut, lower.tail = FALSE) *
          stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = TRUE)
        }, lower = 0, upper = 1)$value
    }
    ccmat <- base::matrix(nrow = 2, ncol = 2, dimnames = list(c("i", "j"), c("i", "j")))
    ccmat["i", "i"] <- p.ii
    ccmat["i", "j"] <- p.ij
    ccmat["j", "i"] <- p.ji
    ccmat["j", "j"] <- p.jj
    out[["consistencymatrix"]] <- ccmat / sum(ccmat)
    out[["classification.consistency"]] <- ccStats(ccmat["i", "i"], ccmat["i", "j"], ccmat["j", "i"], ccmat["j", "j"])
  }
  return(out)
}

#' Confusion matrix
#'
#' @description Organizes supplied values of true and false positives and negatives into a confusion matrix.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @param output Whether the returned output reflects frequencies or proportions. Defaults to returning frequencies.
#' @return A confusion matrix organizing the input values of true and false positive and negatives.
#' @examples
#' # Generate some true and observed conditions.
#' set.seed(1234)
#' true.ability <- rbeta(50, 4, 4)
#' true.category <- ifelse(true.ability < 0.5, 0, 1)
#' observed.score <- rbinom(50, 50, true.ability)
#' observed.category <- ifelse(observed.score < 25, 0, 1)
#' # Calculate the frequencies of true and false positives and negatives based on the true and
#' # observed conditions.
#' TP <- sum(ifelse(observed.category == 0 & true.category == 0, 1, 0))
#' FP <- sum(ifelse(observed.category == 0 & true.category != 0, 1, 0))
#' TN <- sum(ifelse(observed.category == 1 & true.category == 1, 1, 0))
#' FN <- sum(ifelse(observed.category == 1 & true.category != 1, 1, 0))
#' # Organize the above values in a confusion matrix using the confmat function:
#' confmat(tp = TP, fp = FP, tn = TN, fn = FN)
#' @export
confmat <- function(tp, tn, fp, fn, output = "freq") {
  mat <- matrix(nrow = 3, ncol = 3)
  rownames(mat) <- c("True", "False", "Total")
  colnames(mat) <- c("Positive", "Negative", "Total")
  tot <- sum(tp, tn, fp, fn)
  mat[1, 1] <- tp
  mat[1, 2] <- tn
  mat[2, 1] <- fp
  mat[2, 2] <- fn
  mat[, 3] <- rowSums(mat[, -3])
  mat[3, ] <- colSums(mat[-3, ])
  if (output != "freq") {
    mat <- mat / sum(tp, tn, fp, fn)
  }
  mat
}

#' Classification Accuracy Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the certainty of classifications based on tests.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @return A list of diagnostic performance statistics based on true/false positive/negative statistics. Specifically, the sensitivity, specificity, positive likelihood ratio (LR.pos), negative likelihood ratio (LR.neg), positive predictive value (PPV), negative predictive value (NPV), Youden's J. (Youden.J), and Accuracy.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. First, compute the
#' # estimated confusion matrix using LL.CA():
#' cmat <- LL.CA(x = testdata, reliability = .7, cut = 50, min = 0,
#' max = 100)$confusionmatrix
#'
#' # To estimate and retrieve diagnostic performance statistics using caStats(),
#' # feed it the appropriate entries of the confusion matrix.
#' caStats(tp = cmat["True", "Positive"], tn = cmat["True", "Negative"],
#' fp = cmat["False", "Positive"], fn = cmat["False", "Negative"])
#' @references Glas et al. (2003). The Diagnostic Odds Ratio: A Single Indicator of Test Performance, Journal of Clinical Epidemiology, 1129-1135, 56(11). doi: 10.1016/S0895-4356(03)00177-X
#' @export
caStats <- function(tp, tn, fp, fn) {
  sensitivity <-  tp / (tp + fn)
  specificity <-  tn / (tn + fp)
  plr <-          sensitivity / (1 - specificity)
  nlr <-          (1 - sensitivity) / specificity
  ppv <-          tp / (tp + fp)
  npv <-          tn / (tn + fn)
  accuracy <-     (tp + tn) / (tp + tn + fp + fn)
  J <-            (sensitivity + specificity) - 1
  base::list("Sensitivity" = sensitivity, "Specificity" = specificity,
             "LR.pos" = plr, "LR.neg" = nlr,
             "PPV" = ppv, "NPV" = npv,
             "Youden.J" = J, "Accuracy" = accuracy)
}

#' Classification Consistency Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the consistency of classifications based on tests.
#' @param ii The frequency or rate of consistent classifications into category "i".
#' @param ij The frequency or rate of inconsistent classifications into categories "i" and "j".
#' @param ji The frequency or rate of inconsistent classifications into categories "j" and "i".
#' @param jj The frequency or rate of consistent classifications into category "j".
#' @return A list of classification consistency statistics. Specifically, the coefficient of consistent classification (p), the coefficient of consistent classification by chance (p_c), and Cohen's Kappa coefficient.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. First, compute the
#' # estimated consistency matrix using LL.CA():
#' cmat <- LL.CA(x = testdata, reliability = .7, cut = 50, min = 0,
#' max = 100)$consistencymatrix
#'
#' # To estimate and retrieve consistency statistics using ccStats(),
#' # feed it the appropriate entries of the consistency matrix.
#' ccStats(ii = cmat["i", "i"], ij = cmat["i", "j"],
#' ji = cmat["j", "i"], jj = cmat["j", "j"])
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @export
ccStats <- function(ii, ij, ji, jj) {
  p <- (ii + jj) / (ii + ij + ji + jj)
  p_c <- (ii + ij) * (ii + ji) + (ij + jj) * (ji + jj)
  Kappa <- (p - p_c) / (1 - p_c)
  base::list("p" = p, "p_c" = p_c, "Kappa" = Kappa)
}

#' ROC curves for the Livingston and Lewis approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed proportion-score scale.
#' @param x A vector of observed results.
#' @param min The minimum possible value to attain on the observed-score scale. Default is 0 (assuming \code{x} represent proportions).
#' @param max The maximum possible value to attain on the observed-score scale. Default is 1 (assuming \code{x} represent proportions).
#' @param reliability The reliability coefficient of the test.
#' @param truecut The true point along the x-scale that marks the categorization-threshold.
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param error.model The probability distribution to be used for producing the sampling distributions at different points of the true-score scale. Options are \code{binomial} and \code{beta}. The binomial distribution is discrete, and is the distribution used originally by Livingston and Lewis. Use of the binomial distribution involves a rounding of the effective test length to the nearest integer value. The Beta distribution is continuous, and does not involve rounding of the effective test length.
#' @param failsafe If true-model == "4P": Whether to engage a fail-safe reverting to a two-parameter true-score distribution solution should the four-parameter fitting procedure produce impermissible results.
#' @param l If true-model == "2P" or failsafe == TRUE: The lower-bound location parameter of the two-parameter true-score distribution solution.
#' @param u If true-model == "2P" or failsafe == TRUE: The upper-bound location parameter of the two-parameter true-score distribution solution.
#' @param AUC Calculate and include the area under the curve? Default is \code{FALSE}.
#' @param maxJ Mark the point along the curve where Youden's J statistic is maximized? Default is \code{FALSE}.
#' @param raw.out Give raw coordinates as output rather than plot? Default is \code{FALSE}.
#' @param grainsize Specify the number of cutoff-points for which the ROC curve is to be calculated. The greater this number the greater the accuracy. Default is 100 points.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is \code{TRUE}.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To produce a plot
#' # with an ROC curve using LL.ROC(), along with the AUC statistics and the
#' # points at which Youden's J. is maximized:
#' LL.ROC(x = testdata, reliability = .7, truecut = 50, min = 0, max = 100,
#' AUC = TRUE, maxJ = TRUE)
#' @export
LL.ROC <- function(x = NULL, reliability, min = 0, max = 1, truecut, true.model = "4P", error.model = "Binomial", failsafe = FALSE, l = 0, u = 1, AUC = FALSE, maxJ = FALSE, raw.out = FALSE, grainsize = 100) {
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  for (i in 1:(grainsize + 1)) {
    if (i == 1) {
      cuts <- seq(min, max, (max - min) / grainsize)
      outputmatrix <- matrix(nrow = grainsize + 1, ncol = 4)
      outputmatrix[, 4] <- cuts
    }
    axval <- LL.CA(x = x, min = min, max = max, reliability = reliability, cut = cuts[i],
                   truecut = truecut, true.model = true.model, error.model = error.model,
                   output = "a", l = l, u = u)$classification.accuracy
    outputmatrix[i, 1] <- 1 - axval$Specificity
    outputmatrix[i, 2] <- axval$Sensitivity
    outputmatrix[i, 3] <- axval$Youden.J
    colnames(outputmatrix) <- c("FPR", "TPR", "Youden.J", "Cutoff")
    outputmatrix[which(is.na(outputmatrix[, 1])), 1] <- 0
    outputmatrix[which(is.na(outputmatrix[, 2])), 2] <- 1
  }
  if (raw.out) {
    return(outputmatrix)
  }
  graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
  graphics::abline(h = seq(0, 1, .1), v = seq(0, 1, .1), col = "lightgrey", lty = "dotted")
  graphics::par(new = TRUE)
  graphics::plot(outputmatrix[, 1], outputmatrix[, 2], type = "l",
       xlab = "False-Positive Rate (1 - Specificity)",
       ylab = "True-Positive Rate (Sensitivity)",
       main = paste("ROC curve for true-cut equal to", truecut), lwd = 2,
       xlim = c(0, 1), ylim = c(0, 1))
  if (AUC) {
    graphics::legend("bottomright", bty = "n", cex = 1.5,
                     legend = paste("AUC =", round(AUC(outputmatrix[, 1], outputmatrix[, 2]), 3)))
  }
  if (maxJ) {
    graphics::points(outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 1],
           outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 2], cex = 1.5, pch = 19)
    graphics::text(outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 1] + .025,
         outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 2] - .025,
         labels = paste("Maximum Youden's J. at cutoff = ",
                        round(outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3]))[1], 4], 3),
                        "\n(Max. Youden's J. = ", round(max(outputmatrix[, 3]), 3), ").", sep = ""),
         adj = c(0, 1))
  }
}

#' Area Under the ROC Curve.
#'
#' @description Given a vector of false-positive rates and a vector of true-positive rates, calculate the area under the Receiver Operator Characteristic (ROC) curve.
#' @param FPR Vector of False-Positive Rates.
#' @param TPR Vector of True-Positive Rates.
#' @return A value representing the area under the ROC curve.
#' @note Script originally retrieved and modified from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To calculate the
#' # necessary (x, y) coordinates to compute the area under the curve statistic
#' # one can use the LL.ROC() function with the argument
#' # raw.out = TRUE.
#' coords <- LL.ROC(x = testdata, reliability = .7, truecut = 50, min = 0,
#' max = 100, raw.out = TRUE)
#'
#' # To calculate and retrieve the Area Under the Curve (AUC) with the AUC()
#' # function, feed it the raw coordinates calculated above.
#' AUC(coords[, "FPR"], coords[, "TPR"])
#' @export
AUC <- function(FPR, TPR) {
  dFPR <- base::c(diff(FPR), 0)
  dTPR <- base::c(diff(TPR), 0)
  base::sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

#' Calculate Cronbach's Alpha from supplied variables.
#'
#' @description Calculates Cronbach's Alpha, a very commonly used index for assessing the reliability / internal consistency of a sum-score. Often interpreted as the mean correlation across all possible split-half alternate forms of the test.
#' @param x A data-frame or matrix of numerical values where rows are across-items within-respondent observation vectors, and columns are within-item across-respondents observation vectors.
#' @note Missing values are treated by passing \code{na.rm = TRUE} to the \code{var} function call.
#' @note Be aware that this function does not issue a warning if there are negative correlations between variables in the supplied data-set.
#' @return Cronbach's Alpha for the sum-score of supplied variables.
#' @references Cronbach, L.J. (1951). Coefficient alpha and the internal structure of tests. Psychometrika 16, 297--334. doi: 10.1007/BF02310555
#' @examples
#' # Generate some fictional data. Say 100 students take a 50-item long test
#' # where all items are equally difficult.
#' set.seed(1234)
#' p.success <- rBeta.4P(100, .25, .75, 5, 3)
#' for (i in 1:50) {
#'   if (i == 1) {
#'     rawdata <- matrix(nrow = 100, ncol = 50)
#'   }
#'   rawdata[, i] <- rbinom(100, 1, p.success)
#' }
#' # To calculate Cronbach's Alpha for this test:
#' cba(rawdata)
#' @export
cba <- function(x) {
  (base::ncol(x) / (base::ncol(x) - 1)) *
    (1 - (base::sum(base::diag(stats::var(x, na.rm = TRUE))) /
            base::sum(stats::var(x, na.rm = TRUE))))
}

#' Estimate Beta true-score distribution based on observed-score raw-moments and the effective test length.
#'
#' @description Estimator for the Beta true-score distribution shape-parameters from the observed-score distribution and Livingston and Lewis' effective test length. Returns a list with entries representing the lower- and upper shape parameters (l and u), and the shape parameters (alpha and beta) of the four-parameters beta distribution.
#' @param x Vector of observed-scores.
#' @param min The minimum possible score to attain on the test.
#' @param max The maximum possible score to attain on the test.
#' @param etl The value of Livingston and Lewis' effective test length. See ?ETL().
#' @param reliability Optional specification of the test-score reliability coefficient. If specified, overrides the input of the \code{etl} argument.
#' @param true.model The type of Beta distribution which is to be fit to the moments of the true-score distribution. Options are \code{"4P"} and \code{"2P"}, where "4P" refers to the four-parameter (with the same mean, variance, skewness, and kurtosis), and "2P" the two-parameter solution where both location-parameters are specified (with the same mean and variance).
#' @param failsafe Logical. Whether to revert to a failsafe two-parameter solution should the four-parameter solution contain invalid parameter estimates.
#' @param l If \code{failsafe = TRUE} or \code{true.model = "2P"}: The lower-bound of the Beta distribution. Default is 0 (i.e., the lower-bound of the Standard, two-parameter Beta distribution).
#' @param u If \code{failsafe = TRUE} or \code{true.model = "2P"}: The upper-bound of the Beta distribution. Default is 1 (i.e., the upper-bound of the Standard, two-parameter Beta distribution).
#' @param alpha If \code{failsafe = TRUE} or \code{true.model = "2P"}: The Alpha shape-parameter of the Beta distribution. Default is NA (i.e., estimate).
#' @param beta If \code{failsafe = TRUE} or \code{true.model = "2P"}: The Beta shape-parameter of the Beta distribution. Default is NA (i.e., estimate).
#' @param output Option to specify true-score distribution moments as output if the value of the output argument does not equal \code{"parameters"}.
#' @return A list with the parameter values of a four-parameter Beta distribution. "l" is the lower location-parameter, "u" the upper location-parameter, "alpha" the first shape-parameter, and "beta" the second shape-parameter.
#' @references Hanson, B. A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series. Retrieved from https://files.eric.ed.gov/fulltext/ED344945.pdf
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @references Rogosa, D. &  Finkelman, M. (2004). How Accurate Are the STAR Scores for Individual Students? – An Interpretive Guide. Retrieved from http://statweb.stanford.edu/~rag/accguide/guide04.pdf
#' @examples
#' # Generate some fictional data. Say 1000 individuals take a 100-item test
#' # where all items are equally difficult, and the true-score distribution
#' # is a four-parameter Beta distribution with location parameters l = .25,
#' # u = .75, alpha = 5, and beta = 3:
#' set.seed(12)
#' testdata <- rbinom(1000, 100, rBeta.4P(1000, .25, .75, 5, 3))
#'
#' # Since this test contains items which are all equally difficult, the true
#' # effective test length (etl) is the actual test length. I.e., etl = 100.
#' # To estimate the four-parameter Beta distribution parameters underlying
#' # the draws from the binomial distribution:
#' Beta.tp.fit(testdata, 0, 100, 100)
#'
#' # Imagine a case where the fitting procedure produces an impermissible
#' # estimate (e.g., l < 0 or u > 1).
#' set.seed(1234)
#' testdata <- rbinom(1000, 50, rBeta.4P(1000, .25, .75, 5, 3))
#' Beta.tp.fit(testdata, 0, 50, 50)
#'
#' # This example produced an l-value estimate less than 0. One way of
#' # dealing with such an occurrence is to revert to a two-parameter
#' # model, specifying the l and u parameters and estimating the
#' # alpha and beta parameters necessary to produce a Beta distribution
#' # with the same mean and variance as the estimated true-score distribution.
#'
#' # Suppose you have good theoretical reasons to fix the l parameter at a
#' # value of 0.25 (e.g., the test is composed of multiple-choice questions
#' # with four response-options, resulting in a 25% chance of guessing the
#' # correct answer). The l-parameter could be specified to this theoretically
#' # justified value, and the u-parameter could be specified to be equal to the
#' # estimate above (u = 0.7256552) as such:
#' Beta.tp.fit(testdata, 0, 50, 50, true.model = "2P", l = 0.25, u = 0.7256552)
#' @export
Beta.tp.fit <- function(x, min, max, etl, reliability = NULL, true.model = "4P", failsafe = FALSE, l = 0, u = 1, alpha = NA, beta = NA, output = "parameters") {
  if(output != "parameters") {
    moments <- list()
  }
  l.save <- l
  u.save <- u
  alpha.save <- alpha
  beta.save <- beta
  if (!is.null(reliability)) {
    etl <- ETL(base::mean(x), stats::var(x), min, max, reliability)
  }
  x <- (x - min) / (max - min) * etl
  m1 <- mean(x)
  m2 <- mean(x^2)
  m3 <- mean(x^3)
  m4 <- mean(x^4)
  tp.m1 <- mean(x) / etl
  tp.m2 <- (m2 - m1) / (etl * (etl - 1))
  tp.m3 <- (m3 - 3*m2 + 2*m1) / (etl * (etl - 1) * (etl - 2))
  tp.m4 <- (m4 - 6*m3 + 11*m2 - 6*m1) / (etl * (etl - 1) * (etl - 2) * (etl - 3))
  tp.s2 <- tp.m2 - tp.m1^2
  if (output != "parameters") {
    tp.s3 <- (tp.m3 - 3 * tp.m1 * tp.m2 + 2 * tp.m1^3)
    tp.s4 <- (tp.m4 - 4 * tp.m1 * tp.m3 + 6 * tp.m1^2 * tp.m2 - 3 * tp.m1^4)
  }
  tp.g3 <- (tp.m3 - 3 * tp.m1 * tp.m2 + 2 * tp.m1^3) / (sqrt(tp.s2)^3)
  tp.g4 <- (tp.m4 - 4 * tp.m1 * tp.m3 + 6 * tp.m1^2 * tp.m2 - 3 * tp.m1^4) / (sqrt(tp.s2)^4)
  if (output == "parameters") {
    if (true.model == "4P" | true.model == "4p") {
      params <- Beta.4p.fit(mean = tp.m1, variance = tp.s2, skewness = tp.g3, kurtosis = tp.g4)
      l <- params$l
      u <- params$u
      alpha <- params$alpha
      beta <- params$beta
    }
    if ((true.model == "2P" | true.model == "2p") | (failsafe & (any(is.na(c(l, u, alpha, beta))) | (l < 0 | u > 1 | alpha <= 0 | beta <= 0)))) {
      if ((failsafe & (any(is.na(c(l, u, alpha, beta))) | (l < 0 | u > 1 | alpha <= 0 | beta <= 0)))) {
        warning(paste("Fail-safe engaged: l = ", l, ", u = ", u, ", alpha = ", alpha, ", beta = ", beta,
                      ". Finding permissible solution for the true-score distribution in accordance with specifications.", sep = ""))
      }
      if ((true.model != "2p" & true.model != "2P") & is.na(l.save)) {
        l <- 0
        } else {
          l <- l.save
        }
      if ((true.model != "2p" & true.model != "2P") & is.na(u.save)) {
        u <- 1
        } else {
          u <- u.save
        }
      alpha <- alpha.save
      beta <- beta.save
      if (!is.na(alpha) & !is.na(beta) & is.na(l) & is.na(u)) {
        l <- LABMSU(alpha = alpha, beta = beta, mean = tp.m1, variance = tp.s2)
        u <- UABMSL(alpha = alpha, beta = beta, mean = tp.m1, variance = tp.s2)
      }
      if (!is.na(alpha) & !is.na(beta) & is.na(l) & !is.na(u)) {
        l <- LABMSU(alpha = alpha, beta = beta, mean = tp.m1, variance = tp.s2, u = u)
      }
      if (!is.na(alpha) & !is.na(beta) & !is.na(l) & is.na(u)) {
        u <- UABMSL(alpha = alpha, beta = beta, mean = tp.m1, variance = tp.s2, l = l)
      }
      if (!is.na(alpha) & is.na(beta) & !is.na(l) & !is.na(u)) {
        beta <- BMS(mean = tp.m1, variance = tp.s2, l = l, u = u, alpha = alpha)
      }
      if (is.na(alpha) & !is.na(beta) & !is.na(l) & !is.na(u)) {
        alpha <- AMS(mean = tp.m1, variance = tp.s2, l = l,u = u, beta = beta)
      }
      if (!is.na(alpha) & is.na(beta) & !is.na(l) & !is.na(u)) {
        beta <- BMS(mean = tp.m1, variance = tp.s2, l = l, u = u, alpha = alpha)
      }
      if (is.na(alpha) & is.na(beta) & !is.na(l) & !is.na(u)) {
        alpha <- AMS(mean = tp.m1, variance = tp.s2, l = l, u = u, beta = NULL)
        beta <- BMS(mean = tp.m1, variance = tp.s2, l = l, u = u, alpha = NULL)
      }
    }
    return(list("l" = l, "u" = u, "alpha" = alpha, "beta" = beta, "etl" = etl))
  } else {
    moments[["Raw"]] <- list(tp.m1, tp.m2, tp.m3, tp.m4)
    moments[["Central"]] <- list(0, tp.s2, tp.s3, tp.s4)
    moments[["Standardized"]] <- list(0, 1, tp.g3, tp.g4)
    return(moments)
  }
}

#' Descending (falling) factorial.
#'
#' @description Calculate the descending (or falling) factorial of a value \code{x} of order \code{r}.
#' @param x A value for which the descending factorial is to be calculated.
#' @param r The power \code{x} is to be raised to.
#' @return The descending factorial of value \code{x} raised to the \code{r} power.
#' @export
#' @note This function implements the descending factorial by means of the Gamma distribution. As such, \code{x} does not have to be an integer. However, \code{x} cannot be a negative integer.
#' @examples
#' # To calculate the 4th descending factorial for a value (e.g., 3.14):
#' dfac(x = 3.14, r = 4)
#'
#' # To calculate the 5th descending factorial for values 3.14, 2.72, and 0.58:
#' dfac(x = c(3.14, 2.72, 0.58), r = 5)
dfac <- function(x, r) {
  gamma(x + 1) / gamma(x - r + 1)
}

#' Ascending (rising) factorial.
#'
#' @description Calculate the ascending (or rising) factorial of a value \code{x} of order \code{r}.
#' @param x A value for which the ascending factorial is to be calculated.
#' @param r The power \code{x} is to be raised to.
#' @return The ascending factorial of value \code{x} raised to the \code{r} power.
#' @note This function implements the ascending factorial by means of the Gamma distribution. As such, \code{x} does not have to be an integer. However, \code{x} cannot be a negative integer.
#' @export
#' @examples
#' # To calculate the 4th ascending factorial for a value (e.g., 3.14):
#' afac(x = 3.14, r = 4)
#'
#' # To calculate the 5th ascending factorial for values 3.14, 2.72, and 0.58:
#' afac(x = c(3.14, 2.72, 0.58), r = 5)
afac <- function(x, r) {
  gamma(x + r) / gamma(x)
}

#' Proportional true-score distribution raw moments from Livingston and Lewis' effective test-score and effective test-length.
#'
#' @description An implementation of Lords (1965, p. 265) equation 37 for estimating the raw moments of the true-score distribution, modified to function for the Livingston and Lewis approach.
#' @param x The effective test-score of test-takers.
#' @param r The moment-order that is to be calculated (where 1 is the mean, 2 is the raw variance, 3 is the raw skewness, etc.).
#' @param n The effective test-length.
#' @references Lord, F. M. (1965). A strong true-score theory, with applications. Psychometrika. 30(3). pp. 239--270. doi: 10.1007/BF02289490
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
#' @examples
#' # Examine the raw moments of the underlying Beta distribution that is to provide the basis for
#' # observed-scores:
#' betamoments(alpha = 5, beta = 3, l = .25, u = .75, types = "raw")
#'
#' # Generate observed-scores from true-scores by passing the true-scores as binomial probabilities
#' # for the rbinom function.
#' set.seed(1234)
#' obs.scores <- rbinom(1000, 100, rBeta.4P(1000, .25, .75, 5, 3))
#' # Examine the raw moments of the observed-score distribution.
#' observedmoments(obs.scores, type = "raw")
#'
#' # First four estimated raw moment of the proportional true-score distribution from the observed-
#' # score distribution. As all items are equally difficult, the effective test-length is equal to
#' # the actual test-length.
#' tsm(x = obs.scores, r = 1, n = 100)
#' tsm(x = obs.scores, r = 2, n = 100)
#' tsm(x = obs.scores, r = 3, n = 100)
#' tsm(x = obs.scores, r = 4, n = 100)
#' # Which is fairly close to the true raw moments of the proportional true-score distribution
#' # calculated above.
tsm <- function(x, r, n) {
  mean(dfac(x + r, r)) / mean(dfac(n + r-2, r-2)) / dfac(n + r, 2)
}
