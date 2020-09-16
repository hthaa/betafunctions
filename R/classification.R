#' Livingston and Lewis' "Effective Test Length".
#'
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param l The lower-bound of the observed-score distribution. Default is 0 (assuming observed scores represent proportions).
#' @param u The upper-bound of the observed-score distribution. Default is 1 (assuming observed scores represent proportions).
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
#' ETL(mean = mean(testdata), variance = var(testdata), l = 0, u = 100,
#' reliability = .7)
#' @export
ETL <- function(mean, variance, l = 0, u = 1, reliability) {
  ((mean - l) * (u - mean) - (reliability * variance)) / (variance * (1 - reliability))
}

#' An Implementation of the Livingston and Lewis (1995) Approach to Estimate Classification Consistency and Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification consistency and accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to the four-parameter beta distribution, and errors of measurement to the binomial distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores for which a beta-distribution is to be fitted, or a list of pre-defined true-score distribution parameter values. If a list is provided, the list entries must be named after the parameters: \code{l} and \code{u} for the location parameters, and \code{alpha} and \code{beta} for the shape parameters.
#' @param reliability The observed-score squared correlation (i.e., proportion of shared variance) with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0 (assuming \code{x} represent proportions).
#' @param max The maximum value possible to attain on the test. Default is 1 (assuming \code{x} represent proportions).
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @param error.model The probability distribution to be used for producing the sampling distributions at different points of the true-score scale. Options are \code{binomial} and \code{beta}. The binomial distribution is discrete, and is the distribution used originally by Livingston and Lewis. Use of the binomial distribution involves a rounding of the effective test length to the nearest integer value. The Beta distribution is continuous, and does not involve rounding of the effective test length.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves.
#' @param output Character vector indicating which types of statistics (i.e, accuracy and/or consistency) are to be computed and included in the output. Permissible values are \code{"accuracy"} and \code{"consistency"}.
#' @param override Logical value indicating whether to override the automatic default to the two-parameter Beta true-score distribution if the four-parameter fitting procedure produces impermissible parameter estimates. Default is \code{FALSE}.
#' @return A list containing the estimated parameters necessary for the approach (i.e., the effective test-length and the beta distribution parameters), the confusion matrix containing estimated proportions of true/false pass/fail categorizations for a test, diagnostic performance statistics, and / or a classification consistency matrix and indices. Accuracy output includes a confusion matrix and diagnostic performance indices, and consistency output includes a consistency matrix and consistency indices \code{p} (expected proportion of agreement between two independent test administrations), \code{p_c} (proportion of agreement on two independent administrations expected by chance alone), and \code{Kappa} (Cohen's Kappa).
#' @note It should be noted that this implementation differs from the original articulation of Livingston and Lewis (1995) in some respects. First, the procedure includes a number of diagnostic performance (accuracy) indices which the original procedure enables but that were not included. Second, the possibility of employing a two-parameter Beta error distribution in place of the binomial error distribution is not part of the original procedure. Third, the way consistency is calculated differs substantially from the original articulation of the procedure, which made use of a split-half approach. Rather, this implementation uses the approach to calculating classification consistency outlined by Hanson (1991).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and
#' # that the reliability of this test was estimated to 0.7. To estimate and
#' # retrieve the estimated parameters, confusion matrix, consistency and
#' # accuracy statistics using LL.CA():
#' LL.CA(x = testdata, reliability = .7, cut = 50, min = 0, max = 100)
#'
#' # Alternatively to supplying scores to which a true-score distribution is
#' # to be fit, a list with true-score distribution parameter values can be
#' # supplied manually, foregoing the need for actual data. The list entries
#' # must be named. "l" is the lower-bound and "u" the upper-bound location
#' # parameters of the true-score distribution, and "alpha" and "beta" the
#' # shape parameters.
#' trueparams <- list("l" = 0.25, "u" = 0.75, "alpha" = 5, "beta" = 3)
#' LL.CA(x = trueparams, reliability = .7, cut = 50, min = 0, max = 100)
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @export
LL.CA <- function(x = NULL, reliability, cut, min = 0, max = 1, error.model = "binomial", truecut = NULL, output = c("accuracy", "consistency"), override = FALSE) {
  out <- base::list()
  if (class(x) != "list") {
    if ((base::min(x) < min) | (base::max(x) > max)) {
      warning(paste("Observed values not within the specified [", min, ", ", max, "] bounds (observed min = ", min(x), ", observed max = ", max(x), ").", sep = ""))
    }
    x <- (x - min) / (max - min)
    params <- Beta.4p.fit(x)
  } else {
    params <- x
  }
  if (override == FALSE & class(x) != "list") {
    if ((params$l < 0) | (params$u > 1) | any(is.na(params))) {
      warning(paste("Improper solution for true-score distribution location parameters (l = ", params$l, ", u = ", params$u, "). Reverting to two-parameter solution.", sep = ""))
      params <- Beta.2p.fit(x)
    }
  }

  #Calculate mean and variance for the true-score distribution.
  TpMoments <- betamoments(params$alpha, params$beta, params$l, params$u, types = c("raw", "central"), orders = 2)
  #Estimate the effective test-length based on true-score distribution and reliability.
  N <- ETL(TpMoments[["raw"]][[1]], TpMoments[["central"]][[2]], reliability = reliability)
  #Establish cutscores on the proportional scale.
  if (base::is.null(truecut)) {
    truecut <- cut
  }
  cut <- (cut - min) / (max - min)
  truecut <- (truecut - min) / (max - min)
  #Calculate probabilities of producing passing and failing scores along the true-score distribution.
  if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
    N <- base::round(N)
  }
  out[["effectivetestlength"]] <- N
  out[["parameters"]] <- params
  if (any(output == "accuracy") | any(output == "Accuracy") | any(output == "ca") | any(output == "CA") | any(output == "a") | any(output == "A")) {
    if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
      p.tp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(cut * N, N, x, lower.tail = FALSE) }, lower = truecut, upper = 1)$value
      p.fp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(cut * N, N, x, lower.tail = FALSE) }, lower = 0, upper = 1)$value - p.tp
      p.ff <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(cut * N, N, x) }, lower = truecut, upper = 1)$value
      p.tf <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(cut * N, N, x) }, lower = 0, upper = 1)$value - p.ff
    }
    if (error.model == "beta" | error.model == "Beta") {
      p.tp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x), lower.tail = FALSE) }, lower = truecut, upper = 1)$value
      p.fp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x), lower.tail = FALSE) }, lower = 0, upper = 1)$value - p.tp
      p.ff <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x)) }, lower = truecut, upper = 1)$value
      p.tf <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x)) }, lower = 0, upper = 1)$value - p.ff
    }
    # Calculate confusion matrix.
    camat <- base::matrix(nrow = 2, ncol = 2, dimnames = list(c("True", "False"), c("Fail", "Pass")))
    camat["True", "Fail"] <- p.tf
    camat["True", "Pass"] <- p.tp
    camat["False", "Fail"] <- p.ff
    camat["False", "Pass"] <- p.fp
    out[["confusionmatrix"]] <- camat
    # Calculate classification accuracy indices.
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(output == "consistency") | any(output == "Consistency" ) | any(output == "cc") | any(output == "CC") | any(output == "c") | any(output == "C")) {
    # Calculate classification consistency indices.
    if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
      p.ii <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (1 - stats::pbinom(truecut * N, N, x, lower.tail = FALSE))^2 }, lower = 0, upper = 1)$value
      p.ij <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (1 - stats::pbinom(truecut * N, N, x, lower.tail = FALSE)) * stats::pbinom(truecut * N, N, x, lower.tail = FALSE) }, lower = 0, upper = 1)$value
      p.jj <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(truecut * N, N, x, lower.tail = FALSE)^2 }, lower = 0, upper = 1)$value
      p.ji <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(truecut * N, N, x, lower.tail = FALSE) * (1 - stats::pbinom(truecut * N, N, x, lower.tail = FALSE)) }, lower = 0, upper = 1)$value
    }
    if (error.model == "beta" | error.model == "Beta") {
      p.ii <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (1 - stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = FALSE))^2 }, lower = 0, upper = 1)$value
      p.ij <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (1 - stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = FALSE)) * stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = FALSE) }, lower = 0, upper = 1)$value
      p.jj <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = FALSE)^2 }, lower = 0, upper = 1)$value
      p.ji <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = FALSE) * (1 - stats::pbeta(truecut, N * x, N * (1 - x), lower.tail = FALSE)) }, lower = 0, upper = 1)$value
    }
    ccmat <- base::matrix(nrow = 2, ncol = 2, dimnames = list(c("i", "j"), c("i", "j")))
    ccmat["i", "i"] <- p.ii
    ccmat["i", "j"] <- p.ij
    ccmat["j", "i"] <- p.ji
    ccmat["j", "j"] <- p.jj
    out[["consistencymatrix"]] <- ccmat
    # Calculate classification consistency indices.
    out[["classification.consistency"]] <- ccStats(ccmat["i", "i"], ccmat["i", "j"], ccmat["j", "i"], ccmat["j", "j"])
  }
  return(out)
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
#' caStats(tp = cmat["True", "Fail"], tn = cmat["True", "Pass"],
#' fp = cmat["False", "Fail"], fn = cmat["False", "Pass"])
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
#' # To estimate and retrieve consistency statistics using caStats(),
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
LL.ROC <- function(x = NULL, reliability, min = 0, max = 1, truecut, AUC = FALSE, maxJ = FALSE, raw.out = FALSE, grainsize = 100) {
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  for (i in 1:(grainsize + 1)) {
    if (i == 1) {
      cuts <- seq(min, max, (max - min) / grainsize)
      outputmatrix <- matrix(nrow = grainsize + 1, ncol = 4)
      outputmatrix[, 4] <- cuts
    }
    axval <- LL.CA(x = x, min = min, max = max, reliability = reliability, cut = cuts[i], truecut = truecut, output = "a")$classification.accuracy
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
    graphics::legend("bottomright", bty = "n", cex = 1.5, legend = paste("AUC =", round(AUC(outputmatrix[, 1], outputmatrix[, 2]), 3)))
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
