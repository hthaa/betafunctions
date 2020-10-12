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
#' @param true.model The probability distribution to be fitted to the moments of the true-score distribution. Options are \code{"4P"} (default) and \code{"2P"}, referring to four- and two-parameter Beta distributions. The "4P" method produces a four-parameter Beta distribution with the same first four moments (mean, variance, skewness, and kurtosis) as the estimated true-score distribution, while the "2P" method produces a two-parameter Beta distribution with the first two moments (mean and variance) as the estimated true-score distribution.
#' @param error.model The probability distribution to be used for producing the sampling distributions at different points of the true-score scale. Options are \code{binomial} and \code{beta}. The binomial distribution is discrete, and is the distribution used originally by Livingston and Lewis. Use of the binomial distribution involves a rounding of the effective test length to the nearest integer value. The Beta distribution is continuous, and does not involve rounding of the effective test length.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves (see documentation for the \code{LL.ROC()} function).
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
#' # supplied manually along with the effective test length (see documentation
#' # for the ETL() function), foregoing the need for actual data. The list
#' # entries must be named. "l" is the lower-bound and "u" the upper-bound
#' # location parameters of the true-score distribution, and "alpha" and "beta"
#' # the shape parameters, and "etl" for the effective test-length..
#' trueparams <- list("l" = 0.25, "u" = 0.75, "alpha" = 5, "beta" = 3, "etl" = 50)
#' LL.CA(x = trueparams, reliability = .7, cut = 50, min = 0, max = 100)
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing.
#' @export
LL.CA <- function(x = NULL, reliability, cut, min = 0, max = 1, true.model = "4P", error.model = "binomial", truecut = NULL, output = c("accuracy", "consistency"), override = FALSE) {
  out <- base::list()
  if (class(x) != "list") {
    if ((base::min(x) < min) | (base::max(x) > max)) {
      warning(paste("Observed values not within the specified [", min, ", ", max, "] bounds (observed min = ", min(x), ", observed max = ", max(x), ").", sep = ""))
    }
    N <- ETL(base::mean(x), stats::var(x), l = min, u = max, reliability = reliability)
    params <- Beta.tp.fit(x, min = min, max = max, etl = N, true.model = true.model, failsafe = !override)
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
  if (any(output == "accuracy") | any(output == "Accuracy") | any(output == "ca") | any(output == "CA") | any(output == "a") | any(output == "A")) {
    if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
      p.tp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(floor(cut * N), N, x, lower.tail = FALSE) }, lower = truecut, upper = 1)$value
      p.fp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(floor(cut * N), N, x, lower.tail = FALSE) }, lower = 0, upper = 1)$value - p.tp
      p.ff <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(floor(cut * N), N, x) }, lower = truecut, upper = 1)$value
      p.tf <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(floor(cut * N), N, x) }, lower = 0, upper = 1)$value - p.ff
    }
    if (error.model == "beta" | error.model == "Beta") {
      p.tp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x), lower.tail = FALSE) }, lower = truecut, upper = 1)$value
      p.fp <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x), lower.tail = FALSE) }, lower = 0, upper = 1)$value - p.tp
      p.ff <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x)) }, lower = truecut, upper = 1)$value
      p.tf <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbeta(cut, N * x, N * (1 - x)) }, lower = 0, upper = 1)$value - p.ff
    }
    camat <- base::matrix(nrow = 2, ncol = 2, dimnames = list(c("True", "False"), c("Fail", "Pass")))
    camat["True", "Fail"] <- p.tf
    camat["True", "Pass"] <- p.tp
    camat["False", "Fail"] <- p.ff
    camat["False", "Pass"] <- p.fp
    out[["confusionmatrix"]] <- camat
    out[["classification.accuracy"]] <- caStats(camat[1, 1], camat[1, 2], camat[2, 1], camat[2, 2])
  }
  if (any(output == "consistency") | any(output == "Consistency" ) | any(output == "cc") | any(output == "CC") | any(output == "c") | any(output == "C")) {
    if (error.model == "binomial" | error.model == "Binomial" | error.model == "binom" | error.model == "Binom") {
      p.ii <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (1 - stats::pbinom(floor(truecut * N), N, x, lower.tail = FALSE))^2 }, lower = 0, upper = 1)$value
      p.ij <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * (1 - stats::pbinom(floor(truecut * N), N, x, lower.tail = FALSE)) * stats::pbinom(truecut * N, N, x, lower.tail = FALSE) }, lower = 0, upper = 1)$value
      p.jj <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(floor(truecut * N), N, x, lower.tail = FALSE)^2 }, lower = 0, upper = 1)$value
      p.ji <- stats::integrate(function(x) { dBeta.4P(x, params$l, params$u, params$alpha, params$beta) * stats::pbinom(floor(truecut * N), N, x, lower.tail = FALSE) * (1 - stats::pbinom(truecut * N, N, x, lower.tail = FALSE)) }, lower = 0, upper = 1)$value
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
LL.ROC <- function(x = NULL, reliability, min = 0, max = 1, truecut, true.model = "4P", error.model = "Binomial", AUC = FALSE, maxJ = FALSE, raw.out = FALSE, grainsize = 100) {
  oldpar <- graphics::par(no.readonly = TRUE)
  base::on.exit(graphics::par(oldpar))
  for (i in 1:(grainsize + 1)) {
    if (i == 1) {
      cuts <- seq(min, max, (max - min) / grainsize)
      outputmatrix <- matrix(nrow = grainsize + 1, ncol = 4)
      outputmatrix[, 4] <- cuts
    }
    axval <- LL.CA(x = x, min = min, max = max, reliability = reliability, cut = cuts[i], truecut = truecut, true.model = true.model, error.model = error.model, output = "a")$classification.accuracy
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

#' Estimate Beta true-score distribution based on observed-score raw-moments and the effective test length.
#'
#' @description Estimator for the Beta true-score distribution shape-parameters from the observed-score distribution and Livingston and Lewis' effective test length. Returns a list with entries representing the lower- and upper shape parameters (l and u), and the shape parameters (alpha and beta) of the four-parameters beta distribution.
#' @param x Vector of observed-scores.
#' @param min The minimum possible score to attain on the test.
#' @param max The maximum possible score to attain on the test.
#' @param etl The value of Livingston and Lewis' effective test length. See ?ETL().
#' @param true.model The type of Beta distribution which is to be fit to the moments of the true-score distribution. Options are \code{"4P"} and \code{"2P"}, where "4P" refers to the four-parameter (with the same mean, variance, skewness, and kurtosis) and "2P" the two-parameter solution (with the same mean and variance).
#' @param failsafe Logical. Whether to revert to a failsafe two-parameter solution should the four-parameter solution contain invalid parameter estimates.
#' @return A list with the parameter values of a four-parameter Beta distribution. "l" is the lower location-parameter, "u" the upper location-parameter, "alpha" the first shape-parameter, and "beta" the second shape-parameter.
#' @note This estimator is based on the S-Plus code provided by Rogosa and Finkelman (2004). It includes an option for implementing a failsafe should the four-parameter solution be invalid (e.g., l < 0 or u > 1, alpha < 1 or beta < 1).
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
#' @export
Beta.tp.fit <- function(x, min, max, etl, true.model = "4P", failsafe = FALSE, output = "parameters") {
  x <- (x - min) / (max - min) * etl
  mu1 <- sum(x) / length(x)
  mu2 <- sum(x^2) / length(x)
  mu3 <- sum(x^3) / length(x)
  mu4 <- sum(x^4) / length(x)
  tp.m1 <- mu1 / etl
  tp.m2 <- (mu2 - mu1)/(etl * (etl - 1))
  tp.m3 <- (mu3 - 3 * mu2 + 2 * mu1)/(etl * (etl - 1) * (etl - 2))
  tp.m4 <-  (mu4 - 6 * mu3 + 11 * mu2 - 6 * mu1) / (etl * (etl - 1) * (etl - 2) * (etl - 3))
  tp.s2 <- tp.m2 - tp.m1^2
  tp.g3 <- (tp.m3 - 3 * tp.m1 * tp.m2 + 2 * tp.m1^3) / (sqrt(tp.s2)^3)
  tp.g4 <- (tp.m4 - 4 * tp.m1 * tp.m3 + 6 * tp.m1^2 * tp.m2 - 3 * tp.m1^4) / (sqrt(tp.s2)^4)
  if (output == "parameters") {
    if (true.model == "4P" | true.model == "4p") {
      phi <- (6 * (tp.g4 - tp.g3^2 - 1)) / (6 + 3 * tp.g3^2 - 2 * tp.g4)
      if (tp.g3 < 0) {
        alpha <- (phi / 2) * (1 + sqrt(1 - (24 * (phi + 1) / (tp.g4 * (phi + 2) * (phi + 3) - 3 * (phi - 6) * (phi + 1)))))
        beta <- (phi / 2) * (1 - sqrt(1 - (24 * (phi + 1) / (tp.g4 * (phi + 2) * (phi + 3) - 3 * (phi - 6) * (phi + 1)))))
      } else {
        beta <- (phi / 2) * (1 + sqrt(1 - (24 * (phi + 1) / (tp.g4 * (phi + 2) * (phi + 3) - 3 * (phi - 6) * (phi + 1)))))
        alpha <- (phi / 2) * (1 - sqrt(1 - (24 * (phi + 1) / (tp.g4 * (phi + 2) * (phi + 3) - 3 * (phi - 6) * (phi + 1)))))
      }
      l <- tp.m1 - (alpha * sqrt(tp.s2 * (alpha + beta + 1)) / sqrt(alpha * beta))
      u <- tp.m1 + (beta * sqrt(tp.s2 * (alpha + beta + 1)) / sqrt(alpha * beta))
      if (failsafe & (any(is.na(c(l, u, alpha, beta))) | (l < 0 | u > 1 | (alpha < 1 & beta < 1)))) {
        warning(paste("Fail-safe engaged: l = ", l, ", u = ", u, ", alpha = ", alpha, ", beta = ", beta, ". Reverting to a two-parameter solution for the true-score distribution.", sep = ""))
        alpha <- AMS(tp.m1, tp.s2)
        beta <- BMS(tp.m1, tp.s2)
        l <- 0
        u <- 1
      }
    }
    if (true.model == "2P" | true.model == "2p") {
      alpha <- AMS(tp.m1, tp.s2)
      beta <- BMS(tp.m1, tp.s2)
      l <- 0
      u <- 1
    }
    return(list("l" = l, "u" = u, "alpha" = alpha, "beta" = beta))
  } else {
    return(list("mean" = m1, "variance" = s2, "skewness" = g3, "kurtosis" = g4))
  }
}
