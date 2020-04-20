#' Livingston and Lewis' "Effective Test Length".
#'
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param l The lower-bound of the observed-score distribution.
#' @param u The upper-bound of the observed-score distribution.
#' @param reliability The reliability of the observed scores (correlation with true-score distribution).
#' @return An estimate of the effective length of a test, given the stability of the observations it produces.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
ETL <- function(mean, variance, l = 0, u = 1, reliability) {
  ((mean - l) * (u - mean) - (reliability * variance)) / (variance * (1 - reliability))
}

#' An Implementation of the Livingston and Lewis (1995) Approach to Estimate Classification Accuracy based on Observed Test Scores and Test Reliability.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification accuracy, which by employing a compound beta-binomial distribution assumes that true-scores conform to four-parameter beta distributions, and errors of measurement binomial distribution distribution. Under these assumptions, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and test reliability.
#' @param x A vector of observed scores for which a beta-distribution is to be fitted.
#' @param reliability The observed-score squared correlation with the true-score.
#' @param min The minimum value possible to attain on the test. Default is 0.
#' @param max The maximum value possible to attain on the test. Default is 1.
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @param error.model The probability distribution to be used for producing the sampling distributions at different points of the true-score scale. Options are \code{beta} and \code{binomial}. The binomial distribution is discrete, and is the distribution used originally by Livingston and Lewis. Use of the binomial distribution involves a rounding of the effective test length to the nearest integer value. The Beta distribution is continuous, and does not involve rounding of the effective test length..
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curves.
#' @param grainsize The size of the steps for which probabilities along the score distribution are to be calculated. Default is .001 (1001 points).
#' @return A list containing the estimated parameters necessary for the approach, as well as the confusion matrix estimating the proportion of true/false pass/fail categorizations for a test, given a specific distribution of observed scores.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a maximum score of 100 and a minimum score of 0.
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and that the reliability of this test was estimated to 0.7. To estimate and retrieve the necessary parameters and the confusion matrix with LL.CA():
#' LL.CA(x = testdata, reliability = .7, cut = 50, min = 0, max = 100)
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
LL.CA <- function(x = NULL, reliability, cut, min = 0, max = 1, error.model = "binomial", truecut = NULL, grainsize = .001) {
  x <- (x - min) / (max - min)
  params <- Beta.4p.fit(x)
  if (params$l < 0) {
    warning("Improper solution for lower-bound estimate of true-score distribution (< 0). Reverting to two-parameter solution.")
    params$alpha <- AMS(base::mean(x), stats::var(x))
    params$beta <- BMS(base::mean(x), stats::var(x))
    params$l <- 0
    params$u <- 1
  }
  if (params$u > 1) {
    warning("Improper solution for upper-bound estimate of true-score distribution (> 1). Reverting to two-parameter solution.")
    params$alpha <- AMS(base::mean(x), stats::var(x))
    params$beta <- BMS(base::mean(x), stats::var(x))
    params$l <- 0
    params$u <- 1
  }

  #Calculate mean and variance for the true-score distribution.
  TpMoments <- betamoments(params$alpha, params$beta, params$l, params$u, types = c("raw", "central"), orders = 2)

  #Estimate the effective test-length based on true-score distribution and reliability.
  N <- ETL(TpMoments[["raw"]][[1]], TpMoments[["central"]][[2]], reliability = reliability)

  #Establish cutscores on the proportional scale.
  if (base::is.null(truecut)) {
    truecut <- cut
  }
  cut <- cut / max
  truecut <- truecut / max

  #Generate density distribution along the proportional scale.
  xaxis <- base::seq(0, 1, grainsize)
  sumdens <- base::sum(dBeta.4P(xaxis, params$l, params$u, params$alpha, params$beta))
  density <- dBeta.4P(xaxis, params$l, params$u, params$alpha, params$beta) / sumdens


  #Calculate probabilities of producing passing and failing scores along the true-score distribution.
  if (error.model == "beta") {
    p.pass <- stats::pbeta(cut, xaxis * N, (1 - xaxis) * N, lower.tail = FALSE)
  }
  if (error.model == "binomial") {
    N <- base::round(N)
    p.pass <- stats::pbinom(cut * N, N, xaxis, lower.tail = FALSE)
  }
  p.fail <- 1 - p.pass

  # Calculate proportions of true and false positives and negatives.
  p.tf <- p.fail[which(xaxis < truecut)] * density[which(xaxis < truecut)]
  p.fp <- p.pass[which(xaxis < truecut)] * density[which(xaxis < truecut)]
  p.ff <- p.fail[which(xaxis >= truecut)] * density[which(xaxis >= truecut)]
  p.tp <- p.pass[which(xaxis >= truecut)] * density[which(xaxis >= truecut)]

  # Calculate confusion matrix.
  cmat <- base::matrix(nrow = 2, ncol = 2)
  rownames(cmat) <- base::c("True", "False")
  colnames(cmat) <- base::c("Fail", "Pass")
  cmat["True", "Fail"] <- base::sum(p.tf)
  cmat["True", "Pass"] <- base::sum(p.tp)
  cmat["False", "Fail"] <- base::sum(p.ff)
  cmat["False", "Pass"] <- base::sum(p.fp)
  return(base::list("effectivetestlength" = N, "parameters" = params, "confusionmatrix" = cmat))
}

#' Classification Accuracy Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the certainty of classifications based on tests.
#' @param tp The frequency or rate of true-positive classifications.
#' @param tn The frequency or rate of true-negative classifications.
#' @param fp The frequency or rate of false-positive classifications.
#' @param fn The frequency or rate of false-negative classifications.
#' @return A list of diagnostic performance statistics based on true/false positive/negative statistics. Specifically, the sensitivity, specificity, positive likelihood ratio (LR.pos), negative likelihood ratio (LR.neg), diagnostic odds ratio (DOR), positive predictive value (PPV), negative predictive value (NPV), Youden's J. (Youden.J), and Accuracy.
#' @references Glas et al. (2003). The Diagnostic Odds Ratio: A Single Indicator of Test Performance, Journal of Clinical Epidemiology, 1129-1135, 56(11). doi: 10.1016/S0895-4356(03)00177-X
#' @export
caStats <- function(tp, tn, fp, fn) {
  sensitivity <-  tp / (tp + fn)
  specificity <-  tn / (tn + fp)
  plr <-          sensitivity / (1 - specificity)
  nlr <-          (1 - sensitivity) / specificity
  dor <-          plr / nlr
  ppv <-          tp / (tp + fp)
  npv <-          tn / (tn + fn)
  accuracy <-     (tp + tn) / (tp + tn + fp + fn)
  J <-            (sensitivity + specificity) - 1
  base::list("Sensitivity" = sensitivity, "Specificity" = specificity,
             "LR.pos" = plr, "LR.neg" = nlr, "DOR" = dor,
             "PPV" = ppv, "NPV" = npv,
             "Youden.J" = J, "Accuracy" = accuracy)
}

#' ROC curves for the Livingston and Lewis approach.
#'
#' @description Generate a ROC curve plotting the false-positive rate against the true-positive rate at different cut-off values across the observed proportion-score scale.
#' @param x A vector of observed results.
#' @param min The minimum possible value to attain on the observed-score scale. Default is 0 (assuming \code{x} represent proportions).
#' @param max The maximum possible value to attain on the observed-score scale. Default is 1 (assuming \code{x} represent proportions).
#' @param reliability The reliability coefficient of the test.
#' @param truecut The true point along the x-scale that marks the categorization-threshold.
#' @param AUC Calculate and include the area under the curve? Default is FALSE.
#' @param maxJ Mark the point along the curve where Youden's J statistic is maximized? Default is FALSE.
#' @param raw.out Give raw coordinates as output rather than plot? Default is FALSE.
#' @return A plot tracing the ROC curve for the test, or matrix of coordinates if raw.out is TRUE.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a maximum score of 100 and a minimum score of 0.
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # Suppose the cutoff value for attaining a pass is 50 items correct, and that the reliability of this test was estimated to 0.7. To produce a plot with an ROC curve using LL.ROC(), along with the AUC statistics and the points at which Youden's J. is maximized:
#' LL.ROC(x = testdata, reliability = .7, truecut = 50, min = 0, max = 100, AUC = TRUE, maxJ = TRUE)
#' @export
LL.ROC <- function(x = NULL, reliability, min = 0, max = 1, truecut, AUC = FALSE, maxJ = FALSE, raw.out = FALSE) {
  for (i in 1:length(seq(0, 1, .001))) {
    if (i == 1) {
      cuts <- seq(min, max, (max - min) / 1000)
      outputmatrix <- matrix(nrow = length(seq(0, 1, .001)), ncol = 4)
    }
    cmat <- LL.CA(x = x, min = min, max = max, reliability = reliability, cut = cuts[i], truecut = truecut)$confusionmatrix
    axval <- caStats(cmat[1, 1], cmat[1, 2], cmat[2, 1], cmat[2, 2])
    outputmatrix[i, 1] <- 1 - axval$Specificity
    outputmatrix[i, 2] <- axval$Sensitivity
    outputmatrix[i, 3] <- axval$Youden.J
    outputmatrix[, 4] <- cuts
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
#' @param TPR Vector of True-Posiitive Rates.
#' @return A value representing the area under the ROC curve.
#' @note Script originally retrieved and modified from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html.
#' @export
AUC <- function(FPR, TPR) {
  dFPR <- base::c(diff(FPR), 0)
  dTPR <- base::c(diff(TPR), 0)
  base::sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}
