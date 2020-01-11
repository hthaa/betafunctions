
k <- function(alpha, beta, reliability) {
  K <- alpha + beta
  mu <- alpha / K
  sigma2 <- alpha * beta / (alpha + beta + 1) * (alpha + beta)^2
  error <- sigma2 * (1 - reliability)
  K*((K - 1)*(sigma2 - error) - K*sigma2 + mu*(K - mu)) / 2*((mu*K - mu) - (sigma2 - error))
}

#' Livingston and Lewis' "Effective Test Length"
#'
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param l The lower-bound of the observed-score distribution.
#' @param u The upper-bound of the observed-score distribution.
#' @param reliability The reliability of the observed scores (correlation with true-score distribution).
#' @return An estimate of the effective length of a test, given the stability of the observations it produces.
#' @references Livinston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
ETL <- function(mean, variance, l = 0, u = 1, reliability) {
  ((mean - l) * (u - mean) - (reliability * variance)) / (variance * (1 - reliability))
}

#' An implementation of the Livingston and Lewis (1995) Method to Classification Accuracy based on Test Scores.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification accuracy, which assumes that observed-scores, true-scores, and errors of measurement follow the four-parameter beta distribution. Under this assumption, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and estimated test reliability.
#' @param x A vector of observed scores for which a beta-distribution is to be fitted.
#' @param reliability The observed-score correlation with the true-score.
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @param truecut Optional specification of a "true" cutoff. Useful for producing ROC curve values.
#' @param pdist The probability distribution to be used for pruducing the sampling distributions at different points of the true-score scale. Options are \code{beta} and \code{binomial}. The Beta distribution is continuous, while the binomial distribution is discrete. Use of the binomial distribution involves a rounding of the effective test length to the nearest integer value.
#' @return A confusion matrix estimating the proportion of true/false pass/fail categorizations for a test, given a specific distribution of observed scores.
#' @references Livinston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
LL.CA <- function(x = NULL, min = 0, max = 1, reliability, cut, truecut = NULL, pdist = "beta") {
  x <- (x - min) / (max - min)
  params <- Beta.4p.fit(x)
  if (params$l < 0) {
    warning("Improper solution for lower-bound estimate of true-score distribution (< 0). Reverting to two-parameter solution.")
    params$alpha <- AMS(mean(x), var(x))
    params$beta <- BMS(mean(x), var(x))
    params$l <- 0
    params$u <- 1
  }
  if (params$u > 1) {
    warning("Improper solution for upper-bound estimate of true-score distribution (> 1). Reverting to two-parameter solution.")
    params$alpha <- AMS(mean(x), var(x))
    params$beta <- BMS(mean(x), var(x))
    params$l <- 0
    params$u <- 1
  }
  mean <- mean(x)
  variance <- var(x)
  N <- ETL(mean, variance, reliability = reliability)
  if (is.null(truecut)) {
    truecut <- cut
  }
  xaxis <- seq(0, 1, .001)
  density <- dBeta.4P(xaxis, params$l, params$u, params$alpha, params$beta) /
    sum(dBeta.4P(xaxis, params$l, params$u, params$alpha, params$beta))

  #Calculate probabilities of producing passing and failing scores along the true-score distribution.
  if (pdist == "beta") {
    p.pass <- pbeta(cut, xaxis * N, (1 - xaxis) * N, lower.tail = FALSE)
  }
  if (pdist == "binomial") {
    N <- round(N)
    p.pass <- pbinom(cut * N, N, xaxis, lower.tail = FALSE)
  }
  p.fail <- 1 - p.pass

  # Calculate proportions of true and false positives and negatives.
  p.tf <- p.fail[which(xaxis < truecut)] * density[which(xaxis < truecut)]
  p.fp <- p.pass[which(xaxis < truecut)] * density[which(xaxis < truecut)]
  p.ff <- p.fail[which(xaxis >= truecut)] * density[which(xaxis >= truecut)]
  p.tp <- p.pass[which(xaxis >= truecut)] * density[which(xaxis >= truecut)]

  # Calculate confusion matrix.
  cmat <- matrix(nrow = 2, ncol = 2)
  rownames(cmat) <- c("True", "False")
  colnames(cmat) <- c("Fail", "Pass")
  cmat["True", "Fail"] <- sum(p.tf)
  cmat["True", "Pass"] <- sum(p.tp)
  cmat["False", "Fail"] <- sum(p.ff)
  cmat["False", "Pass"] <- sum(p.fp)
  return(list("effectivetestlength" = N, "parameters" = params, "confusionmatrix" = cmat))
}

#' Classification Accuracy Statistics.
#'
#' @description Provides a set of statistics often used for conveying information regarding the certainty of classifications based on tests.
#' @param tp The number or rate of true-positive classifications.
#' @param tn The number or rate of true-negative classifications.
#' @param fp The number or rate of false-positive classifications.
#' @param fn The number or rate of false-negative classifications.
#' @return A list of classification accuracy statistics based on true/false positive/negative statistics. Specifically, the sensitivity, specificity, positive likelihood ratio, negative likelihood ratio, positive predictive value, negative predictive value, and Youden's J.
#' @export
caStats <- function(tp, tn, fp, fn) {
  sensitivity <-  tp / (tp + fn)
  specificity <-  tn / (tn + fp)
  plr <-          sensitivity / (1 - specificity)
  nlr <-          (1 - sensitivity) / specificity
  ppv <-          tp / (tp + fp)
  npv <-          tn / (tn + fn)
  J <-            (sensitivity + specificity) - 1
  list("Sensitivity" = sensitivity, "Specificity" = specificity,
       "LR.pos" = plr, "LR.neg" = nlr,
       "PPV" = ppv, "NPV" = npv,
       "Youden.J" = J)
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
#' @return A plot tracing the ROC curve for the test.
#' @export
LL.ROC <- function(x = NULL, min = 0, max = 1, reliability, truecut, AUC = FALSE, maxJ = FALSE, raw.out = FALSE) {
  for (i in 1:length(seq(0, 1, .001))) {
    if (i == 1) {
      cuts <- seq(0, 1, .001)
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
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
  abline(h = seq(0, 1, .1), v = seq(0, 1, .1), col = "lightgrey", lty = "dotted")
  par(new = TRUE)
  plot(outputmatrix[, 1], outputmatrix[, 2], type = "l",
       xlab = "False-Positive Rate (1 - Specificity)",
       ylab = "True-Positive Rate (Sensitivity)",
       main = paste("ROC curve for true-cut equal to", truecut), lwd = 2,
       xlim = c(0, 1), ylim = c(0, 1))
  if (AUC) {
    legend("bottomright", bty = "n", cex = 1.5, legend = paste("AUC =", round(AUC(outputmatrix[, 1], outputmatrix[, 2]), 3)))
  }
  if (maxJ) {
    points(outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 1],
           outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 2], cex = 1.5, pch = 19)
    text(outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 1] + .025,
         outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 2] - .025,
         labels = paste("Maximum Youden's J. at cutoff = ",
                        round(outputmatrix[which(outputmatrix[, 3] == max(outputmatrix[, 3])), 4], 3),
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
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}
