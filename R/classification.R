
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

#' An implementation of the Livingston and Lewis (1995) Method to Classification Accuracy and Consistency based on Test Scores.
#'
#' @description An implementation of what has been come to be known as the "Livingston and Lewis approach" to classification accuracy and consistency, which assumes that observed-scores, true-scores, and errors of measurement follow the four-parameter beta distribution. Under this assumption, the expected classification consistency and accuracy of tests can be estimated from observed outcomes and estimated test reliability.
#' @param x A vector of observed scores for which a beta-distribution is to be fitted.
#' @param reliability The observed-score correlation with the true-score.
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @return A confusion matrix estimating the proportion of true/false pass/fail categorizations for a test, given a specific distribution of observed scores.
#' @references Livinston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @export
LL.cac <- function(x = NULL, min = 0, max = 1, reliability, cut) {
  x <- (x - min) / (max - min)
  params <- Beta.4p.fit(x)
  x.moments <- observedmoments(x)
  mean <- x.moments[[1]][[1]]
  variance <- x.moments[[2]][[2]]
  N <- ETL(mean, variance, reliability = reliability)

  xaxis <- seq(params$l + .001, params$u - .001, .001)
  density <- dBeta.4P(xaxis, params$l, params$u, params$alpha, params$beta) /
    sum(dBeta.4P(xaxis, params$l, params$u, params$alpha, params$beta))

  #Calculate probabilities of producing passing and failing scores along the true-score distribution.
  p.pass <- pbeta(cut, xaxis * N, (1 - xaxis) * N, lower.tail = FALSE)
  p.fail <- 1 - p.pass

  p.tf <- p.fail[which(xaxis < cut)] * density[which(xaxis < cut)]
  p.fp <- p.pass[which(xaxis < cut)] * density[which(xaxis < cut)]

  p.ff <- p.fail[which(xaxis >= cut)] * density[which(xaxis >= cut)]
  p.tp <- p.pass[which(xaxis >= cut)] * density[which(xaxis >= cut)]

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
