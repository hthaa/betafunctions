#' Compute Moments of Two-to-Four Parameter Beta Probability Density Distributions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Standard Beta probability density distributions.
#' @param a The Alpha shape parameter of the PDD.
#' @param b The Beta shape parameter of the PDD.
#' @param l The first (lower) location parameter of a four-paramteer distribution.
#' @param u The second (upper) location parameter of a four-parameter distribution.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @references Hanson, B. A (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes. American College Testing Research Report Series.
#' @return A list of moment types, each a list of moment orders.
#' @export
betamoments <- function(a, b, l = 0, u = 1, mean = NULL, var = NULL, sd = NULL, types = c("raw", "central", "standardized"), orders = 4) {
  if (!is.null(mean) & (!is.null(var) | !is.null(sd))) {
    if (!is.null(var) & !is.null(sd)) {
      if (var != sd^2) {
        warning("Nonequivalent values of VAR and SD specified. Using VAR.")
      }
    }
    if (is.null(var) & !is.null(sd)) {
      var <- sd
    }
    a <- AMS(mean, var)
    b <- BMS(mean, var)
  }
  BETAMOMENTS <- rep(list(rep(list(NULL), orders)), length(types))
  TYPE <- 1
  if (any(types == "raw")) {
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- integrate(function(x) { l + (u - l) * dbeta(x, a, b) * x^i }, lower = 0, upper = 1)$value
    }
    names(BETAMOMENTS)[TYPE] <- "raw"
    TYPE <- TYPE + 1
  }
  if (any(types == "central")) {
    Mu <- integrate(function(x) { l + (u - l) * dbeta(x, a, b) * x }, lower = 0, upper = 1)$value
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- integrate(function(x) { (u - l)^2 * dbeta(x, a, b) * (x - Mu)^i },
                                            lower = 0, upper = 1)$value
    }
    names(BETAMOMENTS)[TYPE] <- "central"
    TYPE <- TYPE + 1
  }
  if (any(types == "standardized")) {
    Mu <- integrate(function(x) { l + (u - l) * dbeta(x, a, b) * x }, lower = 0, upper = 1)$value
    Sigma <- integrate(function(x) { (u - l)^2 * dbeta(x, a, b) * (x - Mu)^2 }, lower = 0, upper = 1)$value
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- integrate(function(x) { dbeta(x, a, b) * ((x - Mu)^i / sqrt(Sigma)^i) },
                                            lower = 0, upper = 1)$value
    }
    names(BETAMOMENTS)[TYPE] <- "standardized"
  }
  return(BETAMOMENTS)
}

observedmoments <- function(x, type = c("raw", "central", "standardized"),  orders = 4, correct = TRUE) {
  x <- na.omit(x)
  types <- 1
  momentorders <- list()
  if (any(type == "raw")) {
    mu <- list(rep(vector(length = 1), orders))
    for (i in 1:orders) {
      mu[i] <- sum(x^i)/length(x)
    }
    momentorders[[length(momentorders) + 1]] <- mu
    names(momentorders)[types] <- "raw"
    types <- types + 1
  }
  if (any(type == "central")) {
    sigma <- list(rep(vector(length = 1), orders))
    for (i in 1:orders) {
      if (correct) {
        sigma[i] <- sum((x - mean(x))^i)/(length(x) -
                                            1)
      }
      else {
        sigma[i] <- sum((x - mean(x))^i)/(length(x))
      }
    }
    momentorders[[length(momentorders) + 1]] <- sigma
    names(momentorders)[types] <- "central"
    types <- types + 1
  }
  if (any(type == "standardized")) {
    gamma <- list(rep(vector(length = 1), orders))
    for (i in 1:orders) {
      if (correct) {
        gamma[i] <- sum(((x - mean(x))^i)/sqrt(var(x))^i)/(length(x) -  2)
      }
      else {
        gamma[i] <- sum(((x - mean(x))^i)/sqrt(var(x))^i)/length(x)
      }
    }
    momentorders[[length(momentorders) + 1]] <- gamma
    names(momentorders)[types] <- "standardized"
  }
  return(momentorders)
}

#' Alpha shape parameter given mean and variance of a Standard Beta PDD.
#'
#' @description Calculates the Alpha value required to produce a Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Alpha shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
AMS <- function(mean, var, sd = NULL) {
  if ((!is.null(var) & !is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  if(((mean^2 - mean^3) / var) - mean <= 0) {
    warning("Parameter out of bounds (Alpha <= 0).")
  }
  return(((mean^2 - mean^3) / var) - mean)
}

#' Beta shape parameter given mean and variance of a Standard Beta PDD.
#'
#' @description Calculates the Beta value required to produce a Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
BMS <- function(mean, var, sd =NULL) {
  if ((!is.null(var) & !is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  if((mean * (1 - mean)^2) / var + mean - 1 <= 0) {
    warning("Parameter out of bounds (Beta <= 0).")
  }
  return(((mean * (1 - mean)^2) / var) + mean - 1)
}

#' Probability of some specific observation under the Standard Beta PDD with specific mean and variance.
#'
#' @description Calculates the probability of some specific observation falling under a specified interval  ([0, x] or [x, 1]) under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param q A specific point on the x-axis of the Standard Beta probability density distribution with a defined mean and variance.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lt Whether the density that should be considered is between the lower-end (i.e., [0 -> x]) or the higher-end of the distribution (i.e., [x -> 1]).
#' @return A value representing the probability of a random draw from the Standard Beta probability density distribution with a defined mean and variance being from one of two defined intervals (i.e., [0 -> x] or [x -> 1]).
#' @export
pBetaMS <- function(q, mean, var = NULL, sd = NULL, lt = TRUE) {
  if ((!is.null(var) & !is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  pbeta(q, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1, lower.tail = lt)
}

#' Density under a specific point of the Standard Beta PDD with specific mean and variance or standard deviation.
#'
#' @description Calculates the density under specific points of the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param x A specific point on the x-axis of the Standard Beta PDD.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta Shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
dBetaMS <- function(x, mean, var = NULL, sd = NULL) {
  if ((!is.null(var) & !is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  dbeta(x, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1)
}

#' Quantile containing specific proportion of the distribution, given a specific probability of the Standard Beta PDD with specific mean and variance or standard deviation.
#'
#' @description Calculates the quantile corresponding to a specific probability of some observation falling within the [0, x] (LT = TRUE) or [x, 1] (LT = FALSE) interval under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param p A value of probability marking the point of the Y-axis to correspond to the X-axis.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lt Logical. Specifies which end of the tail for which to calculate quantile. Default is TRUE (meaning, find q for lower tail.)
#' @return A numeric value representing the quantile for which the specified proportion of observations fall within.
#' @export
qBetaMS <- function(p, mean, var = NULL, sd = NULL, lt = TRUE) {
  if ((!is.null(var) & !is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  qbeta(p, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1, lower.tail = lt)
}

#' Random draw from the Standard Beta PDD with specific mean and variance.
#'
#' @description Draws random samples of observations from the Standard Beta probability density distribution with defined mean and variance.
#' @param n Number of observations to be drawn from under the Standard Beta PDD.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard probability density distribution.
#' @return A vector of length \code{n}, each value representing a random draw from the Standard Beta probability density distribution with defined mean and variance.
#' @export
rBetaMS <- function(n, mean, var = NULL, sd = NULL) {
  if ((!is.null(var) & !is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  rbeta(n, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1)
}

#' Coordinate generation for marking an area under the curve for the Standard Beta probability density distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Standard Beta probability density function.
#' @param from The point of the x-axis from where to start producing y-density values.
#' @param to The point of the x-axis to where y-density values are to be produced.
#' @param by The resolution (or specing) at which to produce y-density values.
#' @param alpha The Alpha shape-parameter value for the Standard Beta probability density distribution.
#' @param beta The Beta shape-parameter fort he Standard Beta probability density distribution.
#' @return A two-column matrix with density-values of y to plot against corresponding location values of x.
#' @export
Beta.gfx.poly.pdf <- function(from, to, by, alpha, beta) {
  x <- c(from, seq(from, to, by), to)
  for (i in 1:length(x)) {
    if (i == 1) y <- vector(length = length(x))
    if (i == 1 | i == length(x)) {
      y[i] <- 0
    } else {
      y[i] <- dbeta(x[i], alpha, beta)
    }
  }
  return(cbind(x, y))
}

#' Coordinate generation for marking an area under the curve for the Standard Beta cumulative probability density distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Standard Beta probability density function.
#' @param from The point of the x-axis from where to start producing y-density values.
#' @param to The point of the x-axis to where y-density values are to be produced.
#' @param by The resolution (or specing) at which to produce y-density values.
#' @param alpha The Alpha shape-parameter value for the Standard Beta probability density distribution.
#' @param beta The Beta shape-parameter fort he Standard Beta probability density distribution.
#' @return A two-column matrix with density-values of y to plot against corresponding location values of x.
#' @export
Beta.gfx.poly.cdf <- function(from, to, by, alpha, beta) {
  x <- c(from, seq(from, to, by), to)
  for (i in 1:length(x)) {
    if (i == 1) y <- vector(length = length(x))
    if (i == 1 | i == length(x)) {
      y[i] <- 0
    } else {
      y[i] <- pbeta(x[i], alpha, beta)
    }
  }
  return(cbind(x, y))
}

#' Most likely true alpha value given observed outcome.
#'
#' @description Given a fitted Standard Distribution, return the Alpha value where the observed mean becomes the mode.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @return The Alpha shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @export
MLA <- function(a, b, x = NULL, N = NULL) {
  if (is.null(x) | is.null(N)) {
    n <- a + b
    A <- (a*(n - 2) + n) / n
    A
  } else {
    x*(N - 2) + 1
  }
}

#' Most likely true beta value given observed outcome.
#'
#' @description Given a fitted Standard Beta Distribution, return the Beta value where the observed mean becomes the mode.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The Beta shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @export
MLB <- function(a, b, x = NULL, n = NULL) {
  if (is.null(x) | is.null(n)) {
    n <- a + b
    A <- (a*(n - 2) + n) / n
    B <- n - A
    B
  } else {
    n - (x*(n - 2) + 1)
  }
}

#' Most likely mean of the Standard Beta PDD, given that the observation is considered the most likely observation of the Standard Beta PDD (i.e., mode).
#'
#' @description Given a fitted Standard Beta Distribution, returns the expected mean of the distribution under the assumption that the observed value is the most likely value of the distribution.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param N Test-length.
#' @return The expected mean of the Standard Beta probability density distribution, for which the observed mean is the most likely value.
#' @export
MLM <- function(a, b, x = NULL, n = NULL) {
  if (is.null(x) | is.null(n)) {
    n <- a + b
    A <- (a*(n - 2) + n) / n
    B <- (n - A)
    A / (A + B)
  } else {
    (x*(n - 2) + 1) / n
  }
}

#' Probability Density under the Four-Parameter Beta PDD.
#'
#' @description Gives the density at desired values of X under the Four-Parameter Beta PDD.
#' @param x Value of X.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location paraeter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return The value for the probability density at specified values of X.
#' @export
dBeta.4P <- function(x, l, u, alpha, beta) {
  betafunc <- integrate(function(y) { y^(alpha - 1)*(1 - y)^(beta - 1) }, lower = 0, upper = 1)$value
  sapply(x, function(x) {
    if (x < l | x > u) {
      0
    } else {
      (1 / betafunc) * ((x - l)^(alpha - 1) * (u - x)^(beta - 1)) / (u - l)^(alpha + beta -1)
    }
    }
  )
}

#' Random Number Generation under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for generating random numbers from a specified four-parameter beta distribution.
#' @param n Number of draws.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return A vector with length \code{n} of random values drawn from the four-parameter beta distribution.
#' @export
rBeta.4P <- function(n, l, u, alpha, beta) {
  rbeta(n, alpha, beta) * (u - l) + l
}

#' Cumulative Probability Function under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under the four-parameter beta distribution.
#' @param q The quantile or a vector of quantiles for which the proportion is to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lt Whether the proportion to be calculated is to be under the lower or upper tail. Default is TRUE (lower tail).
#' @return A vector of proportions of observations falling under specified quantiles under the four-parameter beta distribution.
#' @export
pBeta.4P <- function(q, l, u, alpha, beta, lt = TRUE) {
  sapply(q, function(x) {
    num <- integrate(function(y) { dBeta.4P(y, l, u, alpha, beta) }, lower = l, upper = x)$value
    den <- integrate(function(y) { dBeta.4P(y, l, u, alpha, beta) }, lower = l, upper = u)$value
    if (lt) {
      num/den
      } else {
        1 - num/den
      }
    }
  )
}

#' Quantile given probability under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for calculating the quantile (i.e., value of x) for a given proportion (i.e., the value of y) under the four-parameter beta distribution.
#' @param p A vector (or single value) of proportions or probabilities for which the corresponding value of x (i.e., the quantiles) are to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lt Whether the quantile(s) to be calculated is to be under the lower or upper tail. Default is TRUE (lower tail).
#' @return A vector of quantiles for specified probabilities or proportions of observations under the four-parameter beta distribution.
#' @export
qBeta.4P <- function(p, l, u, alpha, beta, lt = TRUE) {
  if (lt) {
    qbeta(p, alpha, beta) * (u - l) + l
  } else {
    (1 - qbeta(p, alpha, beta)) * (u - l) + l
  }

}


#' Livingston and Lewis' "Effective Test Length"
#' @description  According to Livingston and Lewis (1995), "The effective test length corresponding to a test score is the number of discrete, dichotomously scored, locally independent, equally difficult items required to produce a total score of the same reliability."
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param l The lower-bound of the observed-score distribution.
#' @param u The upper-bound of the observed-score distribution.
#' @param reliability The reliability of the observed scores (correlation with true-score distribution).
#' @return An estimate of the effective length of a test, given the stability of the observations it produces.
ETL <- function(mean, variance, l = 0, u = 1, reliability) {
  ((mean - l) * (u - mean) - (reliability * variance)) / (variance * (1 - reliability))
}

k <- function(alpha, beta, rho) {
  K <- alpha + beta
  mu <- alpha / K
  sigma2 <- alpha * beta / (alpha + beta + 1) * (alpha + beta)^2
  error <- sigma2 * (1 - rho)
  K*((K - 1)*(sigma2 - error) - K*sigma2 + mu*(K - mu)) / 2*((mu*K - mu) - (sigma2 - error))
}

#' Beta Confusion Matrix.
#'
#' @description Function for calculating the confusion matrix for expected classification accuracy observed-scores, true-scores, and errors are distributed as Beta with shape-parameters alpha and beta.
#' @param x A vector of observed scores for which a beta-distribution is to be fitted.
#' @param alpha The first shape parameter of the observed-score distribution.
#' @param beta The second shape parameter of the observed-score distribution.
#' @param mean The mean of the observed-score distribution.
#' @param variance The variance of the observed-score distribution.
#' @param reliability The observed-score correlation with the true-score.
#' @param cut The cutoff value for classifying observations into pass or fail categories.
#' @return A confusion matrix estimating the proportion of true/false pass/fail categorizations for a test, given a specific distribution of observed scores.
#' @export

BCM <- function(x = NULL, alpha = NULL, beta = NULL, mean = NULL, variance = NULL, reliability, cut) {
  if (!is.null(x)) {
    alpha <- AMS(mean(x), var(x))
    beta <- BMS(mean(x), var(x))
  } else {
    if (!is.null(mean) & !is.null(variance)) {
      alpha <- AMS(mean, variance)
      beta <- BMS(mean, variance)
    }
  }
  N <- ETL(mean, variance, reliability = reliability)

  xaxis <- seq(.001, .999, .001)
  density <- dbeta(xaxis, alpha, beta) / sum(dbeta(xaxis, alpha, beta))

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
  return(list("EffectiveTestLength" = N, "ShapeParameters" = list("Alpha" = alpha, "Beta" = beta), "Confusionmatrix" = cmat))
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

#' Method of Moment Estimates of Shape- and Location Parameters of the Four-Parameter Beta Distribution.
#'
#' @description An implementation of the method of moments estimation of four-parameter beta distribution parameters presented by Hanson (1991). Given a string of values, calculates the shape- and location parameters required to produce a four-parameter beta distribution with the same mean, variance, skewness and kurtosis (i.e., the first four moments) as the observed-score distribution.
#' @param scores A vector of values to which the four-parameter beta distribution is to be fitted.
#' @return A list of parameter-values required to produce a four-parameter beta distribution with the same first four moments as the observed distribution.
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @export
Beta.4p.fit <- function(scores) {
  m <- observedmoments(scores)
  m1 <- m$raw[[1]]
  s2 <- m$central[[2]]
  g3 <- m$standardized[[3]]
  g4 <- m$standardized[[4]]
  r <- 6 * (g4 - g3^2 - 1) / (6 + 3 * g3^2 - 2 * g4)
  a <- r / 2 * (1 + sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
  b <- r / 2 * (1 - sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
  l <- m1 - ((a * sqrt(s2 * (a + b + 1))) / sqrt(a * b))
  u <- m1 + ((b * sqrt(s2 * (a + b + 1))) / sqrt(a * b))
  return(list("alpha" = a, "beta" = b, "l" = l, "u" = u))
}

