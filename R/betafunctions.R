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
betamoments <- function(a, b, l = 0, u = 1, types = c("raw", "central", "standardized"), orders = 4) {
  BETAMOMENTS <- base::rep(base::list(base::rep(base::list(NULL), orders)), base::length(types))
  TYPE <- 1
  if (any(types == "raw")) {
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- stats::integrate(function(x) { dBeta.4P(x, l, u, a, b) * x^i },
                                            lower = l, upper = u)$value
    }
    base::names(BETAMOMENTS)[TYPE] <- "raw"
    TYPE <- TYPE + 1
  }
  if (any(types == "central")) {
    Mu <- stats::integrate(function(x) { dBeta.4P(x, l, u, a, b) * x^1 }, lower = l, upper = u)$value
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- stats::integrate(function(x) { dBeta.4P(x, l, u, a, b) * (x - Mu)^i },
                                            lower = l, upper = u)$value
    }
    base::names(BETAMOMENTS)[TYPE] <- "central"
    TYPE <- TYPE + 1
  }
  if (base::any(types == "standardized")) {
    Mu <- stats::integrate(function(x) { dBeta.4P(x, l, u, a, b) * x^1 }, lower = l, upper = u)$value
    SigmaSquared <- stats::integrate(function(x) { dBeta.4P(x, l, u, a, b) * (x - Mu)^2 },
                       lower = l, upper = u)$value
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- stats::integrate(function(x) { dBeta.4P(x, l, u, a, b) * ((x - Mu)^i / sqrt(SigmaSquared)^i) },
                                            lower = l, upper = u)$value
    }
    base::names(BETAMOMENTS)[TYPE] <- "standardized"
  }
  return(BETAMOMENTS)
}

#' Compute Moments of Observed Value Distribution.
#'
#' @description Computes Raw, Central, or Standardized moment properties of a vector of observed scores.
#' @param x A vector of values, the distribution of which moments are to be calculated.
#' @param type A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @param correct Whether to include bias correction in estimation of orders. Default is TRUE.
#' @return A list of moment types, each a list of moment orders.
#' @export
observedmoments <- function(x, type = c("raw", "central", "standardized"),  orders = 4, correct = TRUE) {
  x <- stats::na.omit(x)
  types <- 1
  momentorders <- base::list()
  if (base::any(type == "raw")) {
    mu <- list(rep(vector(length = 1), orders))
    for (i in 1:orders) {
      mu[i] <- base::sum(x^i)/base::length(x)
    }
    momentorders[[base::length(momentorders) + 1]] <- mu
    base::names(momentorders)[types] <- "raw"
    types <- types + 1
  }
  if (base::any(type == "central")) {
    sigma <- base::list(base::rep(base::vector(length = 1), orders))
    for (i in 1:orders) {
      if (correct) {
        sigma[i] <- stats::var(x)
      }
      else {
        sigma[i] <- base::sum((x - base::mean(x))^i)/(base::length(x))
      }
    }
    momentorders[[length(momentorders) + 1]] <- sigma
    names(momentorders)[types] <- "central"
    types <- types + 1
  }
  if (base::any(type == "standardized")) {
    gamma <- base::list(base::rep(base::vector(length = 1), orders))
    for (i in 1:orders) {
      if (correct) {
        gamma[i] <- 1 / base::length(x) * base::sum(((x - base::mean(x))^i)/sqrt(stats::var(x))^i)
      }
      else {
        gamma[i] <- 1 / base::length(x) * base::sum(((x - base::mean(x))^i)/sqrt(stats::var(x))^i)
      }
    }
    momentorders[[base::length(momentorders) + 1]] <- gamma
    names(momentorders)[types] <- "standardized"
  }
  return(momentorders)
}

#' Alpha Shape Parameter Given Mean and Variance of a Standard Beta PDD.
#'
#' @description Calculates the Alpha value required to produce a Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Alpha shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
AMS <- function(mean, var, sd = NULL) {
  if ((!base::is.null(var) & !base::is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (base::is.null(var) & !base::is.null(sd)) var <- sd^2
  if(((mean^2 - mean^3) / var) - mean <= 0) {
    warning("Parameter out of bounds (Alpha <= 0).")
  }
  return(((mean^2 - mean^3) / var) - mean)
}

#' Beta Shape Parameter Given Mean and Variance of a Standard Beta PDD.
#'
#' @description Calculates the Beta value required to produce a Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
BMS <- function(mean, var, sd = NULL) {
  if ((!base::is.null(var) & !base::is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (base::is.null(var) & !base::is.null(sd)) var <- sd^2
  if((mean * (1 - mean)^2) / var + mean - 1 <= 0) {
    warning("Parameter out of bounds (Beta <= 0).")
  }
  return(((mean * (1 - mean)^2) / var) + mean - 1)
}

#' Probability of Some Specific Observation under the Standard Beta PDD with Specific Mean and Variance.
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
  if (base::is.null(var) & !base::is.null(sd)) var <- sd^2
  stats::pbeta(q, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1, lower.tail = lt)
}

#' Density Under a Specific Point of the Standard Beta PDD with Specific Mean and Variance or Standard Deviation.
#'
#' @description Calculates the density under specific points of the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param x A specific point on the x-axis of the Standard Beta PDD.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta Shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
dBetaMS <- function(x, mean, var = NULL, sd = NULL) {
  if ((!base::is.null(var) & !base::is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (base::is.null(var) & !base::is.null(sd)) var <- sd^2
  stats::dbeta(x, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1)
}

#' Quantile Containing Specific Proportion of the Distribution, Given a Specific Probability of the Standard Beta PDD with Specific Mean and Variance or Standard Deviation.
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
  if ((!base::is.null(var) & !base::is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (base::is.null(var) & !base::is.null(sd)) var <- sd^2
  stats::qbeta(p, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1, lower.tail = lt)
}

#' Random Draw from the Standard Beta PDD With Specific Mean and Variance.
#'
#' @description Draws random samples of observations from the Standard Beta probability density distribution with defined mean and variance.
#' @param n Number of observations to be drawn from under the Standard Beta PDD.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard probability density distribution.
#' @return A vector of length \code{n}, each value representing a random draw from the Standard Beta probability density distribution with defined mean and variance.
#' @export
rBetaMS <- function(n, mean, var = NULL, sd = NULL) {
  if ((!base::is.null(var) & !base::is.null(sd))) {
    if (var != sd^2) {
      warning("Nonequivalent values of VAR and SD specified. Using VAR.")
    }
  }
  if (is.null(var) & !is.null(sd)) var <- sd^2
  stats::rbeta(n, ((mean^2 - mean^3) / var) - mean, (mean * (1 - mean)^2) / var + mean - 1)
}

#' Coordinate Generation for Marking an Area Under the Curve for the Standard Beta Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Standard Beta probability density function.
#' @param from The point of the x-axis from where to start producing y-density values.
#' @param to The point of the x-axis to where y-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce y-density values.
#' @param alpha The Alpha shape-parameter value for the Standard Beta probability density distribution.
#' @param beta The Beta shape-parameter fort he Standard Beta probability density distribution.
#' @return A two-column matrix with density-values of y to plot against corresponding location values of x.
#' @export
Beta.gfx.poly.pdf <- function(from, to, by, alpha, beta) {
  x <- base::c(from, base::seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- 0
    } else {
      y[i] <- stats::dbeta(x[i], alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Standard Beta Cumulative Probability Density Distribution.
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
  x <- base::c(from, seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- 0
    } else {
      y[i] <- stats::pbeta(x[i], alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Most Likely True Alpha Value Given Observed Outcome.
#'
#' @description Given a fitted Standard Distribution, return the Alpha value where the observed mean becomes the mode.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The Alpha shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @export
MLA <- function(a, b, x = NULL, n = NULL) {
  if (base::is.null(x) | base::is.null(n)) {
    n <- a + b
    A <- (a*(n - 2) + n) / n
    A
  } else {
    x*(n - 2) + 1
  }
}

#' Most Likely True Beta Value Given Observed Outcome.
#'
#' @description Assuming a prior standard Beta distribution, return the Beta value where the observed mean becomes the mode.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The Beta shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @export
MLB <- function(a, b, x = NULL, n = NULL) {
  if (base::is.null(x) | base::is.null(n)) {
    n <- a + b
    A <- (a*(n - 2) + n) / n
    B <- n - A
    B
  } else {
    n - (x*(n - 2) + 1)
  }
}

#' Most Likely Mean of the Standard Beta PDD, Given that the Observation is Considered the Most Likely Observation of the Standard Beta PDD (i.e., Mode).
#'
#' @description Assuming a prior standard Beta distribution, returns the expected mean of the distribution under the assumption that the observed value is the most likely value of the distribution.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The expected mean of the Standard Beta probability density distribution, for which the observed mean is the most likely value.
#' @export
MLM <- function(a, b, x = NULL, n = NULL) {
  if (base::is.null(x) | base::is.null(n)) {
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
 bfunc <- base::beta(alpha, beta)
  base::sapply(x, function(x) {
    if (x < l | x > u) {
      0
    } else {
      (1 / bfunc) * ((x - l)^(alpha - 1) * (u - x)^(beta - 1)) / (u - l)^(alpha + beta -1)
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
  stats::rbeta(n, alpha, beta) * (u - l) + l
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
    num <- stats::integrate(function(y) { dBeta.4P(y, l, u, alpha, beta) }, lower = l, upper = x)$value
    den <- stats::integrate(function(y) { dBeta.4P(y, l, u, alpha, beta) }, lower = l, upper = u)$value
    if (lt) {
      num/den
      } else {
        1 - num/den
      }
    }
  )
}

#' Quantile Given Probability Under the Four-Parameter Beta Probability Density Distribution.
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
    stats::qbeta(p, alpha, beta) * (u - l) + l
  } else {
    (1 - stats::qbeta(p, alpha, beta)) * (u - l) + l
  }
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

