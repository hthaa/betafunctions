#' Compute Moments of Beta Probability Density Distributions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Standard Beta probability density distributions.
#' @param a The Alpha shape parameter of the PDD.
#' @param b The Beta shape parameter of the PDD.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @return A list of moment types, each a list of moment orders.
#' @export
betamoments <- function(a, b, mean = NULL, var = NULL, sd = NULL, types = c("raw", "central", "standardized"), orders = 4) {
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
      BETAMOMENTS[[TYPE]][[i]] <- integrate(function(x) { dbeta(x, a, b) * x^i }, lower = 0, upper = 1)$value
    }
    names(BETAMOMENTS)[TYPE] <- "raw"
    TYPE <- TYPE + 1
  }
  if (any(types == "central")) {
    Mu <- integrate(function(x) { dbeta(x, a, b) * x }, lower = 0, upper = 1)$value
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- integrate(function(x) { dbeta(x, a, b) * (x - Mu)^i },
                                            lower = 0, upper = 1)$value
    }
    names(BETAMOMENTS)[TYPE] <- "central"
    TYPE <- TYPE + 1
  }
  if (any(types == "standardized")) {
    Mu <- integrate(function(x) { dbeta(x, a, b) * x }, lower = 0, upper = 1)$value
    Sigma <- integrate(function(x) { dbeta(x, a, b) * (x - Mu)^2 }, lower = 0, upper = 1)$value
    for (i in 1:orders) {
      BETAMOMENTS[[TYPE]][[i]] <- integrate(function(x) { dbeta(x, a, b) * ((x - Mu)^i / sqrt(Sigma)^i) },
                                            lower = 0, upper = 1)$value
    }
    names(BETAMOMENTS)[TYPE] <- "standardized"
  }
  return(BETAMOMENTS)
}

#' Alpha shape parameter given mean and variance of a Standard Beta PDD.
#'
#' @description Calculates the Alpha value required to produce a Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Alpha shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @export
AMS <- function(mean, var, sd) {
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
BMS <- function(mean, var, sd) {
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
