#' Compute Moments of Two-to-Four Parameter Beta Probability Density Distributions.
#'
#' @description Computes Raw, Central, or Standardized moment properties of defined Standard Beta probability density distributions.
#' @param a The Alpha shape parameter of the PDD.
#' @param b The Beta shape parameter of the PDD.
#' @param l The first (lower) location parameter of a four-parameter distribution.
#' @param u The second (upper) location parameter of a four-parameter distribution.
#' @param types A character vector determining which moment-types are to be calculated. Permissible values are "raw", "central", and "standardized".
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @examples
#' # Assume some variable follows a four-parameter beta distribution with
#' # location parameters l = 0.25 and u = .75, and shape
#' # parameters a = 5 and b = 3. To compute the first four
#' # raw, central, and standardized moments of this distrubution using
#' # betamoments():
#' betamoments(a = 5, b = 3, l = .25, u = .75,
#' types = c("raw", "central", "standardized"), orders = 4)
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
#' @param type A character vector determining which moment-types are to be calculated. Permissible values are \code{"raw"}, \code{"central"}, and \code{"standardized"}.
#' @param orders The number of moment-orders to be calculated for each of the moment-types.
#' @param correct Whether to include bias correction in estimation of orders. Default is \code{TRUE}.
#' @return A list of moment types, each a list of moment orders.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100))
#'
#' # To compute the first four raw, central, and standardized moments for this
#' # distribution of observed scores using observedmoments():
#' observedmoments(x = testdata, type = c("raw", "central", "standardized"),
#' orders = 4, correct = TRUE)
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
        sigma[i] <- base::sum((x - base::mean(x))^i)/(base::length(x) - 1)
      }
      else {
        sigma[i] <- base::sum((x - base::mean(x))^i)/(base::length(x))
      }
    }
    momentorders[[base::length(momentorders) + 1]] <- sigma
    base::names(momentorders)[types] <- "central"
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
#' @description Calculates the Alpha value required to produce a Standard  (two-parameter) Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Alpha shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0, rescaled to proportion
#' # of maximum.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 100))
#'
#' # To find the alpha shape-parameter of a Standard (two-parameter) Beta
#' # distribution with the same mean and variance as the observed-score
#' # distribution using AMS():
#' AMS(mean(testdata), var(testdata))
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
#' @description Calculates the Beta value required to produce a Standard (two-parameter) Beta probability density distribution with defined mean and variance or standard deviation.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @return A numeric value representing the required value for the Beta shape-parameter in order to produce a Standard Beta probability density distribution with the target mean and variance.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0, rescaled to proportion
#' # of maximum.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 100))
#'
#' # To find the beta shape-parameter of a Standard (two-parameter) Beta
#' # distribution with the same mean and variance as the observed-score
#' # distribution using BMS():
#' BMS(mean(testdata), var(testdata))
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
#' @examples
#' # To compute the proportion of the density under the lower-end tail of a
#' # point along the Standard (two-parameter) PDD (e.g., .5) with mean of .6
#' # and variance of .04:
#' pBetaMS(q = .5, mean = .6, var = .04)
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
#' @examples
#' # To compute the density at a specific point (e.g., .5) along the Standard
#' # (two-parameter) PDD with mean of .6 and variance of .04:
#' dBetaMS(x = .5, mean =.6, var = .04)
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
#' @description Calculates the quantile corresponding to a specific probability of some observation falling within the [0, x] (\code{lt = TRUE}) or [x, 1] (\code{lt = FALSE}) interval under the Standard Beta probability density distribution with defined mean and variance or standard deviation.
#' @param p A value of probability marking the point of the Y-axis to correspond to the X-axis.
#' @param mean The mean of the target Standard Beta probability density distribution.
#' @param var The variance of the target Standard Beta probability density distribution.
#' @param sd The standard deviation of the target Standard Beta probability density distribution.
#' @param lt Logical. Specifies which end of the tail for which to calculate quantile. Default is \code{TRUE} (meaning, find q for lower tail.)
#' @return A numeric value representing the quantile for which the specified proportion of observations fall within.
#' @examples
#' # To compute the quantile at a specific point (e.g., .5) along the Standard
#' # (two-parameter) PDD with mean of .6 and variance of .04:
#' qBetaMS(p = .5, mean =.6, var = .04)
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

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Standard Beta probability density function.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-density values.
#' @param to The point of the x-axis to where y-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce y-density values.
#' @param alpha The Alpha shape-parameter value for the Standard Beta probability density distribution.
#' @param beta The Beta shape-parameter fo rhe Standard Beta probability density distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with density-values of y to plot against corresponding location values of x.
#' @examples
#' # To box in an area under a four-parameter beta distribution with location
#' # parameters l = .25 and u = .75, and shape parameters
#' # alpha = 5 and beta = 3, from .4 to .6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 7))
#' coords <- Beta.gfx.poly.pdf(from = .4, to = .6, by = .001, alpha = 5,
#' beta = 3, l = .25, u = .75)
#' polygon(coords)
#' @export
Beta.gfx.poly.pdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  x <- base::c(from, base::seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- 0
    } else {
      y[i] <- dBeta.4P(x[i], l, u, alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Quantile Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Standard Beta probability quantile function.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-quantile values.
#' @param to The point of the \code{x}-axis to where \code{y}-quantile values are to be produced.
#' @param by The resolution (or spacing) at which to produce \code{y}-density values.
#' @param alpha The Alpha shape-parameter value for the Standard Beta probability distribution.
#' @param beta The Beta shape-parameter for the Standard Beta probability distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with quantile-values of \code{y} to plot against corresponding location values of \code{x}.
#' @examples
#' # To box in an area under a four-parameter beta quantile distribution with
#' # location parameters l = .25 and u = 75, and shape parameters
#' # alpha = 5 and beta = 3, from .4 to .6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1))
#' coords <- Beta.gfx.poly.qdf(from = .4, to = .6, by = .001, alpha = 5,
#' beta = 3, l = .25, u = .75)
#' polygon(coords)
#' @export
Beta.gfx.poly.qdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  x <- base::c(from, base::seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- l
    } else {
      y[i] <- qBeta.4P(x[i], l, u, alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Coordinate Generation for Marking an Area Under the Curve for the Beta Cumulative Probability Density Distribution.
#'
#' @description Plotting tool, producing a two-column matrix with values of \code{y} corresponding to locations on \code{x}. Useful for shading areas under the curve when tracing the line for the Standard Beta cumulative probability function.
#' @param from The point of the \code{x}-axis from where to start producing \code{y}-density values.
#' @param to The point of the \code{x}-axis to where \code{y}-density values are to be produced.
#' @param by The resolution (or spacing) at which to produce \code{y}-density values.
#' @param alpha The Alpha shape-parameter value for the Standard Beta cumulative probability distribution.
#' @param beta The Beta shape-parameter for the Standard Beta cumulative probability distribution.
#' @param l The lower-bound location parameter of the Beta distribution.
#' @param u The upper-bound location parameter of the Beta distribution.
#' @return A two-column matrix with cumulative probability-values of y to plot against corresponding location values of \code{x}.
#' @examples
#' # To box in an area under a four-parameter Beta cumulative distribution with
#' # location parameters l = .25 and u = 75, and shape parameters
#' # alpha = 5 and beta = 3, from .4 to .6:
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1))
#' coords <- Beta.gfx.poly.cdf(from = .4, to = .6, by = .001, alpha = 5,
#' beta = 3, l = .25, u = .75)
#' polygon(coords)
#' @export
Beta.gfx.poly.cdf <- function(from, to, by, alpha, beta, l = 0, u = 1) {
  x <- base::c(from, seq(from, to, by), to)
  for (i in 1:base::length(x)) {
    if (i == 1) y <- base::vector(length = base::length(x))
    if (i == 1 | i == base::length(x)) {
      y[i] <- 0
    } else {
      y[i] <- pBeta.4P(x[i], l, u, alpha, beta)
    }
  }
  return(base::cbind(x, y))
}

#' Most Likely True Alpha Value Given Observed Outcome.
#'
#' @description Given a fitted Standard (two-parameter) Beta Distribution, return the alpha shape-parameter value where the observed mean becomes the mode.
#' @param a Observed alpha-parameter value for fitted Standard Beta PDD.
#' @param b Observed beta-parameter value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @return The Alpha shape-parameter value for the Standard Beta probability density distribution where the observed mean is the expected mode.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-alpha parameter most likely to have produced the observations:
#' MLA(a = 10, b = 8)
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
#' @description Assuming a prior standard (two-parameter) Beta Distribution, return the beta shape-parameter value where the observed mean becomes the mode.
#' @param a Observed alpha-parameter value for fitted Standard Beta PDD.
#' @param b Observed beta-parameter value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-beta parameter most likely to have produced the observations:
#' MLB(a = 10, b = 8)
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
#' @description Assuming a prior Standard (two-parameter) Beta Distribution, returns the expected mean of the distribution under the assumption that the observed value is the most likely value of the distribution.
#' @param a Observed alpha value for fitted Standard Beta PDD.
#' @param b Observed beta value for fitted Standard Beta PDD.
#' @param x Observed proportion-correct outcome.
#' @param n Test-length.
#' @examples
#' # Assuming a prior Standard (two-parameter) Beta distribution is fit, which
#' # yield an alpha parameter of 10 and a beta parameter of 8, calculate the
#' # true-mean most likely to have produced the observations:
#' MLM(a = 10, b = 8)
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
#' @description Gives the density at desired values of \code{x} under the Four-Parameter Beta PDD.
#' @param x Value of \code{x}.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return The value for the probability density at specified values of \code{x}.
#' @examples
#' # Assume some variable follows a four-parameter beta distribution with
#' # location parameters l = 0.25 and u = .75, and shape
#' # parameters alpha = 5 and beta = 3. To compute the
#' # probability density at a specific point of the distribution (e.g., .5)
#' # using dBeta.4P():
#' dBeta.4P(x = .5, l = .25, u = .75, alpha = 5, beta = 3)
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
#' @description Function for generating random numbers from a specified Four-Parameter Beta Distribution.
#' @param n Number of draws.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @return A vector with length \code{n} of random values drawn from the Four-Parameter Beta Distribution.
#' @examples
#' # Assume some variable follows a four-parameter beta distribution with
#' # location parameters l = 0.25 and u = .75, and shape
#' # parameters alpha = 5 and beta = 3. To draw a random
#' # value from this distribution using rBeta.4P():
#' rBeta.4P(n = 1, l = .25, u = .75, alpha = 5, beta = 3)
#' @export
rBeta.4P <- function(n, l, u, alpha, beta) {
  stats::rbeta(n, alpha, beta) * (u - l) + l
}

#' Cumulative Probability Function under the Four-Parameter Beta Probability Density Distribution.
#'
#' @description Function for calculating the proportion of observations up to a specifiable quantile under the Four-Parameter Beta Distribution.
#' @param q The quantile or a vector of quantiles for which the proportion is to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lt Whether the proportion to be calculated is to be under the lower or upper tail. Default is \code{TRUE} (lower tail).
#' @return A vector of proportions of observations falling under specified quantiles under the four-parameter beta distribution.
#' @examples
#' # Assume some variable follows a four-parameter beta distribution with
#' # location parameters l = 0.25 and u = .75, and shape
#' # parameters alpha = 5 and beta = 3. To compute the
#' # cumulative probability at a specific point of the distribution (e.g., .5)
#' # using pBeta.4P():
#' pBeta.4P(q = .5, l = .25, u = .75, alpha = 5, beta = 3)
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

#' Quantile Given Probability Under the Four-Parameter Beta Distribution.
#'
#' @description Function for calculating the quantile (i.e., value of \code{x}) for a given proportion (i.e., the value of \code{y}) under the Four-Parameter Beta Distribution.
#' @param p A vector (or single value) of proportions or probabilities for which the corresponding value of \code{x} (i.e., the quantiles) are to be calculated.
#' @param l The first (lower) location parameter.
#' @param u The second (upper) location parameter.
#' @param alpha The first shape parameter.
#' @param beta The second shape parameter.
#' @param lt Logical. Whether the quantile(s) to be calculated is to be under the lower or upper tail. Default is \code{TRUE} (lower tail).
#' @return A vector of quantiles for specified probabilities or proportions of observations under the four-parameter beta distribution.
#' @examples
#' # Assume some variable follows a four-parameter beta distribution with
#' # location parameters l = 0.25 and u = .75, and shape
#' # parameters alpha = 5 and beta = 3. To compute the
#' # quantile at a specific point of the distribution (e.g., .5)
#' # using qBeta.4P():
#' qBeta.4P(p = .5, l = .25, u = .75, alpha = 5, beta = 3)
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
#' @description An implementation of the method of moments estimation of four-parameter beta distribution parameters presented by Hanson (1991). Given a vector of values, calculates the shape- and location parameters required to produce a four-parameter beta distribution with the same mean, variance, skewness and kurtosis (i.e., the first four moments) as the observed-score distribution.
#' @param scores A vector of values to which the four-parameter beta distribution is to be fitted.
#' @return A list of parameter-values required to produce a four-parameter beta distribution with the same first four moments as the observed distribution.
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3))
#' hist(testdata, xlim = c(0, 100), freq = FALSE)
#'
#' # To fit and retrieve the parameters for a four-parameter beta distribution
#' # to the observed-score distribution using Beta.4p.fit():
#' (params.4p <- Beta.4p.fit(testdata))
#' curve(dBeta.4P(x, params.4p$l, params.4p$u, params.4p$alpha, params.4p$beta), add = TRUE)
#' @export
Beta.4p.fit <- function(scores) {
  m <- observedmoments(scores)
  m1 <- m$raw[[1]]
  s2 <- m$central[[2]]
  g3 <- m$standardized[[3]]
  g4 <- m$standardized[[4]]
  r <- 6 * (g4 - g3^2 - 1) / (6 + 3 * g3^2 - 2 * g4)
  if (g3 < 0) {
    a <- r / 2 * (1 + base::sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
    b <- r / 2 * (1 - base::sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
  } else {
    b <- r / 2 * (1 + sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
    a <- r / 2 * (1 - sqrt(1 - ((24 * (r + 1)) / ((r + 2) * (r + 3) * g4 - 3 * (r - 6) * (r + 1)))))
  }
  l <- m1 - ((a * base::sqrt(s2 * (a + b + 1))) / base::sqrt(a * b))
  u <- m1 + ((b * base::sqrt(s2 * (a + b + 1))) / base::sqrt(a * b))
  return(base::list("alpha" = a, "beta" = b, "l" = l, "u" = u))
}

#' Method of Moment Estimates of Shape-Parameters of the Two-Parameter (Standard) Beta Distribution.
#'
#' @description An implementation of the method of moments estimation of two-parameter beta distribution parameters. Given a vector of values, calculates the shape parameters required to produce a two-parameter beta distribution with the same mean and variance (i.e., the first two moments) as the observed-score distribution.
#' @param scores A vector of values to which the two-parameter beta distribution is to be fitted. The values ought to fall within the [0, 1] interval.
#' @return A list of parameter-values required to produce a Standard two-parameter beta distribution with the same first two moments as the observed distribution.
#' @examples
#' # Generate some fictional data. Say, 100 individuals take a test with a
#' # maximum score of 100 and a minimum score of 0.
#' set.seed(1234)
#' testdata <- rbinom(100, 100, rBeta.4P(100, .25, .75, 5, 3)) / 100
#' hist(testdata, xlim = c(0, 1), freq = FALSE)
#'
#' # To fit and retrieve the parameters for a two-parameter beta distribution
#' # to the observed-score distribution using Beta.2p.fit():
#' (params.2p <- Beta.2p.fit(testdata))
#' curve(dbeta(x, params.2p$alpha, params.2p$beta), add = TRUE)
#' @export
Beta.2p.fit <- function(scores) {
  if (max(scores) > 1 | min(scores) < 0) {
    warning("Input values outside the range of the Standard Beta distribution (i.e., there are values falling outside the [0, 1] interval).")
  }
  m <- observedmoments(scores)
  m1 <- m$raw[[1]]
  s2 <- m$central[[2]]
  a <- AMS(m1, s2)
  b <- BMS(m1, s2)
  return(base::list("alpha" = a, "beta" = b, "l" = 0, "u" = 1))
}

#' An implementation of the Beta-density Compound Cumulative-Binomial Distribution.
#'
#' @description The Beta Compound Binomial distribution: The product of the four-parameter Beta probability density function and the binomial cumulative probability mass function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param a The alpha shape-parameter of the Beta distribution.
#' @param b The beta shape-parameter of the Beta distribution.
#' @param n The number of trials for the Binomial distribution.
#' @param c The "true-cut" (proportion) of on the Binomial distribution.
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Binomial distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Binomial error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0 can be calculated as:
#' dBeta.pBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pBinom(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5, lower.tail = TRUE)
#'
#' #By integration, the population proportion of (e.g.) passing scores in some
#' #region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' #calculated as:
#' integrate(function(x) { dBeta.pBinom(x, 0.25, .75, 5, 3, 10, 0.5) }, lower = 0.25, upper = 0.5)
#' @export
dBeta.pBinom <- function(x, l, u, a, b, n, c, lower.tail = FALSE) {
  if (!lower.tail) {
    dBeta.4P(x, l, u, a, b) * stats::pbinom(n * c, n, x, lower.tail = lower.tail)
  } else {
    dBeta.4P(x, l, u, a, b) * (1 - stats::pbinom(n * c, n, x, lower.tail = lower.tail))
  }
}

#' An implementation of the Beta-density Compound Cumulative-Beta Distribution.
#'
#' @description The Beta Compound Beta distribution: The product of the four-parameter Beta probability density function and the beta cumulative probability function. Used in the Livingston and Lewis approach to classification accuracy and consistency, the output can be interpreted as the population density of passing scores produced at "x" (a value of true-score).
#' @param x x-axis input for which \code{p} (proportion or probability) is to be computed.
#' @param l The lower-bound of the four-parameter Beta distribution.
#' @param u The upper-bound of the four-parameter Beta distribution.
#' @param a The alpha shape-parameter of the Beta density distribution.
#' @param b The beta shape-parameter of the Beta density distribution.
#' @param n The number of trials for the Beta cumulative probability distribution.
#' @param c The "true-cut" (proportion) of on the Beta cumulative probability distribution.
#' @param lower.tail Logical. Whether to compute the lower or upper tail of the Beta cumulative probability distribution. Default is \code{FALSE} (i.e., upper tail).
#' @references Hanson, Bradley A. (1991). Method of Moments Estimates for the Four-Parameter Beta Compound Binomial Model and the Calculation of Classification Consistency Indexes.American College Testing Research Report Series.
#' @references Livingston, Samuel A. and Lewis, Charles. (1995). Estimating the Consistency and Accuracy of Classifications Based on Test Scores. Journal of Educational Measurement, 32(2).
#' @references Lord, Frederic M. (1965). A Strong True-Score Theory, With Applications. Psychometrika, 30(3).
#' @examples
#' # Given a four-parameter Beta distribution with parameters l = 0.25, u = 0.75,
#' # alpha = 5, and beta = 3, and a Beta error distribution with number of
#' # trials (n) = 10 and a cutoff-point (c) at 50% correct (i.e., proportion correct
#' # of 0.5), the population density of passing scores produced at true-score
#' # (x) = 0.5 can be calculated as:
#' dBeta.pBeta(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5)
#'
#' # Conversely, the density of failing scores produced at x can be calculated
#' # by passing the additional argument "lower.tail = TRUE" to the function.
#' # That is:
#' dBeta.pBeta(x = 0.5, l = 0.25, u = 0.75, a = 5, b = 3, n = 10, c = 0.5, lower.tail = TRUE)
#'
#' # By integration, the population proportion of (e.g.) passing scores in some
#' # region of the true-score distribution (e.g. between 0.25 and 0.5) can be
#' # calculated as:
#' integrate(function(x) { dBeta.pBeta(x, 0.25, .75, 5, 3, 10, 0.5) }, lower = 0.25, upper = 0.5)
#' @export
dBeta.pBeta <- function(x, l, u, a, b, n, c, lower.tail = FALSE) {
  if(!lower.tail) {
    dBeta.4P(x, l, u, a, b) * stats::pbeta(c, x * n, (1 - x) * n, lower.tail = FALSE)
  } else {
    dBeta.4P(x, l, u, a, b) * (1 - stats::pbeta(c, x * n, (1 - x) * n, lower.tail = FALSE))
  }
}
