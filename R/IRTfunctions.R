#' Heavy-Handed Maximum Likelihood Estimation of Theta.
#' @description Estimates Theta given supplied response-pattern and item properties by evaluating the log-likelihood of Theta across the Theta spectrum and choosing the value of Theta corresponding to the largest value.
#' @param rp A vector representing a specific response-pattern.
#' @param ip A table of item properties (columns representing difficulty, discrimination, and guessing parameters, respectively).
#' @param ThetaRange The range across which Theta values are to be evaluated with respect to their probability of producing the response pattern.
#' @param resoluation The resolution at which Theta is to be evaluated (i.e., desired decimal-level of precision).
#' @returns A matrix-row with the first column representing the Theta estimate, the second its standard error, and the third the natural logarithm of its likelihood.
HHMLE <- function(rp, ip, ThetaRange = c(-6, 6), resolution = .001) {
  Theta <- seq(ThetaRange[1], ThetaRange[2], resolution)
  lnL.Theta <- sapply(Theta, function(x) { LNLRV(rp, ip, x) })
  Theta_hat <- Theta[which(lnL.Theta == max(lnL.Theta))]
  Theta_hat.SE <- 1 / sqrt(sum(PCR(ip, Theta_hat) * (1 - PCR(ip, Theta_hat))))
  lnL.Theta_hat <- max(lnL.Theta)
  output <- matrix(c(Theta_hat, Theta_hat.SE, lnL.Theta_hat), ncol = 3)
  colnames(output) <- c("Theta_hat", "Theta_hat.SE", "lnL.Theta_hat")
  return(output)
}

#' Newton's Method of Maximum Likelihood Estimation applied to Theta.
#' @description An implementation of Newton's method for MLE of the person ability IRT parameter (aka. Theta). Finds the value of Theta which maximizes the probability of the observed response-pattern occuring, by substracting the derivative of the estimate from the estimate until successive improvements in the Theta estimates are inconsequentially small.
#' @param rv The response-vector for which Theta is to be estimated.
#' @param ips A table of item properties, where collumns represent the difficulty, discrimination, and guessing parameters (respectively).
#' @param theta_hat the initial guestimate of Theta.
#' @param ccrit The convergece criterion of Theta (i.e., desired degree of accuracy for the Theta estimate).
#' @return An estimate for Theta which maximizes the likelihood of the observations.
#' @export
NMLE <- function(rv, ips, Theta = 0, ccrit = .0001, max.iter = 2000) {
  lnL.history <- vector("numeric")
  Theta.history <- vector("numeric")
  Step.history <- vector("numeric")
  for (i in 1:max.iter) {
    lnL.history[i] <- LNLRV(rv, ips, Theta)
    Theta.history[i] <- Theta
    dy <- (LNLRV(rv, ips, Theta) - LNLRV(rv, ips, Theta + .0001)) / .0001
    Theta2 <- Theta - dy
    Step.history[i] <- Theta - Theta2
    if (abs(Theta - Theta2) < ccrit) {
      break
    }
    Theta <- Theta2
  }
  if (i == max.iter) warning("Solution did not converge.")
  history <- list("lnL.history" = lnL.history, "Theta.history" = Theta.history, "Step.history" = Step.history, "Theta_hat" = Theta)
  return(history)
}

#' The Natural-Log of the Likelihood of a Response Vector given a value of Theta.
#' @description Calculates the natoral logarithm of the likelihood of an individual with ability-value Theta attaining a specific response vector given a set of item-parameters, assuming local independence.
#' @param x A vector of dichotomous values (i.e., 0/1) with length equal to the number of items.
#' @param iptable A table of item-parameters, where the first collumn represents the difficulty parameters, the second the discrimination parameters, and the thirt the guessing parameters.
#' @param theta A specific value of Theta ("Ability") for which the probability of the supplied response vector is to be calculated.
#' @return The log of the likelihood of producing the supplied response vector at a given value of Theta.
#' @export
LNLRV <- function(x, iptable, theta) {
  for (i in 1:length(x)) {
    if (i == 1) p.Theta <- vector(length = length(x))
    if (x[i] == 0) p.Theta[i] <- log(1 - PCR(matrix(iptable[i, ], ncol = 3), theta))
    if (x[i] == 1) p.Theta[i] <- log(PCR(matrix(iptable[i, ], ncol = 3), theta))
  }
  return(sum(p.Theta))
}

#' Probability of Correct Response.
#' @description Calculates the probability of providing a response of "1" to a particular item with particular properties under 1, 2, and 3 parameter logistic IRT models, given a specific value of Theta.
#' @param iptable A table of item properties, where column 1 represents the difficulty parameters, column 2 the discrimination parameters, and column 3 the guessing parameters of the items.
#' @param theta The value of Theta for which the probabilities for responding "1" is to be calculated.
#' @return The probability for responding "1" on one or more items.
#' @export
PCR <- function(iptable, theta = 0) {
  for(i in 1:nrow(iptable)) {
    if (i == 1) pcrv <- vector(length = nrow(iptable))
    pcrv[i] <- iptable[i, 3] + (1 - iptable[i, 3]) * (1 / (1 + (exp(-iptable[i, 2] * (theta - iptable[i, 1])))))
  }
  return(pcrv)
}

#' Item Properties Generator.
#' @description Generates a matrix of dichotomous-item IRT-parameters (discrimination, difficulty, and guessing), where the first column represents the difficulty parameter, the second the discrimination parameter, and the thirt the guessing parameter for the item.
#' @param nitem The number of items for which item parameters are to be generated.
#' @param b.range The range to be covered by the difficulty parameters.
#' @param rndm.b Logical. Whether item difficulty parameters are to be uniformly randomly distributed in the specified range, or whether items are to cover the range with equal spacing between each item. Defaults to FALSE (equal spacing).
#' @param a.range The range to be covered by the discrimination parameters.
#' @param rand.a Logical. Whether item discrimination parameters are to be uniformly randomly distributed in the specified range, or whether items are to cover the range with equal spacing. Defaults to TRUE (uniform random distribution).
#' @param c.range The range to be covered by the discrimination parameters.
#' @param rand.c Logical. Whether item guessing parameters are to be uniformly randomly distributed in the specified range, or whether items are to cover the range with equal spacing. Defaults to FALSE (uniform random distribution).
#' @return A matrix of item properties, where [, 1] represents difficulty, [, 2] discrimination, and [, 3] guessing.
#' @export
IPG <- function(nitem = 10, b.range = c(-3, 3), rndm.b = FALSE, a.range = c(1, 1), rand.a = TRUE,  c.range = c(0, 0), rand.c = TRUE) {
  item.props <- cbind(
    cbind(
      b = if (b.range[1] - b.range[2] != 0) {
        if (rndm.b == FALSE) {
          seq(b.range[1], b.range[2], (b.range[2] - b.range[1]) / (nitem - 1))
        } else {
          runif (nitem, b.range[1], b.range[2])
        }
      } else {
        rep(b.range[1], nitem)
      }
      ,
      a = if (a.range[1] - a.range[2] != 0) {
        if (rndm.a == FALSE) {
          seq(a.range[1], a.range[2], (a.range[2] - a.range[1]) / (nitem - 1))
        } else {
          runif (nitem, a.range[1], a.range[2])
        }
      } else {
        rep(a.range[1], nitem)
      }
    )
    ,
    c = if (c.range[1] - c.range[2] != 0) {
      if (rndm.c == FALSE) {
        seq(c.range[1], c.range[2], (c.range[2] - c.range[1]) / (nitem - 1))
      } else {
        runif (nitem, c.range[1], c.range[2])
      }
    } else {
      rep(c.range[1], nitem)
    }
  )
  return(item.props)
}

#' Dichotomous Data Generator.
#' @description generated a data-frame where ncol = number of items and nrow = number of respondents. Responses take the form of dichotomous values (0/1), where the values of each individual observation is probabilistically generated based on the respondents probability of responding "1" to the specific item with the properties "difficulty" ("b"), "discrimination" ("a"), and "guessing" ("c").
#' @param iptable A table where the number of rows equals the number of items, and the columns represent the items' difficulty ([, 1]), discrimination ([, 2]), and guessing ([, 3]) parameters, respectively.
#' @param theta A vector of values representing respondents standing on the latent variable, where vector length is equal to the number of respondents responding to the items.
#' @return A nRespondent x nItem data-frame of dichotomous responses, probabilistically generated in accordance with the respondents' probability of correct response under some previously specified IRT model.
#' @export
DDG <- function(ipt, theta) {
  nitem <- nrow(ipt)
  nresp <- length(theta)
  dichotomous.data <- matrix(data = runif(nitem*nresp, 0, 1), nrow = nresp, ncol = nitem)
  for (i in 1:nitem) {
    pitemr <- sapply(dichotomous.data[, i], function(x) { PCR(ipt, theta) })
    itemresponse <- ifelse(runif(nresp, 0, 1) > pitemr, 0, 1)
    dichotomous.data[, i] <- itemresponse
  }
  return(dichotomous.data)
}
