ffunc <- function(x, ahat1, bhat1, alphahat1, betahat1, lowbinom, highbinom, n1) {
  (dbeta(((x - ahat1)) / ( bhat1 - ahat1), alphahat1, betahat1) * (pbinom(highbinom, n1, x) - pbinom(lowbinom, n1, x))) / (bhat1 - ahat1)
}

llclassify <- function(observedy, r, maxval, minval = 0, cutscore, adjust = 1) {
  observed <- observedy - minval
  maxval <- maxval - minval
  cutscore <- cutscore - minval
  x <- round(observed)
  p <- x / maxval
  mup <- sum(p) / length(p)
  deviations <- (p - mup)
  squaredeviations <- deviations^2
  sigma2p <- sum(squaredeviations) / length(p)

  effective <- (mup * (1 - mup) - r * sigma2p) / (sigma2p * (1 - r))
  n <- round(effective)

  xprimenotrounded <- (n * x) / maxval
  xprime <- round(xprimenotrounded)
  mprime1x <- sum(xprime) / length(xprime)
  mprime2x <- sum(xprime^2) / length(xprime)
  mprime3x <- sum(xprime^3) / length(xprime)
  mprime4x <- sum(xprime^4) / length(xprime)
  mprime1p <- mprime1x/n
  mprime2p <- (mprime2x - mprime1x) / (n * (n - 1))
  mprime3p <- (mprime3x - 3 * mprime2x + 2 * mprime1x) / (n * (n - 1) * (n - 2))
  mprime4p <- (mprime4x - 6 * mprime3x + 11 * mprime2x - 6 * mprime1x) / (n * (n - 1) * (n - 2) * (n - 3))
  m2p <- mprime2p - mprime1p^2
  m3p <- mprime3p - 3 * mprime1p * mprime2p + 2 * mprime1p^3
  m4p <- mprime4p - 4 * mprime1p * mprime3p + 6 * mprime1p^2 * mprime2p - 3 * mprime1p^4
  k <- 16 * m2p^3 / m3p^2
  l <- m4p * m2p / m3p^2
  phi <- (3 * (k - 16 * (l - 1))) / (16 * (l - 1) - 8 - 3 * k)
  sgnm3p <- m3p / abs(m3p)

  if(as.double(k * (phi + 1) + (phi + 2)^2) <= 0) {
    z <- (2 * (1 - mprime1p) * m2p) / m3p
    new <- m2p / (1 - mprime1p)^2
    betahat <- (z * (1 - new) - 2) / (1 - new + 2 * new * z)
    alphahat <- (new * betahat * (1 + betahat)) / (1 - new * betahat)
    ahat <- ((alphahat + betahat) * mprime1p - alphahat) / betahat
    bhat <- 1
  } else {
    theta <- (sgnm3p * phi * (phi + 2)) / sqrt(k * (phi + 1) + phi + 2)^2
    alphahat <- (phi - theta) / 2
    betahat <- (phi + theta) / 2
    ahat <- (mprime1p - sqrt((m2p * alphahat * (alphahat + betahat + 1)) / betahat))
    bhat <- (mprime1p + sqrt((m2p * betahat * (alphahat + betahat + 1)) / alphahat))
  }
  if ((bhat > 1) | (alphahat < 0) | (betahat < 0)) {
    z <- (2 * (1 - mprime1p) * m2p / m3p)
    new <- m2p / (1 - mprime1p)^2
    betahat <- (z * (1 - new) - 2) / (1 - new + 2 * new * z)
    alphahat <- (new * betahat * (1 + betahat)) / (1 - new * betahat)
    ahat <- ((alphahat + betahat) * mprime1p - alphahat) / betahat
    bhat <- 1
  }
  cutinp <- (cutscore - .5) / maxval
  cutinp2 <- c(ahat, cutinp, bhat)

  cut2 <- floor(((cutscore - .5) * n) / maxval)
  cut3 <- c(0, cut2, n)

  finalmatrix <- matrix(0, ncol = length(cut) + 1, nrow = length(cut) + 1)
  for(i in 1:nrow(finalmatrix)) {
    for(j in 1:nrow(finalmatrix)) {
      finalmatrix[i, j] <- integrate(ffunc, cutinp2[i], cutinp2[i + 1],
                                     ahat1 = ahat, bhat1 = bhat, alphahat1 = alphahat,
                                     betahat1 = betahat, lowbinom = cut3[j], highbinom = cut3[j + 1], n1 = n)$value
    }
  }

  if (adjust == 1) {
    proportions <- c(0, length(cutscore) + 1)
    xlowest <- x[x < cutscore[1]]
    proportions[1] <- length(xlowest) / length(x)
    if(length(cutscore) > 1) {
      for(i in 2:length(cutscore)) {
        temp <- x[x < cut[i]]
        proportions[i] <- (length(temp) / length(x) - sum(proportions[1:(i - 1)]))
      }
    }
    proportions[length(cutscore) + 1] <- 1 - sum(proportions[1:length(cut)])
    for (i in 1:ncol(finalmatrix)) {
      finalmatrix[, i] <- (finalmatrix[, i] * proportions[i]) / sum(finalmatrix[, i])
    }
  }

  list(alphahat = alphahat, betahat = betahat, ahat = ahat, bhat = bhat, n = n, finalmatrix = finalmatrix)
}
