rptmvn <- function(mu, sigma, a, b) {
  omegaMean <- mu[1]
  omegaVar <- sigma[1, 1]
  omegaSD <- sqrt(omegaVar)
  ## Check dimensions
  na <- length(a)
  nb <- length(b)
  if (na == 1 & nb > 1) {
    a <- rep(a, nb)
  } else if (na > 1 & nb == 1) {
    b <- rep(b, na)
  } else if (na != nb) {
    stop("dimensions of a and b are not compatible")
  }
  n <- max(na, nb)
  ## Draw omega based on its marginal distribution
  omega <- rtnorm(1, omegaMean, omegaSD, a, b)
  ## Draw b based on its conditional distribution
  bCondCov <- sigma[-1, -1] - tcrossprod(sigma[-1, 1]) / omegaVar
  bCondMean <- t(mu[-1] + tcrossprod(sigma[-1, 1] / sigma[1, 1], omega - omegaMean))
  b <- bCondMean + mvtnorm::rmvnorm(n, sigma = bCondCov)
  colnames(b) <- paste0("b", 1:ncol(b))
  cbind("omega" = omega, "b" = b)
}

rweibullph <- function(X, beta, shape, scale) {
  eta <- (X %*% beta)[, 1]
  logscale <- eta + log(scale)
  scale <- exp(logscale)
  rexp(nrow(X), rate = scale)^(1 / shape)
}


generate_data <- function(id_sel, n, betatte, scale, shape, time_interval = 0.17) {
  library(dplyr)

  X <- matrix(survdat$Y0SCALE, ncol = 1)
  survdat$PFS_YEARS <- rweibullph(X, betatte, shape, scale)
  survdat$nvisits <- floor(survdat$PFS_YEARS / time_interval)

  survdat <- rbind(survdat %>% filter(nvisits > 0), survdat %>% filter(nvisits == 0))
  survdat$id <- 1:nrow(survdat)

  # Generate random effects
  # ... (rest of the random effects and longitudinal data generation code, using pars.mean, survdat, etc.)

  return(list(survival_data = survdat, longitudinal_data = longdat))
}
