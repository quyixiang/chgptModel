library(TruncatedNormal)
library(randcorr)
library(dplyr)
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



# fmla.tte <- as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ Y0SCALE + Y1SCALE)
# fmla.long <- as.formula(PCHG ~ Y0SCALE + Y2SCALE)
#
# betatte <- c(2, 0.201, 0.5)
# scalette <- 2.41
# shapette <- 1.85
# beta_y <- c(1, -0.025, 1)
# sd_y <- 0.07
# n <- 100

#' Generate Simulation Data for Time-to-Event and Longitudinal Data Analyses
#'
#' This function generates simulated data suitable for time-to-event and longitudinal data analyses,
#' based on the provided parameters and formulae.
#'
#' @param fmla.tte A formula object specifying the structure of the time-to-event model.
#' @param fmla.long A formula object specifying the structure of the longitudinal data model.
#' @param betatte A numeric vector specifying the beta parameters for the TTE model.
#' @param scalette A numeric value specifying the scale parameter for the TTE model.
#' @param shapette A numeric value specifying the shape parameter for the TTE model.
#' @param beta_y A numeric vector specifying the beta parameters for the longitudinal data model.
#' @param sd_y A numeric value specifying the standard deviation of the Y variable in the longitudinal data model.
#' @param randeff.mean A numeric vector specifying the mean of the random effects.
#' @param randeff.sd A numeric vector specifying the standard deviations of the random effects.
#' @param randeff.corr A matrix specifying the correlation structure of the random effects, or NULL to generate a random structure.
#' @param n An integer specifying the number of observations to simulate.
#' @param censor.parameter A numeric value specifying the parameter for the exponential distribution used to generate censoring times.
#' @param time.interval A numeric value specifying the time interval for visits in the longitudinal data.
#' @param seed An integer specifying the seed for random number generation.
#'
#' @return A list containing three elements:
#' \itemize{
#'  \item \code{survdat} A data frame containing the simulated survival data.
#'  \item \code{longdat} A data frame containing the simulated longitudinal data.
#'  \item \code{simulation.para} A list containing the parameters used in the simulation.
#' }
#'
#' @examples
#' simulation.list <- generate_simulation_data(
#'   fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ Y0SCALE + Y1SCALE),
#'   fmla.long = as.formula(PCHG ~ 0 + Y0SCALE + Y2SCALE + Y3),
#'   betatte = c(0.1, 0.05, 0.1), scalette = 2, shapette = 2,
#'   beta_y = c(0.02, -0.02, 0.03), sd_y = 0.1,
#'   randeff.mean = c(0.5, 0, -1, 1), randeff.sd = rep(0.2, 4), randeff.corr = NULL,
#'   n = 1000, censor.parameter = 2, time.interval = 0.1
#' )
#' @export
generate_simulation_data <- function(
    fmla.tte, fmla.long,
    betatte, scalette = NULL, shapette = NULL,
    normal.tte = FALSE, sd_tte = NULL, # sd_tte should not be null if normal.tte=TRUE
    beta_y, sd_y,
    randeff.mean = c(0.5, 0, -1, 1), randeff.sd = rep(0.2, 4), randeff.corr = NULL,
    n = 100, censor.parameter, time.interval = 0.1,
    seed = 1) {
  # Set the seed for random number generation
  set.seed(seed)
  n_final <- n
  n <- 3 * n
  # First generate X
  # Get variable names excluding the response variable
  Xtte.name <- all.vars(fmla.tte)[-c(1:2)]
  Xlong.name <- all.vars(fmla.long)[-1]
  Xtte.response.name <- all.vars(fmla.tte)[1]
  Xtte.indicator.name <- all.vars(fmla.tte)[2]
  Xlong.response.name <- all.vars(fmla.long)[1]
  Xall.name <- unique(c(Xtte.name, Xlong.name))

  # Create a matrix with values from a standard normal distribution
  X <- matrix(rnorm(n * length(Xall.name), mean = 0, sd = 1), nrow = n, ncol = length(Xall.name), byrow = FALSE)

  # Assign column names from Xtte.name to the matrix
  colnames(X) <- Xall.name
  X <- data.frame(X)
  fmla.tte.rhs <- as.formula(paste("~", deparse(fmla.tte[[3]])))
  fmla.long.rhs <- as.formula(paste("~", deparse(fmla.long[[3]])))
  Xtte <- model.matrix(fmla.tte.rhs, data = X)
  Xlong <- model.matrix(fmla.long.rhs, data = X)
  Xlong.name.intercept <- colnames(Xlong)
  X_combined <- cbind(Xtte, Xlong)
  X_combined_df <- as.data.frame(X_combined)
  X_combined_unique <- X_combined_df %>%
    select(unique(colnames(X_combined_df)))
  X.model <- as.matrix(X_combined_unique)

  # Generate TTE data
  id <- c(1:n)
  survdat <- as.data.frame(X.model)
  survdat$id <- id
  survdat <- survdat[, c("id", names(survdat)[!names(survdat) %in% "id"])]
  if (normal.tte) {
    survdat$event_years <- exp(rnorm(Xtte %*% betatte, sd_tte))
  } else {
    survdat$event_years <- rweibullph(Xtte, betatte, scalette, scalette)
  }
  survdat$nvisits <- ceiling(survdat$event_years / time.interval)

  survdat <- survdat %>% filter(nvisits > 0)
  survdat$id <- 1:nrow(survdat)

  # Add censor data (the name of the outcome should be the same with the fmla.tte)
  set.seed(seed)
  survdat$censor_years <- rexp(nrow(survdat), censor.parameter)
  survdat[[Xtte.indicator.name]] <- ifelse(survdat$event_years <= survdat$censor_years, 1, 0)
  survdat[[Xtte.response.name]] <- pmin(survdat$censor_years, survdat$event_years)



  # Generate random effects
  if (!is.null(randeff.corr)) {
    randeff.corr <- matrix(randeff.corr, nrow = 4, ncol = 4)
  } else {
    set.seed(seed)
    randeff.corr <- randcorr(4)
  }

  randeff.cov <- diag(randeff.sd) %*% randeff.corr %*% diag(randeff.sd)
  randeff <- rptmvn(randeff.mean, randeff.cov, 0, survdat$event_years)
  survdat <- cbind(survdat, randeff)

  # Generate longitudinal data
  longdat <- survdat %>%
    filter(nvisits > 0)
  rep.indx <- lapply(1:nrow(longdat), function(i) rep(longdat$id[i], longdat$nvisits[i])) %>% unlist()

  set.seed(seed)
  longdat <- longdat[rep.indx, ] %>%
    group_by(id) %>%
    mutate(visitnum = row_number()) %>%
    rowwise() %>%
    mutate(rand_val = rnorm(1, mean = 0, sd = 0.02)) %>%
    mutate(visittime = time.interval * visitnum + rand_val) %>%
    ungroup() %>%
    mutate(
      delta = visittime - omega,
      eta_re = b1 + if_else(delta < 0, b2 * delta, b3 * delta)
    ) %>%
    rowwise() %>%
    mutate(
      eta_fe = sum(c_across(all_of(Xlong.name.intercept)) * beta_y),
      {{ Xlong.response.name }} := rnorm(1, mean = eta_fe + eta_re, sd = sd_y)
    )


  # Make the maximum of long data less than the event/censor time
  longdat <- longdat %>% filter(visittime <= !!sym(Xtte.response.name))
  longdat <- longdat %>%
    group_by(id) %>%
    mutate(visitnum_aftercensor = row_number())

  # we only keep these subjects
  final_id <- unique(longdat$id)[1:n_final]
  longdat <- longdat %>% filter(id %in% final_id)
  survdat <- survdat %>% filter(id %in% final_id)

  cat(
    paste0(
      "Number of Censor: ", sum(survdat[[Xtte.indicator.name]] == 0), "\nNumber of Observation: ", sum(survdat[[Xtte.indicator.name]] == 1),
      "\nProportion of Censor: ", sum(survdat[[Xtte.indicator.name]] == 0) / nrow(survdat), "\n"
    )
  )

  # We also need to save the initial parameters finally
  rownames(randeff.corr) <- c("omega", "b1", "b2", "b3")
  colnames(randeff.corr) <- c("omega", "b1", "b2", "b3")
  # randeff.corr.upper <- randeff.corr[upper.tri(randeff.corr, diag = TRUE)]
  # names(randeff.corr.upper) <- paste("randeff.corr[", rownames(randeff.corr)[row(randeff.corr)[upper.tri(randeff.corr, diag = TRUE)]],
  #   ",", colnames(randeff.corr)[col(randeff.corr)[upper.tri(randeff.corr, diag = TRUE)]], "]",
  #   sep = ""
  # )

  simulation.para <- c(
    betatte = betatte,
    scalette = scalette,
    shapette = shapette,
    beta_y = beta_y,
    sd_y = sd_y,
    randeff.mean = randeff.mean,
    randeff.sd = randeff.sd,
    randeff.corr = randeff.corr
  )
  return(list(survdat = survdat, longdat = longdat, simulation.para = simulation.para))
}

# # test
# simulation.list <- generate_simulation_data(
#   fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ Y0SCALE + Y1SCALE),
#   fmla.long = as.formula(PCHG ~ 0 + Y0SCALE + Y2SCALE + Y3),
#   betatte = c(0.1, 0.05, 0.1), scalette = 2, shapette = 2,
#   beta_y = c(0.02, -0.02, 0.03), sd_y = 0.1,
#   randeff.mean = c(0.5, 0, -1, 1), randeff.sd = rep(0.2, 4), randeff.corr = NULL,
#   n = 1000, censor.parameter = 2, time.interval = 0.1
# )
#
# survdat.simulation <- simulation.list[["survdat"]]
# longdat.simulation <- simulation.list[["longdat"]]
#
# real.parameter <- simulation.list[["simulation.para"]]
# real.parameter
#
# library(chgptModel)
# # simulation.standat.list <- load_onearm_data(survdat.simulation, longdat.simulation,
# #   as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ Y0SCALE + Y1SCALE),
# #   as.formula(PCHG ~ Y0SCALE + Y2SCALE + Y3),
# #   longdat.time = "visittime"
# # )
#
#
#
# simulation.results <- chgptMCMC(
#   fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ Y0SCALE + Y1SCALE),
#   fmla.long = as.formula(PCHG ~ 0 + Y0SCALE + Y2SCALE + Y3),
#   survdat = survdat.simulation,
#   longdat = longdat.simulation,
#   longdat.time = "visittime",
#   iter_warmup = 1000, iter_sampling = 1000,
#   chains = 4, parallel_chains = 4, adapt_delta = 0.8, refresh = 200
# )
#
# parms <- c("beta_tte", "scale_tte", "shape_tte", "beta_y", "sd_y","chgpt_mean", "b_mean", "chgpt_sd", "b_sd", "randeff_corr")
#
# simulation.results$draws(c("chgpt_mean", "chgpt_sd")) %>% mcmc_trace()
#
# simulation.draws <- simulation.results$draws(parms)
# dim(simulation.draws)
# first_chain_draws <- simulation.draws[,1,]
# mcmc_trace(first_chain_draws)
#
# simulation.results$draws(parms) %>% mcmc_hist()
# simulation.results$draws(c("chgpt_mean", "chgpt_sd")) %>% mcmc_acf()
# compare <- cbind(simulation.results$draws(parms) %>% summarize_draws(), real = unname(real.parameter))
