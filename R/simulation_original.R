library(flexsurv)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(truncnorm)
library(TruncatedNormal)
library(magrittr)

projdir <- "/proj/ibrahimlab/changepoint"
stanfile <- file.path(projdir, "yqu/stan", "chgpt_weibull.stan")
stanmod <- cmdstanr::cmdstan_model(stanfile)
datdir <- "~/Projects/changepoint/"

survdat.original <- data.frame(fread(file.path(datdir, "data_pseudo_surv.csv")))
trt.survdat.original <- survdat.original %>% filter(ARM==1)
trt.cen.prop <- 1-mean(trt.survdat.original$PFS_EVENT)

ctrl.survdat.original <- survdat.original %>% filter(ARM==0)
ctrl.cen.prop <- 1-mean(ctrl.survdat.original$PFS_EVENT)

n <- 1000

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


savedir <- "/proj/ibrahimlab/changepoint/yqu/results"
standat <- readRDS(file.path(savedir, "standat.rds"))

trt.fit <- readRDS(file.path(savedir, "treat_draws_withcen_sd_01.rds"))
ctrl.fit <- readRDS(file.path(savedir, "control_draws_withcen_sd_01.rds"))

trt.standata <- standat[[2]]
ctrl.standata <- standat[[1]]
remove(list = "standat")

ttepars <- c("beta_tte", "shape_tte", "scale_tte")
ypars <- c("beta_y", "sd_y")
randeffpars <- c("chgpt_mean", "b0_mean", "b1_mean", "b2_mean", "chgpt_sd", "b_sd")
modelpars <- c(ttepars, ypars, randeffpars)


# ctrl.fit$draws(modelpars) %>% mcmc_trace()
# ctrl.fit$draws(modelpars) %>% mcmc_hist()
# trt.fit$draws(modelpars) %>% mcmc_trace()
# trt.fit$draws(modelpars) %>% mcmc_hist()
#
# ctrl.fit$draws(modelpars) %>% summarize_draws()
# trt.fit$draws(modelpars) %>% summarize_draws()

ctrl.pars.mean <- ctrl.fit$draws(c(modelpars, "randeff_corr")) %>% summarize_draws("mean")
trt.pars.mean <- trt.fit$draws(c(modelpars, "randeff_corr")) %>% summarize_draws("mean")
set.seed(1)
ctrl.id.sel <- sample(ctrl.standata$idobs, size = n, replace = TRUE)
set.seed(1)
trt.id.sel <- sample(trt.standata$idobs, size = n, replace = TRUE)

trt.survdat <- data.frame("id.orig.dat" = trt.id.sel, "Y0SCALE" = trt.standata$Xtte_obs[trt.id.sel, ])
ctrl.survdat <- data.frame("id.orig.dat" = ctrl.id.sel, "Y0SCALE" = ctrl.standata$Xtte_obs[trt.id.sel, ])

## GENERATE TTE DATA
ctrl.shape <- ctrl.pars.mean$mean[ctrl.pars.mean$variable == "shape_tte"]
ctrl.blscale <- ctrl.pars.mean$mean[ctrl.pars.mean$variable == "scale_tte"]
ctrl.betatte <- ctrl.pars.mean$mean[ctrl.pars.mean$variable == "beta_tte[1]"]
ctrl.scale <- ctrl.blscale * exp(ctrl.survdat$Y0SCALE * ctrl.betatte)
ctrl.survdat <- data.frame("id.orig.dat" = ctrl.id.sel, "Y0SCALE" = ctrl.standata$Xtte_obs[ctrl.id.sel, ])
ctrl.X <- matrix(ctrl.survdat$Y0SCALE, ncol = 1)
ctrl.survdat$PFS_YEARS <- rweibullph(ctrl.X, ctrl.betatte, ctrl.shape, ctrl.blscale)
ctrl.survdat$nvisits <- floor(ctrl.survdat$PFS_YEARS / 0.17)
ctrl.survdat <- rbind(ctrl.survdat %>% filter(nvisits > 0), ctrl.survdat %>% filter(nvisits == 0))
ctrl.survdat$id <- 1:nrow(ctrl.survdat)


#ctrl.e.survtime <-

## GENERATE RANDOM EFFECTS
ctrl.randeff.mean <- ctrl.fit$draws(
  c("chgpt_mean", "b0_mean", "b1_mean", "b2_mean")
) %>%
  summarize_draws("mean") %>%
  select(mean) %>%
  unlist()
#ctrl.randeff.mean[c(1, 3, 4)] <- c(0.5, -1, 1)

ctrl.randeff.sd <- ctrl.fit$draws(
  c("chgpt_sd", "b_sd")
) %>%
  summarize_draws("mean") %>%
  select(mean) %>%
  unlist()
#ctrl.randeff.sd <- rep(0.2, 4)

ctrl.randeff.corr <- ctrl.fit$draws(
  c("randeff_corr")
) %>%
  summarize_draws("mean") %>%
  select(mean) %>%
  unlist()
ctrl.randeff.corr <- matrix(ctrl.randeff.corr, 4, 4)
ctrl.randeff.cov <- diag(ctrl.randeff.sd) %*% ctrl.randeff.corr %*% diag(ctrl.randeff.sd)
ctrl.randeff <- rptmvn(ctrl.randeff.mean, ctrl.randeff.cov, 0, ctrl.survdat$PFS_YEARS)
ctrl.survdat <- cbind(ctrl.survdat, ctrl.randeff)

## GENERATE LONGITUDINAL DATA
ctrl.beta_y <- ctrl.fit$draws("beta_y") %>%
  summarize_draws(mean) %>%
  select(mean) %>%
  unlist()
ctrl.sd_y <- ctrl.fit$draws("sd_y") %>%
  summarize_draws(mean) %>%
  select(mean) %>%
  unlist()
ctrl.longdat <- ctrl.survdat %>% filter(nvisits > 0) %>% select(-c(PFS_YEARS))
rep.indx <- lapply(1:nrow(ctrl.longdat), function(i) rep(ctrl.longdat$id[i], ctrl.longdat$nvisits[i]))
rep.indx <- unlist(rep.indx)
ctrl.longdat <- ctrl.longdat[rep.indx, ] %>%
  group_by(id) %>%
  mutate(visitnum = row_number()) %>%
  mutate(visittime = 60 / 365 * visitnum) %>%
  mutate(
    delta = visittime - omega,
    eta_re = b1 + if_else(delta < 0, b2 * delta, b3 * delta),
    eta_fe = ctrl.beta_y * Y0SCALE,
    y = rnorm(row_number(), mean = eta_fe + eta_re, sd = ctrl.sd_y)
  )

## Remove TTE measurements without longitudinal data
ctrl.survdat <- ctrl.survdat %>% filter(id %in% unique(ctrl.longdat$id))


# ctrl.standata.new <- ctrl.standata
# ctrl.standata.new$nobs <- nrow(ctrl.survdat)
# ctrl.standata.new$Nobs <- nrow(ctrl.longdat)
# ctrl.standata.new$yobs <- ctrl.longdat$y
# ctrl.standata.new$tobs <- ctrl.survdat$PFS_YEARS
# ctrl.standata.new$visitobs <- ctrl.longdat$visittime
# ctrl.standata.new$Xlong_obs <- matrix(ctrl.longdat$Y0SCALE, ncol = 1)
# ctrl.standata.new$Xtte_obs <- matrix(ctrl.survdat$Y0SCALE, ncol = 1)
# ctrl.standata.new$idobs <- ctrl.longdat$id
# ctrl.standata.new$randeff_location <- c(0.5, 0.0, -0.5, 0.0)
# ctrl.standata.new$randeff_scale <- c(0.5, 1.0, 0.5, 1.0)
# ctrl.standata.new$randeff_shape <- c(8.0, 2.0, 8.0, 2.0)
#


## GENERATE TTE DATA for Treatment
trt.shape <- trt.pars.mean$mean[trt.pars.mean$variable == "shape_tte"]
trt.blscale <- trt.pars.mean$mean[trt.pars.mean$variable == "scale_tte"]
trt.betatte <- trt.pars.mean$mean[trt.pars.mean$variable == "beta_tte[1]"]


#trt.shape <- ctrl.shape
trt.blscale <- ctrl.blscale
#trt.betatte <- ctrl.betatte
trt.scale <- trt.blscale * exp(trt.survdat$Y0SCALE * trt.betatte)
trt.X <- matrix(trt.survdat$Y0SCALE, ncol = 1)
trt.survdat$PFS_YEARS <- rweibullph(trt.X, trt.betatte, trt.shape, trt.blscale)
trt.survdat$nvisits <- floor(trt.survdat$PFS_YEARS / 0.17)
trt.survdat <- rbind(trt.survdat %>% filter(nvisits > 0), trt.survdat %>% filter(nvisits == 0))
trt.survdat$id <- max(ctrl.survdat$id) + 1:nrow(trt.survdat)



## GENERATE RANDOM EFFECTS for Treatment
trt.randeff.mean <- trt.fit$draws(
  c("chgpt_mean", "b0_mean", "b1_mean", "b2_mean")
) %>%
  summarize_draws("mean") %>%
  select(mean) %>%
  unlist()
trt.randeff.mean[c(1, 3, 4)] <- c(0.5, -1, 1)


trt.randeff.sd <- trt.fit$draws(
  c("chgpt_sd", "b_sd")
) %>%
  summarize_draws("mean") %>%
  select(mean) %>%
  unlist()
trt.randeff.sd <- rep(0.2, 4)


trt.randeff.corr <- trt.fit$draws(
  c("randeff_corr")
) %>%
  summarize_draws("mean") %>%
  select(mean) %>%
  unlist()
trt.randeff.corr <- matrix(trt.randeff.corr, 4, 4)
trt.randeff.cov <- diag(trt.randeff.sd) %*% trt.randeff.corr %*% diag(trt.randeff.sd)
trt.randeff <- rptmvn(trt.randeff.mean, trt.randeff.cov, 0, trt.survdat$PFS_YEARS)
trt.survdat <- cbind(trt.survdat, trt.randeff)

## GENERATE LONGITUDINAL DATA for Treatment
trt.beta_y <- trt.fit$draws("beta_y") %>%
  summarize_draws(mean) %>%
  select(mean) %>%
  unlist()
trt.sd_y <- trt.fit$draws("sd_y") %>%
  summarize_draws(mean) %>%
  select(mean) %>%
  unlist()
trt.longdat <- trt.survdat %>% filter(nvisits > 0) %>% select(-c(PFS_YEARS))
rep.indx <- lapply(1:nrow(trt.longdat), function(i) rep(trt.longdat$id[i], trt.longdat$nvisits[i]))
rep.indx <- unlist(rep.indx) - max(ctrl.survdat$id)
trt.longdat <- trt.longdat[rep.indx, ] %>%
  group_by(id) %>%
  mutate(visitnum = row_number()) %>%
  mutate(visittime = 60 / 365 * visitnum) %>%
  mutate(
    delta = visittime - omega,
    eta_re = b1 + if_else(delta < 0, b2 * delta, b3 * delta),
    eta_fe = trt.beta_y * Y0SCALE,
    y = rnorm(row_number(), mean = eta_fe + eta_re, sd = trt.sd_y)
  )

## Remove TTE measurements without longitudinal data for Treatment
trt.survdat <- trt.survdat %>% filter(id %in% unique(trt.longdat$id))

## Prepare data for Stan for Treatment
# trt.standata.new <- trt.standata
# trt.standata.new$nobs <- nrow(trt.survdat)
# trt.standata.new$Nobs <- nrow(trt.longdat)
# trt.standata.new$yobs <- trt.longdat$y
# trt.standata.new$tobs <- trt.survdat$PFS_YEARS
# trt.standata.new$visitobs <- trt.longdat$visittime
# trt.standata.new$Xlong_obs <- matrix(trt.longdat$Y0SCALE, ncol = 1)
# trt.standata.new$Xtte_obs <- matrix(trt.survdat$Y0SCALE, ncol = 1)
# trt.standata.new$idobs <- trt.longdat$id
# trt.standata.new$randeff_location <- c(0.5, 0.0, -0.5, 0.0)
# trt.standata.new$randeff_scale <- c(0.5, 1.0, 0.5, 1.0)
# trt.standata.new$randeff_shape <- c(8.0, 2.0, 8.0, 2.0)

# find the mean censor time
mean.ctrl.tcen <- mean(ctrl.standata$tcen)
sd.ctrl.tcem <- sd(ctrl.standata$tcen)
mean.trt.tcen <- mean(trt.standata$tcen)
sd.trt.tcem <- sd(trt.standata$tcen)


# Generate individual censoring times from a normal distribution expontial

set.seed(1)
ctrl.survdat$censor_time <- rtruncnorm(nrow(ctrl.survdat), a = 0, mean = mean.ctrl.tcen, sd = sd.ctrl.tcem)
set.seed(1)
trt.survdat$censor_time <- rtruncnorm(nrow(trt.survdat), a = 0, mean = mean.trt.tcen, sd = sd.trt.tcem)

# Update PFS_YEARS and PFS_EVENT based on the generated censoring times
# For Control Dataset
ctrl.survdat$PFS_EVENT <- ifelse(ctrl.survdat$PFS_YEARS <= ctrl.survdat$censor_time, 1, 0) # 1 for event, 0 for censoring
ctrl.survdat$PFS_YEARS <- pmin(ctrl.survdat$PFS_YEARS, ctrl.survdat$censor_time)

# For Treatment Dataset
trt.survdat$PFS_EVENT <- ifelse(trt.survdat$PFS_YEARS <= trt.survdat$censor_time, 1, 0) # 1 for event, 0 for censoring
trt.survdat$PFS_YEARS <- pmin(trt.survdat$PFS_YEARS, trt.survdat$censor_time)

real_randeff_mean <- c(trt.randeff.mean, ctrl.randeff.mean)

# fit.ctrl.sim <- stanmod$sample(
#   data = ctrl.standata.new, iter_warmup = 2000, iter_sampling = 2000,
#   chains = 9, parallel_chains = 9
# )
#
# fit.trt.sim <- stanmod$sample(
#   data = trt.standata.new, iter_warmup = 2000, iter_sampling = 2000,
#   chains = 9, parallel_chains = 9
# )

ctrl.survdat$ARM <- 0
trt.survdat$ARM <- 1
ctrl.longdat$ARM <- 0
trt.longdat$ARM <- 1

# Append simulated data to the original datasets
survdat <- rbind(ctrl.survdat, trt.survdat)
longdat <- rbind(ctrl.longdat, trt.longdat)

longdat <- left_join(longdat, survdat %>% select(id, PFS_EVENT, PFS_YEARS), by = "id")

# Make the maximum of long data less than the event/censor time
longdat %<>% filter(visittime <= PFS_YEARS)

fmla.tte <- ~ 0 + Y0SCALE
fmla.long <- ~ 0 + Y0SCALE

standat.list <- list()

for (a in 0:1) {
  # Filter the data based on the current arm (a)
  survdat_arm <- survdat %>% filter(ARM == a)
  longdat_arm <- longdat %>% filter(ARM == a)

  # Separate out observed events and censored observations
  longdat.obs <- longdat_arm %>%
    filter(PFS_EVENT == TRUE)
  longdat.obs$id.new = as.numeric(as.factor(longdat.obs$id))
  longdat.obs %<>% arrange(id.new)

  longdat.cen <- longdat_arm %>%
    filter(PFS_EVENT == FALSE)

  # remove some censored individuals !!! do not forget to remove it for real simulation
  # longdat.cen %<>% filter(id < 930)
  longdat.cen$id.new = as.numeric(as.factor(longdat.cen$id))
  longdat.cen %<>% arrange(id.new)

  # Extract unique patient IDs
  iddata.obs <- unique(select(longdat.obs, id, id.new))
  iddata.cen <- unique(select(longdat.cen, id, id.new))

  # Get survival data for these unique patients
  survdat.obs <- merge(survdat_arm, iddata.obs, by = "id")
  survdat.cen <- merge(survdat_arm, iddata.cen, by = "id")

  # Prepare data for Stan
  standat.list[[a + 1]] <- list(
    nobs = nrow(survdat.obs),
    Nobs = nrow(longdat.obs),
    plong = ncol(model.matrix(fmla.long, data = longdat.obs)),
    ptte = ncol(model.matrix(fmla.tte, data = survdat.obs)),
    idobs = longdat.obs$id.new,
    tobs = survdat.obs$PFS_YEARS,
    yobs = longdat.obs$y,
    visitobs = longdat.obs$visittime,
    Xlong_obs = model.matrix(fmla.long, data = longdat.obs),
    Xtte_obs = model.matrix(fmla.tte, data = survdat.obs),
    ncen = nrow(survdat.cen),
    Ncen = nrow(longdat.cen),
    idcen = longdat.cen$id.new,
    tcen = survdat.cen$PFS_YEARS,
    ycen = longdat.cen$y,
    visitcen = longdat.cen$visittime,
    Xlong_cen = model.matrix(fmla.long, data = longdat.cen),
    Xtte_cen = model.matrix(fmla.tte, data = survdat.cen),
    randeff_location = c(0.5, 0.0, -0.5, 0.0),
    randeff_scale = c(0.5, 1.0, 0.5, 1.0),
    randeff_shape = c(8.0, 2.0, 8.0, 2.0)
  )
}

saveRDS(standat.list, file.path(projdir, "yqu/data", "standat.list.RDS"))

ggplot(
  data = longdat %>% filter(ARM==0), aes(x = visittime, y = y, color = factor(id))
) + geom_line() +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('Spaghetti plot of simulated data') +
  ylab('Tumor burden') +
  xlab('Time (years)')


# Renaming for clarity
longdat <- longdat %>% rename(YEAR = visittime, PCHG = y)

library(ggplot2)

# Define a custom labeller
custom_labeller <- function(variable, value){
  if(variable == "ARM"){
    return(ifelse(value == 0, "Control arm", "Treatment arm"))
  }
  if(variable == "PFS_EVENT"){
    return(ifelse(value == 0, "Censor", "Non-censor"))
  }
}

# Plot
ggplot(data = longdat, aes(x = YEAR, y = PCHG, color = factor(id))) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(ARM ~ PFS_EVENT, ncol = 2, labeller = custom_labeller) +
  ggtitle("Spaghetti plot of simulated data") +
  ylab("Tumor burden") +
  xlab("Time (years)")


# Filtering for observed and censored data for both treatment and control groups
longdat.obs.trt <- longdat %>%
  filter(ARM == 1, PFS_EVENT == 1)
longdat.obs.trt$id.new = as.numeric(as.factor(longdat.obs.trt$id))
longdat.obs.trt %<>% arrange(id.new)

longdat.obs.ctrl <- longdat %>%
  filter(ARM == 0, PFS_EVENT == 1)
longdat.obs.ctrl$id.new = as.numeric(as.factor(longdat.obs.ctrl$id))
longdat.obs.ctrl %<>% arrange(id.new)


longdat.cen.trt <- longdat %>%
  filter(ARM == 1, PFS_EVENT == 0)
longdat.cen.trt$id.new = as.numeric(as.factor(longdat.cen.trt$id))
longdat.cen.trt %<>% arrange(id.new)

longdat.cen.ctrl <- longdat %>%
  filter(ARM == 0, PFS_EVENT == 0)
longdat.cen.ctrl$id.new = as.numeric(as.factor(longdat.cen.ctrl$id))
longdat.cen.ctrl %<>% arrange(id.new)

iddata.obs.trt <- with(longdat.obs.trt, data.frame(id = unique(id), id.new = unique(id.new)))
iddata.obs.ctrl <- with(longdat.obs.ctrl, data.frame(id = unique(id), id.new = unique(id.new)))

iddata.cen.trt <- with(longdat.cen.trt, data.frame(id = unique(id), id.new = unique(id.new)))
iddata.cen.ctrl <- with(longdat.cen.ctrl, data.frame(id = unique(id), id.new = unique(id.new)))

survdat.obs.trt <- merge(survdat, iddata.obs.trt)
survdat.obs.ctrl <- merge(survdat, iddata.obs.ctrl)
survdat.cen.trt <- merge(survdat, iddata.cen.trt)
survdat.cen.ctrl <- merge(survdat, iddata.cen.ctrl)

fmla.tte <- ~ 0 + Y0SCALE
fmla.long <- ~ 0 + Y0SCALE + ARM

# Creating the data list for Stan
standat_twoarm <- list(
  nobs_trt = nrow(survdat.obs.trt),
  nobs_ctrl = nrow(survdat.obs.ctrl),
  Nobs_trt = nrow(longdat.obs.trt),
  Nobs_ctrl = nrow(longdat.obs.ctrl),
  plong = ncol(model.matrix(fmla.long, data = rbind(longdat.obs.trt, longdat.obs.ctrl))),
  ptte = ncol(model.matrix(fmla.tte, data = rbind(survdat.obs.trt, survdat.obs.ctrl))),
  idobs_trt = longdat.obs.trt$id.new,
  idobs_ctrl = longdat.obs.ctrl$id.new,
  tobs_trt = survdat.obs.trt$PFS_YEARS,
  tobs_ctrl = survdat.obs.ctrl$PFS_YEARS,
  yobs_trt = longdat.obs.trt$PCHG,
  yobs_ctrl = longdat.obs.ctrl$PCHG,
  visitobs_trt = longdat.obs.trt$YEAR,
  visitobs_ctrl = longdat.obs.ctrl$YEAR,
  Xlong_obs_trt = model.matrix(fmla.long, data = longdat.obs.trt),
  Xlong_obs_ctrl = model.matrix(fmla.long, data = longdat.obs.ctrl),
  Xtte_obs_trt = model.matrix(fmla.tte, data = survdat.obs.trt),
  Xtte_obs_ctrl = model.matrix(fmla.tte, data = survdat.obs.ctrl),
  ncen_trt = nrow(survdat.cen.trt),
  ncen_ctrl = nrow(survdat.cen.ctrl),
  Ncen_trt = nrow(longdat.cen.trt),
  Ncen_ctrl = nrow(longdat.cen.ctrl),
  idcen_trt = longdat.cen.trt$id.new,
  idcen_ctrl = longdat.cen.ctrl$id.new,
  tcen_trt = survdat.cen.trt$PFS_YEARS,
  tcen_ctrl = survdat.cen.ctrl$PFS_YEARS,
  ycen_trt = longdat.cen.trt$PCHG,
  ycen_ctrl = longdat.cen.ctrl$PCHG,
  visitcen_trt = longdat.cen.trt$YEAR,
  visitcen_ctrl = longdat.cen.ctrl$YEAR,
  Xlong_cen_trt = model.matrix(fmla.long, data = longdat.cen.trt),
  Xlong_cen_ctrl = model.matrix(fmla.long, data = longdat.cen.ctrl),
  Xtte_cen_trt = model.matrix(fmla.tte, data = survdat.cen.trt),
  Xtte_cen_ctrl = model.matrix(fmla.tte, data = survdat.cen.ctrl),
  randeff_location = c(0.5, 0.0, -0.5, 0.0),
  randeff_scale = c(0.5, 1.0, 0.5, 1.0),
  randeff_shape = c(8.0, 2.0, 8.0, 2.0)

  # randeff_location = c(0.5, 0.0, -0.5, 0.5),
  # randeff_scale = c(0.5, 1.0, 0.5, 0.5),
  # randeff_shape = c(8.0, 8.0, 8.0, 8.0)
)

saveRDS(standat_twoarm, file.path(projdir, "yqu/data", "standat.list.twoarm.RDS"))


# data.frame(
#   fit.ctrl.sim$draws(c(modelpars, "randeff_corr")) %>% summarize_draws(mean, ~ quantile(.x, probs = c(0.025, 0.975))),
#   "truth" = (fit.ctrl$draws(c(modelpars, "randeff_corr")) %>% summarize_draws("truth" = mean) %>% select(truth))
# ) %>%
#   mutate_if(is.numeric, round, 4)
#
