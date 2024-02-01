# Change Point Model

## Install

1. Install `CmdStanR` following the instructions [here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).

2. Use the following command to install the package in R.

```
library(devtools)
install_github("quyixiang/ChgptModel")
```

## Example code

### Generate simulation datasets

```
simulation_data <- generate_simulation_data(
  fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE),
  fmla.long = as.formula(PCHG ~ 0 + Y0SCALE),
  bootstrap = FALSE, id.name = "id",
  beta.tte = c(0.2),
  scale.tte = 3, shape.tte = 2,
  normal.tte = FALSE,
  beta.y = c(-0.01), sd.y = 0.1,
  randeff.mean = c(1, 0, -0.5, 0.5), randeff.sd = c(0.2, 0.2, 0.2, 1),
  randeff.corr = matrix(
    c(
      1, -0.4, -0.2, -0.3,
      -0.4, 1, 0.5, 0.20,
      -0.2, 0.5, 1, 0.2,
      -0.3, 0.20, 0.2, 1
    ),
    nrow = 4
  ),
  n = 100, censor.parameter = 0.5, time.interval = 0.1, seed = 1
)


survdat <- simulation_data[["survdat"]]
longdat <- simulation_data[["longdat"]]
simulation.para <- simulation_data[["simulation.para"]]
```

### Use MCMC to sample the posterior distribution of the change point model

```
results <- chgptMCMC(
  survdat = survdat, longdat = longdat,
  fmla.tte = as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE), fmla.long = as.formula(PCHG ~ 0 + Y0SCALE),
  longdat.time = "visittime", id.indicator = "id",
  iter_warmup = 1000, iter_sampling = 2000,
  chains = 1, parallel_chains = 1, adapt_delta = 0.8, refresh = 200,
  non_censor_only = FALSE
)
```

### Compare posterior mean and real parameters

```
params <- c("beta_tte", "scale_tte", "shape_tte", "beta_y", "sd_y", "chgpt_mean", "b_mean", "chgpt_sd", "b_sd", "randeff_corr")
results_summary <- results$draws(params) %>% summarize_draws(mean, sd, ~ quantile(.x, probs = c(0.025, 0.975)))
results_summary$truth <- simulation.para
print(results_summary, n = 100)
```
