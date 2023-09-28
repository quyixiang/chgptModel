#' Bayesian Analysis for Change Point Model using CmdStan
#'
#' This function performs Bayesian change point analysis on the given data using CmdStan. It leverages a Weibull distribution to identify change points in data with potential applications in time-to-event and longitudinal analyses.
#'
#' @param survdat A data frame containing the pseudo survival data. If not provided, the function will load the default dataset from the package.
#' @param longdat A data frame containing the pseudo longitudinal data. If not provided, the function will load the default dataset from the package.
#' @param fmla.tte A formula defining the structure for time-to-event analysis.
#' @param fmla.long A formula defining the structure for longitudinal analysis.
#' @param longdat.time A string specifying the name of time variable in the longitudinal dataset.
#' @param id.indicator A string specifying the identifier variable, default is "id".
#' @param non_censor_only A logical value indicating whether to use only non-censored data. Default is FALSE.
#' @param chgpt_location A numerical value defining the initial location for the change point in the Bayesian change point model. Default is 0.5.
#' @param b0_location A numerical value defining the initial location for the parameter `b0` in the model. Default is 0.
#' @param b1_location A numerical value defining the initial location for the parameter `b1` in the model. Default is -0.5.
#' @param b2_location A numerical value defining the initial location for the parameter `b2` in the model. Default is 0.5.
#'
#' @param chgpt_scale A numerical value defining the initial scale for the change point in the Bayesian change point model. Default is 0.5.
#' @param b0_scale A numerical value defining the initial scale for the parameter `b0` in the model. Default is 1.
#' @param b1_scale A numerical value defining the initial scale for the parameter `b1` in the model. Default is 0.5.
#' @param b2_scale A numerical value defining the initial scale for the parameter `b2` in the model. Default is 0.5.
#'
#' @param chgpt_shape A numerical value defining the initial shape for the change point in the Bayesian change point model. Default is 8.
#' @param b0_shape A numerical value defining the initial shape for the parameter `b0` in the model. Default is 8.
#' @param b1_shape A numerical value defining the initial shape for the parameter `b1` in the model. Default is 8.
#' @param b2_shape A numerical value defining the initial shape for the parameter `b2` in the model. Default is 8.
#' @param ... Additional arguments to be passed to the \code{cmdstanr::cmdstan_model$sample} function.
#'
#' @return A \code{cmdstanr::fit} object containing the results of the Bayesian analysis.
#'
#' @examples
#' # Load required datasets
#' data("survdat")
#' data("longdat")
#'
#' # Define a formula for time-to-event analysis
#' example_formula_tte <- as.formula(Surv(PFS_YEARS, PFS_EVENT) ~ 0 + Y0SCALE)
#'
#' # Define a formula for longitudinal analysis
#' example_formula_long <- as.formula(PCHG ~ 0 + Y0SCALE)
#'
#' # Run the stan model
#' results <- chgptMCMC(
#'   survdat = survdat, longdat = longdat,
#'   fmla.tte = example_formula_tte, fmla.long = example_formula_long,
#'   longdat.time = "YEAR", id.indicator = "id",
#'   iter_warmup = 1000, iter_sampling = 1000,
#'   chains = 4, parallel_chains = 4, adapt_delta = 0.8, refresh = 200,
#'   non_censor_only = TRUE
#' )
#' # Posterior plot
#' results$draws(c("chgpt_mean", "chgpt_sd")) %>% mcmc_trace()
#' results$draws(c("chgpt_mean", "chgpt_sd")) %>% summarize_draws()
#' @export
chgptMCMC <- function(survdat,
                      longdat,
                      fmla.tte,
                      fmla.long,
                      longdat.time,
                      id.indicator = "id",
                      non_censor_only = FALSE,
                      chgpt_location = 0.5, b0_location = 0, b1_location = -0.5, b2_location = 0.5,
                      chgpt_scale = 0.5, b0_scale = 1, b1_scale = 0.5, b2_scale = 0.5,
                      chgpt_shape = 8, b0_shape = 8, b1_shape = 8, b2_shape = 8,
                      ...) {
  data.list <- load_onearm_data(
    survdat = survdat,
    longdat = longdat,
    fmla.tte = fmla.tte,
    fmla.long = fmla.long,
    longdat.time = longdat.time,
    id.indicator = id.indicator,
    chgpt_location = chgpt_location, b0_location = b0_location, b1_location = b1_location, b2_location = b2_location,
    chgpt_scale = chgpt_scale, b0_scale = b0_scale, b1_scale = b1_scale, b2_scale = b2_scale,
    chgpt_shape = chgpt_shape, b0_shape = b0_shape, b1_shape = b1_shape, b2_shape = b2_shape
  )
  if (non_censor_only) {
    stan_model_file <- system.file("stan", "chgpt_weibull_nocen.stan", package = "chgptModel")
  } else {
    stan_model_file <- system.file("stan", "chgpt_weibull.stan", package = "chgptModel")
  }
  model <- cmdstanr::cmdstan_model(stan_model_file)
  fit <- model$sample(
    data = data.list, ...
  )
  return(fit)
}
