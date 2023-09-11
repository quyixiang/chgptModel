functions {
  //' Log likelihood of Weibull PH model
  //'
  //' @param y positive event times
  //' @param X design matrix (no intercept term)
  //' @param beta regression coefficients (no intercept)
  //' @param shape shape parameter for Weibull distribution
  //' @param scale scale parameter for Weibull distribution
  //'
  //' @return log density
  real weibull_ph_lpdf(
    vector y, matrix X, vector beta, real shape, real scale
  ) {
    int nobs = size(y);
    vector[nobs] log_y = log(y);
    real log_scale = log(scale);
    real log_shape = log(shape);
    vector[nobs] eta = X * beta;
    return sum(
        -exp( log_scale + eta + shape * log_y )            // log survival = -cumhaz
      + log_scale + eta + log_shape + (shape - 1) * log_y  // log hazard
    );
  }
  //' Generalized normal density
  //'
  //' Density function for the generalized normal distribution:
  //' f(x) \propto exp( [ |x - mu| / scale ]^shape ). If scale = 1 and shape = 2,
  //' corresponds to normal distribution
  //'
  //' @param x the random variable
  //' @param mu location parameter
  //' @param scale scale parameter
  //' @param shape shape parameter
  //'
  //' @return log density
  real generalized_normal_lpdf(real x, real mu, real scale, real shape) {
    return   log(shape) - log(scale) - 0.69314718056 - lgamma(inv(shape))
           - pow( abs(x - mu) * inv(scale) , shape );
  }
  //' Truncated normal PDF with vector of upper bounds (lower bound = 0)
  //'
  //' @param x vector of PTMVN random variables
  //' @param mu mean parameter of MVN distribution
  //' @param sigma covariance parameter of MVN distribution
  //' @parma ub upper bound of MVN distribution
  //'
  //' @return log density
  real tnormal_lpdf(vector x, real mu, real sigma, vector ub) {
    int nobs = size(ub);
    vector[nobs] upsilon = (ub - mu) * inv(sigma);
    real lambda = -mu * inv(sigma);
    real logDiffPhi;
    real res = normal_lpdf(x | mu, sigma);
    real logCDFlambda = std_normal_lcdf(lambda);
    real logCCDFlambda = std_normal_lccdf(lambda);
    for ( i in 1:nobs ) {
      if ( upsilon[i] < 0 )
        res -= log_diff_exp( std_normal_lcdf(upsilon[i]), logCDFlambda );
      else
        res -= log_diff_exp( logCCDFlambda, std_normal_lccdf(upsilon[i]));
    }
    return(res);
  }
}
data {
  int<lower=0>             nobs;             // number of observed event times
  int<lower=0>             Nobs;             // number of longitudinal measurements for observed event times
  int<lower=0>             plong;            // number of covariates for longitudinal outcome
  int<lower=0>             ptte;             // number of covariates for time-to-event outcome
  array[Nobs] int<lower=1> idobs;            // idobs variable for observed event times
  vector[Nobs]             yobs;             // longitudinal measurements for observed event times
  vector<lower=0>[Nobs]    visitobs;         // Visit times for observed individuals
  matrix[Nobs, plong]      Xlong_obs;        // Design matrix for longitudinal outcome (observed)
  vector<lower=0>[nobs]    tobs;             // Event times (observed)
  matrix[nobs, ptte]       Xtte_obs;         // Design matrix for TTE outcome (observed)
  vector[4]                randeff_location; // location parameter for generalized normal prior on random effects
  vector<lower=0>[4]       randeff_scale;    // scale parameter for generalized normal prior on random effects
  vector<lower=0>[4]       randeff_shape;    // shape parameter for generalized normal prior on random effects
}
parameters {
  // TTE parameters
  vector[ptte] beta_tte;              // regression coefficients for TTE outcome
  real<lower=0> shape_tte;            // shape parameter for TTE outcome
  real<lower=0> scale_tte;            // scale parameter for TTE outcome

  // Longitudinal parameters
  vector[plong]                   beta_y;             // regression coefficients for fixed effects
  real<lower=0>                   sd_y;               // standard deviation of longitudinal outcomes | randeff
  real<lower=0>                   chgpt_mean;         // mean of change point
  real<lower=0,upper=10>        chgpt_sd;           // SD of change point
  real                            b0_mean;            // mean of random intercept
  real                            b1_mean;            // mean of preslope
  real                            b2_mean;            // mean of postslope
  vector<lower=0,upper=100>[3]    b_sd;               // SD of random effects assoc. with regression
  cholesky_factor_corr[4]         randeff_corr_chol;  // Cholesky decomposition of correlation of random effects
  vector<lower=0,upper=1>[nobs]   chgpt_obs_raw;      // raw version of change point random effect U(0, 1)
  matrix[3, nobs]                 b_obs_raw;          // raw version of random coefficients (will be std. normal)
}
transformed parameters {
  vector[nobs] chgpt_obs;           // changepoint (observed times)
  matrix[3,nobs] b_obs;             // random intercept and slopes (observed event times)
  matrix[4,4] randeff_cov_chol;     // cholesky factor of covariance matrix
  vector[3]   LHS;                  // sigma[1,2:4] / sigma[1,1] = chol[1,2:4] / chol[1,1]
  vector[3]   b_mean;               // mean of random intercept, preslope, and postslope

  b_mean[1] = b0_mean;
  b_mean[2] = b1_mean;
  b_mean[3] = b2_mean;

  // Transform chgpt to get proper bounds
  chgpt_obs = tobs .* chgpt_obs_raw;       // in [0, tobs]

  // cholesky factor of covariance matrix of MVN prior to truncation
  randeff_cov_chol = diag_pre_multiply(append_row(chgpt_sd, b_sd), randeff_corr_chol);

  // Distribution of b_obs | chgpt_obs is MVN(condMean, condCov)
  LHS = randeff_cov_chol[2:4,1] * inv(randeff_cov_chol[1,1]);
  b_obs   = randeff_cov_chol[2:4,2:4] * b_obs_raw;             // implies b_obs | chgpt_obs ~ MVN(0, condCov)
  for ( i in 1:nobs )
    b_obs[, i] += b_mean + LHS * (chgpt_obs[i] - chgpt_mean);  // implies b_obs | chgpt_obs ~ MVN(condMean, condCov)
}
model {
  vector[Nobs] yobs_mean = Xlong_obs * beta_y;
  real diff_ij;

  // Compute longitudinal means | random effects
  for ( i in 1:Nobs ) {
    diff_ij = visitobs[i] - chgpt_obs[idobs[i]];
    yobs_mean[i] += b_obs[1,idobs[i]] + (diff_ij <= 0 ? b_obs[2,idobs[i]] : b_obs[3,idobs[i]] ) * diff_ij;
  }

  // Priors for TTE model
  beta_tte  ~ std_normal();
  shape_tte ~ std_normal();
  scale_tte ~ std_normal();

  // Priors for longitudinal model
  beta_y            ~ normal(0, 10);
  sd_y              ~ normal(0, 10);

  // Priors for random effects parameters
  randeff_corr_chol ~ lkj_corr_cholesky(1.0);   // implies corr ~ uniform over PD correlation matrices
  chgpt_mean ~ generalized_normal(randeff_location[1], randeff_scale[1], randeff_shape[1]);
  b0_mean    ~ generalized_normal(randeff_location[2], randeff_scale[2], randeff_shape[2]);
  b1_mean    ~ generalized_normal(randeff_location[3], randeff_scale[3], randeff_shape[3]);
  b2_mean    ~ generalized_normal(randeff_location[4], randeff_scale[4], randeff_shape[4]);
  chgpt_sd   ~ std_normal();
  b_sd       ~ std_normal();

  // Complete data likelihood (observed)
  tobs                 ~ weibull_ph(Xtte_obs, beta_tte, shape_tte, scale_tte);
  yobs                 ~ normal(yobs_mean, sd_y);  // conditional distribution of yobs | random effects
  chgpt_obs            ~ tnormal(chgpt_mean, randeff_cov_chol[1,1], tobs);
  to_vector(b_obs_raw) ~ std_normal();

}
generated quantities {
  corr_matrix[4] randeff_corr = multiply_lower_tri_self_transpose(randeff_corr_chol);
}


