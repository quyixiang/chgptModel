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
  int<lower=0>             nobs_trt;             // number of observed event times
  int<lower=0>             nobs_ctrl;             // number of observed event times
  int<lower=0>             Nobs_trt;             // number of longitudinal measurements for observed event times
  int<lower=0>             Nobs_ctrl;             // number of longitudinal measurements for observed event times
  int<lower=0>             plong;            // number of covariates for longitudinal outcome
  int<lower=0>             ptte;             // number of covariates for time-to-event outcome
  array[Nobs_trt] int<lower=1> idobs_trt;            // idobs variable for observed event times
  array[Nobs_ctrl] int<lower=1> idobs_ctrl;            // idobs variable for observed event times
  vector[Nobs_trt]             yobs_trt;             // longitudinal measurements for observed event times
  vector[Nobs_ctrl]             yobs_ctrl;             // longitudinal measurements for observed event times
  vector<lower=0>[Nobs_trt]    visitobs_trt;         // Visit times for observed individuals
  vector<lower=0>[Nobs_ctrl]    visitobs_ctrl;         // Visit times for observed individuals
  matrix[Nobs_trt, plong]      Xlong_obs_trt;        // Design matrix for longitudinal outcome (observed)
  matrix[Nobs_ctrl, plong]      Xlong_obs_ctrl;        // Design matrix for longitudinal outcome (observed)
  vector<lower=0>[nobs_trt]    tobs_trt;             // Event times (observed)
  vector<lower=0>[nobs_ctrl]    tobs_ctrl;             // Event times (observed)
  matrix[nobs_trt, ptte]       Xtte_obs_trt;         // Design matrix for TTE outcome (observed)
  matrix[nobs_ctrl, ptte]       Xtte_obs_ctrl;         // Design matrix for TTE outcome (observed)
  int<lower=0>             ncen_trt;             // number of censored event times
  int<lower=0>             ncen_ctrl;             // number of censored event times
  int<lower=0>             Ncen_trt;             // number of longitudinal measurements for censored event times
  int<lower=0>             Ncen_ctrl;             // number of longitudinal measurements for censored event times
  array[Ncen_trt] int<lower=1> idcen_trt;            // idcen variable for censored event times
  array[Ncen_ctrl] int<lower=1> idcen_ctrl;            // idcen variable for censored event times
  vector[Ncen_trt]             ycen_trt;             // longitudinal measurements for censored event times
  vector[Ncen_ctrl]             ycen_ctrl;             // longitudinal measurements for censored event times
  vector<lower=0>[Ncen_trt]    visitcen_trt;         // Visit times for censored individuals
  vector<lower=0>[Ncen_ctrl]    visitcen_ctrl;         // Visit times for censored individuals
  matrix[Ncen_trt, plong]      Xlong_cen_trt;        // Design matrix for longitudinal outcome (censored)
  matrix[Ncen_ctrl, plong]      Xlong_cen_ctrl;        // Design matrix for longitudinal outcome (censored)
  vector<lower=0>[ncen_trt]    tcen_trt;             // Event times (censored)
  vector<lower=0>[ncen_ctrl]    tcen_ctrl;             // Event times (censored)
  matrix[ncen_trt, ptte]       Xtte_cen_trt;         // Design matrix for TTE outcome (censored)
  matrix[ncen_ctrl, ptte]       Xtte_cen_ctrl;         // Design matrix for TTE outcome (censored)
  vector[4]                randeff_location; // location parameter for generalized normal prior on random effects
  vector<lower=0>[4]       randeff_scale;    // scale parameter for generalized normal prior on random effects
  vector<lower=0>[4]       randeff_shape;    // shape parameter for generalized normal prior on random effects
}
parameters {
  // TTE parameters
  vector[ptte] beta_tte_trt;              // regression coefficients for TTE outcome
  vector[ptte] beta_tte_ctrl;              // regression coefficients for TTE outcome
  real<lower=0> shape_tte_trt;            // shape parameter for TTE outcome
  real<lower=0> shape_tte_ctrl;            // shape parameter for TTE outcome
  real<lower=0> scale_tte_trt;            // scale parameter for TTE outcome
  real<lower=0> scale_tte_ctrl;            // scale parameter for TTE outcome
  vector<lower=tcen_trt>[ncen_trt] tcen_trtimp;   // Imputed event times for censored individuals
  vector<lower=tcen_ctrl>[ncen_ctrl] tcen_ctrlimp;   // Imputed event times for censored individuals
  
  // Longitudinal parameters
  vector[plong]                   beta_y;             // regression coefficients for fixed effects
  real<lower=0>                   sd_y;               // standard deviation of longitudinal outcomes | randeff
  real<lower=0>                   chgpt_mean_trt;         // mean of change point for treatment
  real<lower=0>                   chgpt_mean_ctrl;         // mean of change point for control
  real<lower=0,upper=10>        chgpt_sd;           // SD of change point
  real                            b0_mean_trt;            // mean of random intercept for treatment
  real                            b0_mean_ctrl;            // mean of random intercept for control
  real                            b1_mean_trt;            // mean of preslope for treatment
  real                            b1_mean_ctrl;            // mean of preslope for control
  real                            b2_mean_trt;            // mean of postslope for treatment
  real                            b2_mean_ctrl;            // mean of postslope for control
  vector<lower=0,upper=100>[3]    b_sd;               // SD of random effects assoc. with regression (? different for trt and ctrl)
  cholesky_factor_corr[4]         randeff_corr_chol;  // Cholesky decomposition of correlation of random effects (? this is a matrix)
  vector<lower=0,upper=1>[nobs_trt]   chgpt_obs_trt_raw;      // raw version of change point for trt random effect U(0, 1)
  vector<lower=0,upper=1>[nobs_ctrl]   chgpt_obs_ctrl_raw;      // raw version of change point random effect U(0, 1)
  vector<lower=0,upper=1>[ncen_trt]   chgpt_cen_trt_raw;      // raw version of change point random effect U(0, 1)
  vector<lower=0,upper=1>[ncen_ctrl]   chgpt_cen_ctrl_raw;      // raw version of change point random effect U(0, 1)
  matrix[3, nobs_trt]                 b_obs_trt_raw;          // raw version of random coefficients (will be std. normal)
  matrix[3, nobs_ctrl]                 b_obs_ctrl_raw;          // raw version of random coefficients (will be std. normal)
  matrix[3, ncen_trt]                 b_cen_trt_raw;          // raw version of random coefficients (will be std. normal)
  matrix[3, ncen_ctrl]                 b_cen_ctrl_raw;          // raw version of random coefficients (will be std. normal)
}
transformed parameters {
  vector[nobs_trt] chgpt_obs_trt;           // changepoint (observed times)
  vector[nobs_ctrl] chgpt_obs_ctrl;           // changepoint (observed times)
  vector[ncen_trt] chgpt_cen_trt;           // changepoint (censored times)
  vector[ncen_ctrl] chgpt_cen_ctrl;           // changepoint (censored times)
  matrix[3,nobs_trt] b_obs_trt;             // random intercept and slopes (observed event times)
  matrix[3,nobs_ctrl] b_obs_ctrl;             // random intercept and slopes (observed event times)
  matrix[3,ncen_trt] b_cen_trt;             // random intercept and slopes (censored event times)
  matrix[3,ncen_ctrl] b_cen_ctrl;             // random intercept and slopes (censored event times)
  matrix[4,4] randeff_cov_chol;     // cholesky factor of covariance matrix
  vector[3]   LHS;                  // sigma[1,2:4] / sigma[1,1] = chol[1,2:4] / chol[1,1]
  vector[3]   b_mean_trt;           // mean of random intercept for trt, preslope, and postslope
  vector[3]   b_mean_ctrl;          // mean of random intercept for ctrl, preslope, and postslope
  
  b_mean_trt[1] = b0_mean_trt;
  b_mean_trt[2] = b1_mean_trt;
  b_mean_trt[3] = b2_mean_trt;

  b_mean_ctrl[1] = b0_mean_ctrl;
  b_mean_ctrl[2] = b1_mean_ctrl;
  b_mean_ctrl[3] = b2_mean_ctrl;
  
  // Transform chgpt to get proper bounds
  chgpt_obs_trt = tobs_trt .* chgpt_obs_trt_raw;       // in [0, tobs]
  chgpt_obs_ctrl = tobs_ctrl .* chgpt_obs_ctrl_raw;       // in [0, tobs]
  chgpt_cen_trt = tcen_trtimp .* chgpt_cen_trt_raw;    // in [0, tcenimp] REQUIRES JACOBIAN ADJUSTMENT SINCE tecenimp IS RANDOM
  chgpt_cen_ctrl = tcen_ctrlimp .* chgpt_cen_ctrl_raw;    // in [0, tcenimp] REQUIRES JACOBIAN ADJUSTMENT SINCE tecenimp IS RANDOM
  
  // cholesky factor of covariance matrix of MVN prior to truncation
  randeff_cov_chol = diag_pre_multiply(append_row(chgpt_sd, b_sd), randeff_corr_chol);
  
  // Distribution of b_obs | chgpt_obs is MVN(condMean, condCov)
  LHS = randeff_cov_chol[2:4,1] * inv(randeff_cov_chol[1,1]);
  b_obs_trt   = randeff_cov_chol[2:4,2:4] * b_obs_trt_raw;             // implies b_obs | chgpt_obs ~ MVN(0, condCov)
  b_obs_ctrl   = randeff_cov_chol[2:4,2:4] * b_obs_ctrl_raw;             // implies b_obs | chgpt_obs ~ MVN(0, condCov)
  b_cen_trt   = randeff_cov_chol[2:4,2:4] * b_cen_trt_raw;             // implies b_obs | chgpt_obs ~ MVN(0, condCov)
  b_cen_ctrl   = randeff_cov_chol[2:4,2:4] * b_cen_ctrl_raw;             // implies b_obs | chgpt_obs ~ MVN(0, condCov)
  
  for ( i in 1:nobs_trt )
    b_obs_trt[, i] += b_mean_trt + LHS * (chgpt_obs_trt[i] - chgpt_mean_trt);  // implies b_obs | chgpt_obs ~ MVN(condMean, condCov)
  for ( i in 1:nobs_ctrl )
    b_obs_ctrl[, i] += b_mean_ctrl + LHS * (chgpt_obs_ctrl[i] - chgpt_mean_ctrl);  // implies b_obs | chgpt_obs ~ MVN(condMean, condCov)
  for ( i in 1:ncen_trt ) 
    b_cen_trt[, i] += b_mean_trt + LHS * (chgpt_cen_trt[i] - chgpt_mean_trt);  // implies b_cen | chgpt_cen ~ MVN(condMean, condCov)
  for ( i in 1:ncen_ctrl ) 
    b_cen_ctrl[, i] += b_mean_ctrl + LHS * (chgpt_cen_ctrl[i] - chgpt_mean_ctrl);  // implies b_cen | chgpt_cen ~ MVN(condMean, condCov)
}
model {
  // beta_y includes an indicator for trt
  vector[Nobs_trt] yobs_trt_mean = Xlong_obs_trt * beta_y;
  vector[Nobs_ctrl] yobs_ctrl_mean = Xlong_obs_ctrl * beta_y;
  vector[Ncen_trt] ycen_trt_mean = Xlong_cen_trt * beta_y;
  vector[Ncen_ctrl] ycen_ctrl_mean = Xlong_cen_ctrl * beta_y;
  real diff_ij;
  
  // Compute longitudinal means | random effects
  for ( i in 1:Nobs_trt ) {
    diff_ij = visitobs_trt[i] - chgpt_obs_trt[idobs_trt[i]];
    yobs_trt_mean[i] += b_obs_trt[1,idobs_trt[i]] + (diff_ij <= 0 ? b_obs_trt[2,idobs_trt[i]] : b_obs_trt[3,idobs_trt[i]] ) * diff_ij;
  }
  for ( i in 1:Nobs_ctrl ) {
    diff_ij = visitobs_ctrl[i] - chgpt_obs_ctrl[idobs_ctrl[i]];
    yobs_ctrl_mean[i] += b_obs_ctrl[1,idobs_ctrl[i]] + (diff_ij <= 0 ? b_obs_ctrl[2,idobs_ctrl[i]] : b_obs_ctrl[3,idobs_ctrl[i]] ) * diff_ij;
  }

  for ( i in 1:Ncen_trt ) {
    diff_ij = visitcen_trt[i] - chgpt_cen_trt[idcen_trt[i]];
    ycen_trt_mean[i] += b_cen_trt[1,idcen_trt[i]] + (diff_ij <= 0 ? b_cen_trt[2,idcen_trt[i]] : b_cen_trt[3,idcen_trt[i]] ) * diff_ij;
  }
  for ( i in 1:Ncen_ctrl ) {
    diff_ij = visitcen_ctrl[i] - chgpt_cen_ctrl[idcen_ctrl[i]];
    ycen_ctrl_mean[i] += b_cen_ctrl[1,idcen_ctrl[i]] + (diff_ij <= 0 ? b_cen_ctrl[2,idcen_ctrl[i]] : b_cen_ctrl[3,idcen_ctrl[i]] ) * diff_ij;
  }
  
  // Priors for TTE model
  beta_tte_trt  ~ std_normal();
  beta_tte_ctrl  ~ std_normal();
  shape_tte_trt ~ std_normal();
  shape_tte_ctrl ~ std_normal();
  scale_tte_trt ~ std_normal();
  scale_tte_ctrl ~ std_normal();
  
  // Priors for longitudinal model
  beta_y            ~ normal(0, 10);
  sd_y              ~ normal(0, 10);
  
  // Priors for random effects parameters
  randeff_corr_chol ~ lkj_corr_cholesky(1.0);   // implies corr ~ uniform over PD correlation matrices
  chgpt_mean_trt ~ generalized_normal(randeff_location[1], randeff_scale[1], randeff_shape[1]);
  chgpt_mean_ctrl ~ generalized_normal(randeff_location[1], randeff_scale[1], randeff_shape[1]);
  b0_mean_trt    ~ generalized_normal(randeff_location[2], randeff_scale[2], randeff_shape[2]);
  b0_mean_ctrl    ~ generalized_normal(randeff_location[2], randeff_scale[2], randeff_shape[2]);
  b1_mean_trt    ~ generalized_normal(randeff_location[3], randeff_scale[3], randeff_shape[3]);
  b1_mean_ctrl    ~ generalized_normal(randeff_location[3], randeff_scale[3], randeff_shape[3]);
  b2_mean_trt    ~ generalized_normal(randeff_location[4], randeff_scale[4], randeff_shape[4]);
  b2_mean_ctrl    ~ generalized_normal(randeff_location[4], randeff_scale[4], randeff_shape[4]);
  chgpt_sd   ~ std_normal();
  b_sd       ~ std_normal();
  
  // Complete data likelihood (observed)
  tobs_trt                 ~ weibull_ph(Xtte_obs_trt, beta_tte_trt, shape_tte_trt, scale_tte_trt);
  tobs_ctrl                 ~ weibull_ph(Xtte_obs_ctrl, beta_tte_ctrl, shape_tte_ctrl, scale_tte_ctrl);
  yobs_trt                 ~ normal(yobs_trt_mean, sd_y);  // conditional distribution of yobs | random effects
  yobs_ctrl                 ~ normal(yobs_ctrl_mean, sd_y);  // conditional distribution of yobs | random effects
  chgpt_obs_trt            ~ tnormal(chgpt_mean_trt, randeff_cov_chol[1,1], tobs_trt);
  chgpt_obs_ctrl            ~ tnormal(chgpt_mean_ctrl, randeff_cov_chol[1,1], tobs_ctrl);
  to_vector(b_obs_trt_raw) ~ std_normal();
  to_vector(b_obs_ctrl_raw) ~ std_normal();
  
  // Complete data likelihood (censored)
  tcen_trtimp              ~ weibull_ph(Xtte_cen_trt, beta_tte_trt, shape_tte_trt, scale_tte_trt);
  tcen_ctrlimp              ~ weibull_ph(Xtte_cen_ctrl, beta_tte_ctrl, shape_tte_ctrl, scale_tte_ctrl);
  ycen_trt                 ~ normal(ycen_trt_mean, sd_y);                                // conditional distribution of ycen | random effects
  ycen_ctrl                 ~ normal(ycen_ctrl_mean, sd_y);                                // conditional distribution of ycen | random effects
  chgpt_cen_trt            ~ tnormal(chgpt_mean_trt, randeff_cov_chol[1,1], tcen_trtimp);
  chgpt_cen_ctrl            ~ tnormal(chgpt_mean_ctrl, randeff_cov_chol[1,1], tcen_ctrlimp);
  to_vector(b_cen_trt_raw) ~ std_normal();
  to_vector(b_cen_ctrl_raw) ~ std_normal();
  target += log(tcen_trtimp); // jacobian adjustment for change point
  target += log(tcen_ctrlimp); // jacobian adjustment for change point
  
}
generated quantities {
  corr_matrix[4] randeff_corr = multiply_lower_tri_self_transpose(randeff_corr_chol);
}


