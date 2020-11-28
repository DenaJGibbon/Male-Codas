data {
  int<lower=1> K; // number of features
  int<lower=1> J; // number of unique groups (for random intercepts)
  int<lower=1> M; // number of unique sites
  int<lower=0> N; // number of observations (sample size) 
  vector[K] y[N]; // N * K matrix of observations
  int<lower=1> group[N]; // vector of integer-coded group IDs for each observation
  int<lower=1> site[N]; // vector of integer-coded site IDs for each observation
}

transformed data {
  vector[K] zeros;
  zeros = rep_vector(0, K);
}

parameters {
  cholesky_factor_corr[K] L_Omega; // for residual vectors
  vector<lower=0>[K] L_sigma; // for residual vectors
  
  cholesky_factor_corr[K] L_Phi; // for random intercepts <- for groups
  vector<lower=0>[K] L_theta; // for random intercepts <- for groups
  vector[K] mu_group[J]; // J * K matrix of random intercepts <- for groups
  
  cholesky_factor_corr[K] L_Zeta; // for random intercepts <- for sites
  vector<lower=0>[K] L_beta; // for random intercepts <- for sites
  vector[K] mu_site[M]; // M * K matrix of random intercepts <- for site
  
  real<lower=2>nu_site; // Degrees of freedom parameter for multivariate-t
  real<lower=2>nu_group; // Degrees of freedom parameter for multivariate-t
  real<lower=2>nu_obs; // Degrees of freedom parameter for multivariate-t
  
}


transformed parameters {
  // Definitions and calculations are here, rather than in the model block, so that 
  // these are available for generated quantities. 
  matrix[K, K] L_Sigma; // for residual var/cov matrix
  matrix[K, K] L_Theta; // for female var/cov matrix
  matrix[K, K] L_Beta; // for site var/cov matrix
  vector[K] mu[N]; // N * K matrix of predicted means
  vector[K] resids[N]; // N * K matrix of raw residuals
  
  matrix[K, K] Beta;
  matrix[K, K] Theta;
  matrix[K, K] Sigma;
  
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  L_Theta = diag_pre_multiply(L_theta, L_Phi);
  L_Beta = diag_pre_multiply(L_beta, L_Zeta);
  for (n in 1:N){
    // Predicted mean contains nested reference for the group intercept, 
    // along with nested reference for site. 
    mu[n] = mu_group[group[n]] + mu_site[site[n]];
    resids[n] = y[n] - mu[n];
  }
  
  // All covariances are defined here.
  Beta = tcrossprod(L_Beta); // Var-cov matrix for site
  Theta = tcrossprod(L_Theta); // Female
  Sigma = tcrossprod(L_Sigma); // Observation
  
}

model {
  
  // residual variance components
  // For lkj: prior=1 means uniform prior, 
  // if prior=2 then correlations off diagonal close to zero.
  L_Omega ~ lkj_corr_cholesky(1.5); 
  L_sigma ~ cauchy(0, 5);
  
  // variance components for group
  L_Phi ~ lkj_corr_cholesky(1.5);
  L_theta ~ cauchy(0, 5);
  
  // variance components for site
  L_Zeta ~ lkj_corr_cholesky(1.5);
  L_beta ~ cauchy(0, 5);
  
  mu_site ~ multi_student_t(nu_site, zeros, Beta); // Site contributions. 
  mu_group  ~ multi_student_t(nu_group, zeros, Theta); // Group (female) contributions. 
  y ~ multi_student_t(nu_obs, mu, Sigma); // Observation-level model.  

  // t-distribution degrees of freedom
  nu_site ~ gamma(2, 0.1); // Recommended by Aki on Gelman's blog, 17 May 2015. 
  nu_group ~ gamma(2, 0.1); 
  nu_obs ~ gamma(2, 0.1); 
  
}

generated quantities {
  matrix[K, K] Vcov_obs;
  matrix[K, K] Vcov_group;
  matrix[K, K] Vcov_site;
  matrix[K, K] L_Sigma_inv; 
  matrix[K, K] Sigma_inv;

  vector[K] ICC_group; 
  vector[K] ICC_site;
  vector[N] Maha_sqd;
  vector[K] site_rand_intercept[M];

  real DF_site; 
  real DF_group;
  real DF_obs;  
  real IF_site; 
  real IF_group;
  real IF_obs;
  
  DF_site = nu_site;
  // Inflation factor that multiplies the scale to calculate a variance.  
  IF_site = nu_site/(nu_site-2); 
  DF_group = nu_group;
  IF_group= nu_group/(nu_group-2);
  DF_obs = nu_obs;
  IF_obs= nu_obs/ (nu_obs-2);
  
  // Variance/covariance matrices, from scale matrices of the multivariate t. 
  Vcov_obs= IF_obs * Sigma;
  Vcov_group= IF_group * Theta;
  Vcov_site= IF_site * Beta;
  
  site_rand_intercept = mu_site;

  //Loop over features to get ICCs for group and site.  
  for (k in 1:K){  
    ICC_group[k] = Vcov_group[k, k] / (Vcov_obs[k, k] + Vcov_site[k, k] + Vcov_group[k, k]);   
    ICC_site[k] = Vcov_site[k, k] / (Vcov_obs[k, k] + Vcov_site[k, k] + Vcov_group[k, k]);
  }

  // Mahalanobis distances for goodness-of-fit check. 
  // Inverse Cholesky factor of Sigma.  
  L_Sigma_inv = inverse(L_Sigma);
  // Inverse covariance matrix from inverted Cholesky factor.
  Sigma_inv = crossprod(L_Sigma_inv); 
  // Squared Mahalanobis distance of each observation from its predicted value.  
  for (n in 1:N) Maha_sqd[n] = quad_form_sym(Sigma_inv, resids[n,]); 
  
}
