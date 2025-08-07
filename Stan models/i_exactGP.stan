functions {
  // base covariance fn
  
  // SE
    matrix se(vector x, real alpha_obs, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha_obs = pow(alpha_obs, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha_obs + delta;
      for(j in (i + 1):N) {
        K[i, j] = exp(r * square(x[i] - x[j])) * sq_alpha_obs;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
}

data {
  int<lower=1> N;				// sample size
	int<lower=1> D;       // output dims
	vector[N] inputs;				//matrix of total (training and test) observations
	real x_min;      // lower bound for true_x prior
	real<lower=x_min> x_max;      // upper bound for true_x prior
	matrix[N, D] y_f; // for output 1
	matrix[N, D] y_g; // for output 2
	real<lower=0> s; // measurement sd for latent inputs
	// User input for which covariance function to use
	int<lower=0, upper=1> latent; // indicator if x is latent. 1: latent inputs
	int<lower=0, upper=1> is_vary; // 0 = constant hyperparams for each dims; 1 = varying params
  int<lower=0, upper=1> is_corr; // 0 = no correlation b/w dims; 1 = correlated dims
  int<lower=0, upper=1> rho_prior; // 0 = Normal; 1 = InvGamma
  // prior hyperparams specification (only for two parameter dist families)
  vector[2] ls_param_f; 
  vector[2] ls_param_g; 
  vector[2] msd_param_f;
  vector[2] msd_param_g;
  vector[2] esd_param_f;
  vector[2] esd_param_g;
  // Input data specific mean and sd for the y1 and y2
  vector[2] intc_yf;
  vector[2] intc_yg;
  // covfn diag constant
  //real nugget;
}

transformed data {
// add constant to the diagonal of the covariance matrix for computational stability
  real delta = 1e-6;
}

parameters {
  // GP Length scale parameter (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] rho_temp_f;
  // GP marginal SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] alpha_temp_f;
  // Error SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] sigma_temp_f;
  // GP Length scale parameter (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] rho_temp_g;
  // GP marginal SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] alpha_temp_g;
  // Error SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] sigma_temp_g;
  // only for latent variable models
  vector<lower=x_min,upper=x_max>[latent==1 ? N:0] x_temp;
  //vector[latent==1 ? N:0] z;
  // Between dimension correlation parameter (temporary storage)
  cholesky_factor_corr[is_corr==1 ? D:0] L_omega_temp_f;
  cholesky_factor_corr[is_corr==1 ? D:0] L_omega_temp_g;
	array[D] real intercept_yf;
	array[D] real intercept_yg;
	// Gaussian link parameters
	matrix[N, D] eta_f;
	matrix[N, D] eta_g;
}

transformed parameters {
  // Model condition adjustments (temporary to actual parameters)
  vector<lower=0>[D] rho_f;
  vector<lower=0>[D] alpha_f;
  vector<lower=0>[D] sigma_f;
  vector<lower=0>[D] rho_g;
  vector<lower=0>[D] alpha_g;
  vector<lower=0>[D] sigma_g;
	vector[N] x;
	matrix[N, D] f;
	matrix[N, D] g;
	cholesky_factor_corr[D] L_omega_f;
	cholesky_factor_corr[D] L_omega_g;
	 if (latent) {
	  // x = inputs + z * s;
	  x = x_temp;
	} else{
	  x = inputs;
	}
	// if params are constant for dims, it will be repeated
  if (is_vary) {
    rho_f = rho_temp_f;
    alpha_f = alpha_temp_f;
    sigma_f = sigma_temp_f;
    rho_g = rho_temp_g;
    alpha_g = alpha_temp_g;
    sigma_g = sigma_temp_g;
  } else {
    for (k in 1:D) {
      rho_f[k] = rho_temp_f[1];
      alpha_f[k] = alpha_temp_f[1];
      sigma_f[k] = sigma_temp_f[1];
      rho_g[k] = rho_temp_g[1];
      alpha_g[k] = alpha_temp_g[1];
      sigma_g[k] = sigma_temp_g[1];
    }
  }
  if (is_corr==0) {
    L_omega_f = diag_matrix(rep_vector(1, D));
    L_omega_g = diag_matrix(rep_vector(1, D));
  } else {
    L_omega_f = L_omega_temp_f;
    L_omega_g = L_omega_temp_g;
  }
  // Computing covariance matrix for standard GP
      for (k in 1:D) {
        matrix[N, N] K_f = se(x, alpha_f[k], rho_f[k], delta);
        matrix[N, N] K_g = se(x, alpha_g[k], rho_g[k], delta);
        f[, k] = K_f * eta_f[, k];
        g[, k] = K_g * eta_g[, k];
    }
  // For correlated outputs
  f = f * L_omega_f';
  g = g * L_omega_g';
}


model {
  	// Priors
  if (rho_prior == 0) {
    rho_temp_f ~ normal(ls_param_f[1], ls_param_f[2]); 
    rho_temp_g ~ normal(ls_param_g[1], ls_param_g[2]); 
  } else if (rho_prior == 1) {
    rho_temp_f ~ inv_gamma(ls_param_f[1], ls_param_f[2]);
    rho_temp_g ~ inv_gamma(ls_param_g[1], ls_param_g[2]);
  }
  alpha_temp_f ~ normal(msd_param_f[1], msd_param_f[2]);
  alpha_temp_g ~ normal(msd_param_g[1], msd_param_g[2]);
  sigma_temp_f ~ normal(esd_param_f[1], esd_param_f[2]);
  sigma_temp_g ~ normal(esd_param_g[1], esd_param_g[2]);
  L_omega_temp_f ~ lkj_corr_cholesky(1);
  L_omega_temp_g ~ lkj_corr_cholesky(1);
  if (latent) {
    // z ~ std_normal();
    inputs ~ normal(x_temp, s);
	}
	// set prior on x to match data generating process 
	x_temp ~ uniform(x_min, x_max);
	for(j in 1:D) {
	  intercept_yf[j] ~ normal(intc_yf[1], intc_yf[2]);
	  intercept_yg[j] ~ normal(intc_yg[1], intc_yg[2]);
	  to_vector(eta_f[, j]) ~ std_normal();
	  to_vector(eta_g[, j]) ~ std_normal();
	}
	// Likelihood
	for(j in 1:D) {
	  y_f[,j] ~ normal(intercept_yf[j] + f[,j], sigma_f[j]);  
	  y_g[,j] ~ normal(intercept_yg[j] + g[,j], sigma_g[j]); 
	}
}
