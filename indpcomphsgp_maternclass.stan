functions {
	real lambda(real L, int m) {
		real lam;
		lam = ((m*pi())/(2*L))^2;
				
		return lam;
	}
	vector phi(real L, int m, vector x) {
		vector[rows(x)] fi;
		fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));
				
		return fi;
	}
	// Spectral densities
	real spd_se(real alpha, real rho, real w) {
		real S;
		S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));
				
		return S;
	}
	real spd_m32(real alpha, real rho, real w) {
		real S;
		S = (alpha^2) * ((2 * tgamma(2) * 3^(1.5)) / (0.5 * rho^3)) * ((3/rho^2) + (w^2))^(-2);
		return S;
	}
	real spd_m52(real alpha, real rho, real w) {
		real S;
		S = (alpha^2) * ((2 * tgamma(3) * 5^(2.5)) / (0.75 * rho^5)) * ((5/rho^2) + (w^2))^(-3);
		return S;
	}
}

data {
	real L;						//boundary condition factor
	int<lower=1> M_f;				// no. of basis functions	for f	
	int<lower=1> M_g;				// no. of basis functions	for g	
	int<lower=1> N;				// sample size
	int<lower=1> D;       // output dims
	vector[N] inputs;				//matrix of total (training and test) observations
	real x_min;      // lower bound for true_x prior
	real<lower=x_min> x_max;      // upper bound for true_x prior
	matrix[N, D] y1; // for output 1
	matrix[N, D] y2; // for output 2
	real<lower=0> s; // measurement sd for latent inputs
	// User input for which covariance function to use
	int<lower=0, upper=2> covfn;  //0: se, 1: m32, 2: m52
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
  vector[2] intc_y1;
  vector[2] intc_y2;
}

parameters {
	array[D] vector[M_f] beta_f;
	array[D] vector[M_g] beta_g;
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
  cholesky_factor_corr[is_corr==1 ? D:0] L_omega_temp;
	array[D] real intercept_y1;
	array[D] real intercept_y2;
}

transformed parameters{
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
	vector[M_f] diagSPD_f;
	vector[M_g] diagSPD_g;
	matrix[M_f, D] SPD_beta_f;
	matrix[M_g, D] SPD_beta_g;
	matrix[N,M_f] PHI_f;
	matrix[N,M_g] PHI_g;
	cholesky_factor_corr[D] L_omega;
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
    L_omega = diag_matrix(rep_vector(1, D));
  } else {
    L_omega = L_omega_temp;
  }
	for (m in 1:M_f){
	  PHI_f[, m] = phi(L, m, x);
	}
	for (m in 1:M_g){
	  PHI_g[, m] = phi(L, m, x);
	}
	for (j in 1:D) {
	  for (m in 1:M_f) { 
	    if (covfn==1) {
	      diagSPD_f[m] = sqrt(spd_m32(alpha_f[j], rho_f[j], sqrt(lambda(L, m))));
	    } else if (covfn==2) {
	      diagSPD_f[m] = sqrt(spd_m52(alpha_f[j], rho_f[j], sqrt(lambda(L, m))));
	    } else {
	      diagSPD_f[m] = sqrt(spd_se(alpha_f[j], rho_f[j], sqrt(lambda(L, m)))); 
	    }
	  }
	  for (m in 1:M_g) {
	    // TODO: compute this in transformed data
	    real sqrt_lambda_m = sqrt(lambda(L, m));
	    if (covfn==1) {
	      diagSPD_g[m] = sqrt(spd_m32(alpha_g[j], rho_g[j], sqrt_lambda_m)); 
	    } else if (covfn==2) {
	      diagSPD_g[m] = sqrt(spd_m52(alpha_g[j], rho_g[j], sqrt_lambda_m)); 
	    } else {
	      diagSPD_g[m] = sqrt(spd_se(alpha_g[j], rho_g[j], sqrt_lambda_m)); 
	    }
	  }
	  SPD_beta_f[,j] = diagSPD_f .* beta_f[j];
	  SPD_beta_g[,j] = diagSPD_g .* beta_g[j];
	  f[,j] = PHI_f * SPD_beta_f[,j]; 
	  g[,j] = PHI_g * SPD_beta_g[,j]; 
	}
	// For correlated outputs
  f = f * L_omega';
  g = g * L_omega';
}

model{
  for (j in 1: D) {
    beta_f[j] ~ normal(0,1);
    beta_g[j] ~ normal(0,1);
  }
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
  L_omega_temp ~ lkj_corr_cholesky(1);
  if (latent) {
    // z ~ std_normal();
    inputs ~ normal(x_temp, s);
	}
	// set prior on x to match data generating process 
	x_temp ~ uniform(x_min, x_max);
	for(j in 1:D) {
	  intercept_y1[j] ~ normal(intc_y1[1], intc_y1[2]);
	  intercept_y2[j] ~ normal(intc_y2[1], intc_y2[2]);
	}
	// Likelihood
	for(j in 1:D) {
	  y1[,j] ~ normal(intercept_y1[j] + f[,j], sigma_f[j]);  
	  y2[,j] ~ normal(intercept_y2[j] + g[,j], sigma_g[j]); 
	}
}
