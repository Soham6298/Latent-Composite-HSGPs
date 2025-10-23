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
	real spd_se_deriv(real alpha, real rho, real w) {
	  real S;
    real sq_alpha = alpha^2;
    real sq_rho = rho^2;
    S = sqrt(2*pi()) * sq_alpha * rho * (w^2) * exp(-0.5*(sq_rho)*(w^2));
    return S;
	}
}

data {
	real L;						//boundary condition factor
	int<lower=1> M;				// no. of basis functions		
	int<lower=1> N;				// sample size
	int<lower=1> D;       // output dims
	vector[N] inputs;				//matrix of total (training and test) observations
	real x_min;      // lower bound for true_x prior
	real x_max;      // upper bound for true_x prior
	matrix[N, D] y;
	real<lower=0> s;
	// User input for which covariance function to use
	int<lower=0, upper=3> covfn;  //0: se, 1: m32, 2: m52; 3: se_deriv
	int<lower=0, upper=1> latent; // indicator if x is latent. 1: latent inputs
	int<lower=0, upper=1> is_vary; // 0 = constant hyperparams for each dims; 1 = varying params
  int<lower=0, upper=1> is_corr; // 0 = no correlation b/w dims; 1 = correlated dims
  int<lower=0, upper=1> rho_prior; // 0 = Normal; 1 = InvGamma
  // prior hyperparams specification (only for two parameter dist families)
  vector[2] ls_param; 
  vector[2] msd_param;
  vector[2] esd_param;
  // Input data specific mean and sd for the y
  vector[2] intc_y;
}

parameters {
	array[D] vector[M] beta;
	// GP Length scale parameter (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] rho_temp;
  // GP marginal SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] alpha_temp;
  // Error SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] sigma_temp;
  // only for latent variable models
  vector<lower=x_min,upper=x_max>[latent==1 ? N:0] x_temp;
  // Between dimension correlation parameter (temporary storage)
  cholesky_factor_corr[is_corr==1 ? D:0] L_omega_temp;
	array[D] real intercept;
}

transformed parameters{
	// Model condition adjustments (temporary to actual parameters)
  vector<lower=0>[D] rho;
  vector<lower=0>[D] alpha;
  vector<lower=0>[D] sigma;
	vector[N] x;
	matrix[N, D] f;
	vector[M] diagSPD;
	matrix[M, D] SPD_beta;
	matrix[N,M] PHI;
	cholesky_factor_corr[D] L_omega;
	 if (latent) {
	  x = x_temp;
	} else{
	  x = inputs;
	}
	// if params are constant for dims, it will be repeated
  if (is_vary) {
    rho = rho_temp;
    alpha = alpha_temp;
    sigma = sigma_temp;
  } else {
    for (k in 1:D) {
      rho[k] = rho_temp[1];
      alpha[k] = alpha_temp[1];
      sigma[k] = sigma_temp[1];
    }
  }
  if (is_corr==0) {
    L_omega = diag_matrix(rep_vector(1, D));
  } else {
    L_omega = L_omega_temp;
  }
	for (m in 1:M){ 
	  PHI[,m] = phi(L, m, x); 
	  }
	for(j in 1:D){
	 for(m in 1:M){ 
	  if(covfn==1){
	   diagSPD[m] =  sqrt(spd_m32(alpha[j], rho[j], sqrt(lambda(L, m)))); 
	  } else if(covfn==2){
	   diagSPD[m] =  sqrt(spd_m52(alpha[j], rho[j], sqrt(lambda(L, m)))); 
	  } else if(covfn==3){
	   diagSPD[m] =  sqrt(spd_se_deriv(alpha[j], rho[j], sqrt(lambda(L, m)))); 
	  } else {
	    diagSPD[m] = sqrt(spd_se(alpha[j], rho[j], sqrt(lambda(L, m)))); 
	  }
	}
	  SPD_beta[,j] = diagSPD .* beta[j];
	  f[,j] = PHI * SPD_beta[,j]; 
	}
	// For correlated outputs
    f = f * L_omega';
}

model{
  for (j in 1: D) {
    beta[j] ~ normal(0,1);
  }
	// Priors
  if (rho_prior == 0) {
    rho_temp ~ normal(ls_param[1], ls_param[2]); 
  } else if (rho_prior == 1) {
    rho_temp ~ inv_gamma(ls_param[1], ls_param[2]);
  }
  alpha_temp ~ normal(msd_param[1], msd_param[2]);
  sigma_temp ~ normal(esd_param[1], esd_param[2]);
  L_omega_temp ~ lkj_corr_cholesky(1);
  if (latent) {
    inputs ~ normal(x_temp, s);
	}
	// set prior on x to match data generating process 
	x_temp ~ uniform(x_min, x_max);
	for(j in 1:D) {
	  intercept[j] ~ normal(intc_y[1], intc_y[2]);
	}
	// Likelihood
	for(j in 1:D) {
	 y[,j] ~ normal(intercept[j] + f[,j], sigma[j]);  
	}
}
