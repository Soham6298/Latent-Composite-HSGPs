functions {
  // base covariance fn
  
  // SE
    matrix se(vector x, real alpha, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha = pow(alpha, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha + delta;
      for(j in (i + 1):N) {
        K[i, j] = exp(r * square(x[i] - x[j])) * sq_alpha;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  // Matern 3/2
  matrix m32(vector x, real alpha, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_alpha = pow(alpha, 2);
    real r = -inv(rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha + delta;
      for(j in (i + 1):N) {
        K[i, j] = (1 - (r * sqrt(3) * abs(x[i] - x[j]))) * 
                    exp(r * sqrt(3) * abs(x[i] - x[j])) * sq_alpha;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  //Matern 5/2 
   matrix m52(vector x, real alpha, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha = pow(alpha, 2);
    real r = -inv(2 * sq_rho);
    
    for(i in 1:N) {
      K[i, i] = sq_alpha + delta;
      for(j in (i + 1):N) {
        K[i, j] = (1 + (sqrt(5) * abs(x[i] - x[j])/rho) + 
                    ((5 * square(x[i] - x[j]))/ (3 * sq_rho))) * 
                    exp(- sqrt(5) * abs(x[i] - x[j])/rho) * sq_alpha;
        K[j, i] = K[i, j];
      }
    }
    return cholesky_decompose(K);
  }
  
  // Deriv SE
  matrix se_deriv(vector x, real alpha, real rho, real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_rho = square(rho);
    real rho4 = pow(rho, 4);
    real sq_alpha = pow(alpha, 2);
    real r = -inv(2 * sq_rho);
    
    for (i in 1:N) {
      K[i, i] = (sq_alpha / sq_rho) + delta;
      for (j in (i + 1):N) {
        K[i, j] = exp(r * square(x[i] - x[j])) * 
            (sq_rho - square(x[i] - x[j])) * (sq_alpha / rho4);
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
	matrix[N, D] y; // response
	real<lower=0> s; // measurement sd for latent inputs
	// User input for which covariance function to use
	int<lower=0, upper=3> covfn;  // 0: se, 1: m32, 2: m52; 3: deriv (SE)
	int<lower=0, upper=1> latent; // indicator if x is latent. 1: latent inputs
	int<lower=0, upper=1> is_vary; // 0 = constant hyperparams for each dims; 1 = varying params
  int<lower=0, upper=1> is_corr; // 0 = no correlation b/w dims; 1 = correlated dims
  int<lower=0, upper=1> rho_prior; // 0 = Normal; 1 = InvGamma
  // prior hyperparams specification (only for two parameter dist families)
  vector[2] ls_param;
  vector[2] msd_param;
  vector[2] esd_param;
  // Input data specific mean and sd for the y1 and y2
  vector[2] intc;
  real delta;
}

transformed data {
// add constant to the diagonal of the covariance matrix for computational stability
  //real delta = 1e-10;
}

parameters {
	// Gaussian link parameters
	matrix[N,D] eta;
	// GP Length scale parameter (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] rho_temp;
  // GP marginal SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] alpha_temp;
  // Error SD parameters (temporary storage)
  vector<lower=0>[is_vary==1 ? D:1] sigma_temp;
  // GP Length scale parameter (temporary storage)
  // only for latent variable models
  vector<lower=x_min,upper=x_max>[latent==1 ? N:0] x_temp;
  //vector[latent==1 ? N:0] z;
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
	cholesky_factor_corr[D] L_omega;
	 if (latent) {
	  // x = inputs + z * s;
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
  // Computing covariance matrix for standard GP
   if (covfn == 1) {
    for (j in 1:D) {
      matrix[N, N] K = m32(x, alpha[j], rho[j], delta);
      f[, j] = K * eta[, j];
    } 
  } else if (covfn == 2) {
    for (j in 1:D) {
      matrix[N, N] K = m52(x, alpha[j], rho[j], delta);
      f[, j] = K * eta[, j];
    } 
  } else if (covfn == 3) {
    for (j in 1:D) {
     matrix[N, N] K = se_deriv(x, alpha[j], rho[j], delta);
      f[, j] = K * eta[, j];
      // print(eta);
      //print(eta);
     }
    } else {
    for (j in 1:D) {
      matrix[N, N] K = se(x, alpha[j], rho[j], delta);
      f[,j] = K * eta[,j];
    }
  }
    // For correlated outputs
  f = f * L_omega';
}

model{
  for (j in 1: D) {
    to_vector(eta[,j]) ~ std_normal();
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
    // z ~ std_normal();
    inputs ~ normal(x_temp, s);
	}
	// set prior on x to match data generating process 
	x_temp ~ uniform(x_min, x_max);
	for(j in 1:D) {
	  intercept[j] ~ normal(intc[1], intc[2]);
	}
	// Likelihood
	for(j in 1:D) {
	  y[,j] ~ normal(intercept[j] + f[,j], sigma[j]);  
	}
}
