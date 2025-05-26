#library
library(rstan)
library(posterior)
library(ggplot2)
library(bayesplot)

# Source fns 
source('indpcompgpfns.R')
comphsgpmodel <- stan_model('indpcomphsgp_maternclass.stan')
set.seed(123)
N <- 50  # higher N will take very long. Better to stick to <=30
dims <- 10
# Set true input x
x_min <- 0
x_max <- 10
#x_true <- seq(0.5, by = 0.5, length.out = N) #length.out = N)
s <- 0.3 # latent sd # needed in case of latent GPs
# Generate parameters
marginal_sd_params_f <- c(3, 0.25)
error_sd_params_f <- c(1, 0.25)
ls_params_f <- c(1, 0.05)
marginal_sd_params_g <- c(2, 0.25)
error_sd_params_g <- c(1, 0.25)
ls_params_g <- c(1.5, 0.05)
intc_params_y1 <- c(0, 1)
intc_params_y2 <- c(0, 5)
delta <- 1e-12
covfn <- 'se' # only needed for gp data

# Model specifications
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
marginal_sd_params_model_f <- marginal_sd_params_f
error_sd_params_model_f <- error_sd_params_f 
ls_params_model_f <- ls_params_f 
marginal_sd_params_model_g <- marginal_sd_params_g
error_sd_params_model_g <- error_sd_params_g
ls_params_model_g <- ls_params_g
intc_params_model_y1 <- intc_params_y1
intc_params_model_y2 <- intc_params_y2
latent_model <- 1 # 1 for latent inputs, 0 for manifest
covfn_model <- 0  # 0 - SE; 1 = Matern3/2; 2 = Matern 5/2
adapt_delta_model <- 0.95

params1 <- true_params_vary(n_obs = N, 
                               dims = dims, 
                               msd_pars = marginal_sd_params_f, 
                               esd_pars = error_sd_params_f, 
                               ls_pars = ls_params_f, 
                               intc_pars = intc_params_y1) 
params2 <- true_params_vary(n_obs = N, 
                             dims = dims, 
                             msd_pars = marginal_sd_params_g, 
                             esd_pars = error_sd_params_g, 
                             ls_pars = ls_params_g, 
                             intc_pars = intc_params_y2) 

# Simulate data
dat1 <- gp_sim_data(n_obs = N, 
                     dims = dims,
                     true_x = x_true,
                     covfn = covfn,
                     delta = delta,
                     alpha = params1$alpha, 
                     sigma = params1$sigma,
                     rho = params1$rho, 
                   intc = params1$intc)
dat2 <- gp_sim_data(n_obs = N, 
                    dims = dims, 
                    true_x = x_true,
                    covfn = covfn,
                    delta = delta,
                    alpha = params2$alpha, 
                    sigma = params2$sigma,
                    rho = params2$rho, 
                    intc = params2$intc)
y1 <- dat1[1:N,1:dims]
f <- dat1[1:N,(dims+1):(2*dims)]
y2 <- dat2[1:N,1:dims]
g <- dat2[1:N,(dims+1):(2*dims)]
x_true <- runif(N, x_min, x_max)
input_dat <- sim_input(n_obs = N, true_x = x_true, s_x = s)
x_true <- input_dat$x_true
x_obs <- input_dat$x_obs
# Check generated GP
dim_id <- 1
fplot <- gp_plot(dim_id = dim_id, output = y1, gpfns = f, input = x_true, 
          msd = params1$alpha, esd = params1$sigma, 
          plot_type = 'gpfns', label = 'f')
gplot <- gp_plot(dim_id = dim_id, output = y1, gpfns = g, input = x_true, 
       msd = params2$alpha, esd = params2$sigma, 
       plot_type = 'gpfns', label = 'g')
bayesplot_grid(fplot, gplot)
# fit hsgp model
min_m_f <- ceiling(1.75 * (5/4) * (max(x_true) - min(x_true))/(2 * ls_params_f[1])) # mean(rho) = 0.2, S = 0.5 (half range of x), c = 5/4
min_m_g <- ceiling(1.75 * (5/4) * (max(x_true) - min(x_true))/(2 * ls_params_g[1])) # mean(rho) = 0.2, S = 0.5 (half range of x), c = 5/4
hsgp_fit <- hsgp_model(comphsgpmodel,
                   n_obs = N, 
                   dims = dims,
                   m_obs_f = min_m_f,
                   m_obs_g = min_m_g,
                   x_min = x_min,
                   x_max = x_max,
                   c = 5/4,
                   outputs1 = y1, 
                   outputs2 = y2, 
                   inputs = x_obs, 
                   latent_sd = s, 
                   latent_inputs = latent_model, 
                   covfn = covfn_model, 
                   rho_prior = rho_prior_model,
                   ls_param_f = ls_params_model_f,
                   ls_param_g = ls_params_model_g, 
                   msd_param_f = marginal_sd_params_model_f,
                   msd_param_g = marginal_sd_params_model_g,
                   esd_param_f = error_sd_params_model_f,
                   esd_param_g = error_sd_params_model_g,
                   intc_y1 = intc_params_model_y1,
                   intc_y2 = intc_params_model_y2,
                   is_vary = 1, 
                   is_corr = 1, 
                   iter = 2000, 
                   warmup = 1000, 
                   chains = 2, 
                   cores = 2, 
                   init = 0, 
                   adapt_delta = adapt_delta_model)
# model output summary
print(hsgp_fit, pars = c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))
traceplot(hsgp_fit, pars = c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))
#### post summary
x_names <- sprintf('x[%s]', seq(1:N))
out <- as_draws_matrix(hsgp_fit)
subset <- subset_draws(out, variable = x_names)
summary <- summarise_draws(subset)
sd <- rep(NA, length(x_names))
abs_bias <- rep(NA, length(x_names))
rmse <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd[i] <- sd(subset[,i])
  abs_bias[i] <- abs_bias_draws(subset[,i], x_true[i])
  rmse[i] <- rmse_draws(subset[,i], x_true[i])
}
mean(rmse)
mean(abs_bias)
mean(sd)
