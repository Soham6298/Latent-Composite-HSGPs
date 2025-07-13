#library
library(rstan)
library(posterior)
library(ggplot2)
library(bayesplot)

# Source fns 
source('indcompgpfns.R')
testderivhsgpmodel <- stan_model('testderivhsgp.stan')
testderivgpmodel <- stan_model('testderivgp.stan')
set.seed(2414)
## Specify data conditions
N <- 20  # higher N will take very long. Better to stick to <=30
dims <- 5
# Set true input x
#x_true <- seq(0.5, 9.5, length.out = N)
x_min <- 0
x_max <- 10
s_x <- 0.3 # latent sd # needed in case of latent GPs
# Generate parameters
marginal_sd_params <- c(3, 0.25)
error_sd_params <- c(1, 0.25)
ls_params <- c(1, 0.05)
intc_params <- c(0, 1)
delta <- 1e-12
# Choices of covfn: 'se';'m32';'m52';'se_deriv'
covfn <- 'se_deriv'

# Model specifications
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
marginal_sd_params_model <- marginal_sd_params # alpha_grad ~ N( , )
error_sd_params_model <- error_sd_params    # sigma_grad ~ N( , )
ls_params_model <- ls_params # rho ~ N( , ) or InvGamma( , )
#ls_params_model_deriv <- c(0.7, 0.05)
intc_params_model <- intc_params
# Set covariance function based scaler for HSGPs
# if(covfn_model==1){
#   const <- 3.42
# }else if(covfn_model==2){
#   const <- 2.65
# }else {
#   const <- 1.75
# }
#min_m <- 33 #ceiling(const * (5/4) * (x_max - x_min)/(2*ls_params_model[1])) 
adapt_delta_model <- 0.95
x_true <- runif(N, x_min, x_max)
params <- true_params_vary(n_obs = N, 
                           dims = dims,
                           msd_pars = marginal_sd_params,
                           esd_pars = error_sd_params,
                           ls_pars = ls_params,
                           intc_pars = intc_params)
data <- test_gp_sim_data(n_obs = N, 
                          dims = dims,
                          true_x = x_true, 
                          s_x = s_x, 
                          rho = params$rho,
                          alpha = params$alpha,
                          sigma = params$sigma,
                          intc = params$intc,
                          covfn = covfn,
                          delta = delta)
y <- data[,1:dims]
x_obs <- data$x_obs[1:N]
x_true <- data$x_true[1:N]
f <- data[,(dims+1):(2*dims)]
# Check generated GP
dim_id <- 1
fplot <- gp_plot(dim_id = dim_id, output = y, gpfns = f, input = x_true, 
                 msd = params$alpha, esd = params$sigma, 
                 plot_type = 'gpfns', label = 'f')
fplot

# fit hsgp se model
min_m <- 33 #ceiling(1.75 * (5/4) * (max(x_true) - min(x_true))/(2 * ls_params[1])) # mean(rho) = 0.2, S = 0.5 (half range of x), c = 5/4
hsgp_se_fit <- test_hsgp_model(testhsgpmodel = testderivhsgpmodel,
                                  n_obs = N, 
                                  dims = dims,
                                  m_obs = min_m,
                                  x_min = x_min,
                                  x_max = x_max,
                                  c = 5/4,
                                  y = y,
                                  inputs = x_obs, 
                                  latent_sd = s_x, 
                                  latent_inputs = 1, 
                                  covfn = 0, 
                                  rho_prior = rho_prior_model,
                                  ls_param = ls_params_model,
                                  msd_param = marginal_sd_params_model,
                                  esd_param = error_sd_params_model,
                                  intc = intc_params_model,
                                  is_vary = 1, 
                                  is_corr = 1,
                                  iter = 2000, 
                                  warmup = 1000, 
                                  chains = 2, 
                                  cores = 2, 
                                  init = 0, 
                                  adapt_delta = adapt_delta_model)
# model output summary
print(hsgp_se_fit, pars = c('rho', 'alpha', 'sigma','x'))
traceplot(hsgp_se_fit, pars = 'x')#c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))

# fit hsgp se deriv model
min_m <- 11 #ceiling(1.75 * (5/4) * (max(x_true) - min(x_true))/(2 * ls_params[1])) # mean(rho) = 0.2, S = 0.5 (half range of x), c = 5/4
hsgp_se_deriv_fit <- test_hsgp_model(testhsgpmodel = testderivhsgpmodel,
                       n_obs = N, 
                       dims = dims,
                       m_obs = min_m,
                       x_min = x_min,
                       x_max = x_max,
                       c = 5/4,
                       y = y,
                       inputs = x_obs, 
                       latent_sd = s_x, 
                       latent_inputs = 1, 
                       covfn = 3, 
                       rho_prior = rho_prior_model,
                       ls_param = ls_params_model,
                       msd_param = marginal_sd_params_model,
                       esd_param = error_sd_params_model,
                       intc = intc_params_model,
                       is_vary = 1, 
                       is_corr = 1,
                       iter = 2000, 
                       warmup = 1000, 
                       chains = 2, 
                       cores = 2, 
                       init = 0, 
                       adapt_delta = adapt_delta_model)
# model output summary
print(hsgp_se_deriv_fit, pars = c('rho', 'alpha', 'sigma','x'))
traceplot(hsgp_se_deriv_fit, pars = 'x')#c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))

## Exact SE fit
delta <- 1e-6
gp_se_fit <- test_gp_model(testgpmodel = testderivgpmodel,
                           n_obs = N, 
                           dims = dims,
                           x_min = x_min,
                           x_max = x_max,
                           y = y,
                           inputs = x_obs, 
                           latent_sd = s_x, 
                           latent_inputs = 1, 
                           covfn = 0, 
                           delta = delta,
                           rho_prior = rho_prior_model,
                           ls_param = ls_params_model,
                           msd_param = marginal_sd_params_model,
                           esd_param = error_sd_params_model,
                           intc = intc_params_model,
                           is_vary = 1, 
                           is_corr = 1,
                           iter = 2000, 
                           warmup = 1000, 
                           chains = 2, 
                           cores = 2, 
                           init = 0, 
                           adapt_delta = adapt_delta_model)
# model output summary
print(gp_se_fit, pars = c('rho', 'alpha', 'sigma','x'))
traceplot(gp_se_fit, pars = 'x')#c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))

## Exact SE deriv fit
gp_se_deriv_fit <- test_gp_model(testgpmodel = testderivgpmodel,
                                 n_obs = N, 
                                 dims = dims,
                                 x_min = x_min,
                                 x_max = x_max,
                                 y = y,
                                 inputs = x_obs, 
                                 latent_sd = s_x, 
                                 latent_inputs = 1, 
                                 covfn = 3, 
                                 delta = delta,
                                 rho_prior = rho_prior_model,
                                 ls_param = ls_params_model,
                                 msd_param = marginal_sd_params_model,
                                 esd_param = error_sd_params_model,
                                 intc = intc_params_model,
                                 is_vary = 1, 
                                 is_corr = 1,
                                 iter = 2000, 
                                 warmup = 1000, 
                                 chains = 2, 
                                 cores = 2, 
                                 init = 0, 
                                 adapt_delta = adapt_delta_model)
# model output summary
print(gp_se_deriv_fit, pars = c('rho', 'alpha', 'sigma','x'))
traceplot(gp_se_deriv_fit, pars = 'x')#c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))

#### post summary HSGP
x_names <- sprintf('x[%s]', seq(1:N))
# Summary for hsgp se
out_hsgp_se <- as_draws_matrix(hsgp_se_fit)
subset_hsgp_se <- subset_draws(out_hsgp_se, variable = x_names)
sd_hsgp_se <- rep(NA, length(x_names))
abs_bias_hsgp_se <- rep(NA, length(x_names))
rmse_hsgp_se <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd_hsgp_se[i] <- sd(subset_hsgp_se[,i])
  abs_bias_hsgp_se[i] <- abs_bias_draws(subset_hsgp_se[,i], x_true[i])
  rmse_hsgp_se[i] <- rmse_draws(subset_hsgp_se[,i], x_true[i])
}
# Summary for hsgp se deriv
out_hsgp_se_deriv <- as_draws_matrix(hsgp_se_deriv_fit)
subset_hsgp_se_deriv <- subset_draws(out_hsgp_se_deriv, variable = x_names)
sd_hsgp_se_deriv <- rep(NA, length(x_names))
abs_bias_hsgp_se_deriv <- rep(NA, length(x_names))
rmse_hsgp_se_deriv <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd_hsgp_se_deriv[i] <- sd(subset_hsgp_se_deriv[,i])
  abs_bias_hsgp_se_deriv[i] <- abs_bias_draws(subset_hsgp_se_deriv[,i], x_true[i])
  rmse_hsgp_se_deriv[i] <- rmse_draws(subset_hsgp_se_deriv[,i], x_true[i])
}
# Summary GP se
out_gp_se <- as_draws_matrix(gp_se_fit)
subset_gp_se <- subset_draws(out_gp_se, variable = x_names)
sd_gp_se <- rep(NA, length(x_names))
abs_bias_gp_se <- rep(NA, length(x_names))
rmse_gp_se <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd_gp_se[i] <- sd(subset_gp_se[,i])
  abs_bias_gp_se[i] <- abs_bias_draws(subset_gp_se[,i], x_true[i])
  rmse_gp_se[i] <- rmse_draws(subset_gp_se[,i], x_true[i])
}
# Summary GP se deriv
out_gp_se_deriv <- as_draws_matrix(gp_se_deriv_fit)
subset_gp_se_deriv <- subset_draws(out_gp_se_deriv, variable = x_names)
sd_gp_se_deriv <- rep(NA, length(x_names))
abs_bias_gp_se_deriv <- rep(NA, length(x_names))
rmse_gp_se_deriv <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd_gp_se_deriv[i] <- sd(subset_gp_se_deriv[,i])
  abs_bias_gp_se_deriv[i] <- abs_bias_draws(subset_gp_se_deriv[,i], x_true[i])
  rmse_gp_se_deriv[i] <- rmse_draws(subset_gp_se_deriv[,i], x_true[i])
}

# mean(rmse_gp_se)
# mean(abs_bias_gp_se)
# mean(sd_gp_se)
# mean(rmse_gp_se_deriv)
# mean(abs_bias_gp_se_deriv)
# mean(sd_gp_se_deriv)
# 
# mean(rmse_hsgp_se)
# mean(abs_bias_hsgp_se)
# mean(sd_hsgp_se)
# mean(rmse_hsgp_se_deriv)
# mean(abs_bias_hsgp_se_deriv)
# mean(sd_hsgp_se_deriv)

# compare with naive estimate on the basis of the x_obs prior
n_obs <- 100000
rmse_x_naive <- vector(length = n_obs)
bias_x_naive <- vector(length = n_obs)
s_x <- 0.3
for (i in 1:n_obs) {
  # as if we only used x_obs to infer x
  x_obs <- rnorm(1, 0, s_x)
  draws_x_prior <- rnorm(2000, x_obs, s_x)
  rmse_x_naive[i] <- rmse_draws(draws_x_prior, 0)
  bias_x_naive[i] <- abs_bias_draws(draws_x_prior, 0)
}

summary_results_df <- data.frame(
  Model = c('Exact SE', 'Exact SE Deriv', 'HSGP SE', 'HSGP SE Deriv', 'Prior'),
  SD = c(mean(sd_gp_se), mean(sd_gp_se_deriv), mean(sd_hsgp_se), mean(sd_hsgp_se_deriv), s_x),
  AbsBias = c(mean(abs_bias_gp_se), mean(abs_bias_gp_se_deriv), mean(abs_bias_hsgp_se), mean(abs_bias_hsgp_se_deriv), mean(bias_x_naive)),
  RMSE = c(mean(rmse_gp_se), mean(rmse_gp_se_deriv), mean(rmse_hsgp_se), mean(rmse_hsgp_se_deriv), mean(rmse_x_naive))
)



# Exact Derivative GP
gp_fit <- deriv_gp_model(stanmodel = derivgpmodel,
                         n_obs = N, 
                         dims = dims, 
                         outputs = rbind(y1, y2), 
                         inputs = x_obs, 
                         latent_sd = s_x, 
                         latent_inputs = 1, 
                         covfn = covfn_model, 
                         rho_prior = rho_prior_model,
                         ls_param = ls_params_model, 
                         msd_param = marginal_sd_params_model, 
                         esd_param = error_sd_params_model, 
                         msd_param_grad = marginal_sd_params_grad_model,
                         esd_param_grad = error_sd_params_grad_model,
                         intc_y = intc_params_model_y,
                         intc_yderiv = intc_params_model_yderiv,
                         is_deriv = 1,
                         is_scale = 1,
                         is_vary = 1,
                         is_corr = 1,
                         iter = 2000, 
                         warmup = 1000, 
                         chains = 2, 
                         cores = 2, 
                         init = 0, 
                         adapt_delta = adapt_delta_model)
print(gp_fit, pars = c('rho', 'alpha_obs', 'sigma_obs', 'alpha_grad', 'sigma_grad','x'))
traceplot(gp_fit, pars = 'x')#c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))

