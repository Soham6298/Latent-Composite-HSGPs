#library
library(rstan)
library(posterior)
library(ggplot2)
library(bayesplot)

# Source fns 
source('indcompgpfns.R')
comphsgpmodel <- stan_model('indpcomphsgp_maternclass.stan')
compgpmodel <- stan_model('indpcompexactGP_maternclass.stan')
derivgpmodel <- stan_model('DerivGPmodels.stan')
set.seed(42321)
## Specify data conditions
N <- 50  # higher N will take very long. Better to stick to <=30
dims <-10
# Set true input x
#x_true <- seq(0.5, 9.5, length.out = N)
x_min <- 0
x_max <- 1
s_x <- 0.3 # latent sd # needed in case of latent GPs
# Generate parameters
marginal_sd_params <- c(3, 0.25)
error_sd_params <- c(1, 0.25)
ls_params <- c(0.3, 0.25)
intc_params_y <- c(0, 1)
intc_params_yderiv <- c(0, 1)
delta <- 1e-12
lambda <- 3        # scale between y and y'
covfn <- 'se'

# Model specifications
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
lambda <- 3  # derivative scale factor
marginal_sd_params_grad_model <- marginal_sd_params # alpha_grad ~ N( , )
error_sd_params_grad_model <- error_sd_params    # sigma_grad ~ N( , )
marginal_sd_params_model <- marginal_sd_params_grad_model * lambda  # alpha_obs ~ N( , )
error_sd_params_model <- error_sd_params_grad_model * lambda  # sigma_obs ~ N( , )
ls_params_model <- ls_params # rho ~ N( , ) or InvGamma( , )
#ls_params_model_deriv <- c(0.7, 0.05)
intc_params_model_y <- intc_params_y
intc_params_model_yderiv <- intc_params_yderiv
covfn_model <- 0  # 0 = SE; 1 = Matern3/2; 2 = Matern 5/2
# Set covariance function based scaler for HSGPs
if(covfn_model==1){
  const <- 3.42
}else if(covfn_model==2){
  const <- 2.65
}else {
  const <- 1.75
}
min_m <- ceiling(const * (5/4) * (x_max - x_min)/(2*ls_params_model[1])) 
adapt_delta_model <- 0.95
x_true <- runif(N, x_min, x_max)
params <- true_params_deriv_vary(n_obs = N, 
                                 dims = dims, 
                                 deriv_scale = lambda,
                                 msd_pars = marginal_sd_params,
                                 esd_pars = error_sd_params,
                                 ls_pars = ls_params,
                                 intc_pars = intc_params_y,
                                 intc_deriv_pars = intc_params_yderiv)
data <- deriv_gp_sim_data(n_obs = N, 
                          dims = dims,
                          true_x = x_true, 
                          s_x = s_x, 
                          rho = params$rho,
                          alpha_obs = params$alpha_obs, 
                          alpha_grad = params$alpha_grad,
                          sigma_obs = params$sigma_obs,
                          sigma_grad = params$sigma_grad,
                          intc = params$intc,
                          intc_deriv = params$intc_deriv,
                          covfn = covfn,
                          delta = delta)
y1 <- data[1:N,1:dims]
y2 <- data[(N+1):(2*N),1:dims]
x_obs <- data$x_obs[1:N]
x_true <- data$x_true[1:N]
f <- data[1:N,(dims+1):(2*dims)]
g <- data[(N+1):(2*N),(dims+1):(2*dims)]
# Check generated GP
dim_id <- 1
fplot <- gp_plot(dim_id = dim_id, output = y1, gpfns = f, input = x_true, 
          msd = params$alpha_obs, esd = params$sigma_obs, 
          plot_type = 'gpfns', label = 'f')
gplot <- gp_plot(dim_id = dim_id, output = y2, gpfns = g, input = x_true, 
       msd = params$alpha_grad, esd = params$sigma_grad, 
       plot_type = 'gpfns', label = 'g')
bayesplot_grid(fplot, gplot)
# fit hsgp model
min_m <- ceiling(1.75 * (5/4) * (max(x_true) - min(x_true))/(2 * ls_params[1])) # mean(rho) = 0.2, S = 0.5 (half range of x), c = 5/4
#min_m_g <- ceiling(1.75 * (5/4) * (max(x_true) - min(x_true))/(2 * ls_params[1])) # mean(rho) = 0.2, S = 0.5 (half range of x), c = 5/4
hsgp_fit <- hsgp_model(hsgpmodel = comphsgpmodel,
                       n_obs = N, 
                       dims = dims,
                       m_obs_f = min_m,
                       m_obs_g = min_m,
                       x_min = x_min,
                       x_max = x_max,
                       c = 5/4,
                       y_f = y1, 
                       y_g = y2, 
                       inputs = x_obs, 
                       latent_sd = s_x, 
                       latent_inputs = 1, 
                       covfn = covfn_model, 
                       rho_prior = rho_prior_model,
                       ls_param_f = ls_params_model,
                       ls_param_g = ls_params_model, 
                       msd_param_f = marginal_sd_params_model,
                       msd_param_g = marginal_sd_params_grad_model,
                       esd_param_f = error_sd_params_model,
                       esd_param_g = error_sd_params_grad_model,
                       intc_yf = intc_params_model_y,
                       intc_yg = intc_params_model_yderiv,
                       is_vary = 1, 
                       is_corr = 1,
                       is_deriv_f = 0,
                       is_deriv_g = 0,
                       iter = 2000, 
                       warmup = 1000, 
                       chains = 2, 
                       cores = 2, 
                       init = 0, 
                       adapt_delta = adapt_delta_model)
# model output summary
print(hsgp_fit, pars = c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))
traceplot(hsgp_fit, pars = 'x')#c('rho_f', 'alpha_f', 'sigma_f', 'rho_g', 'alpha_g', 'sigma_g','x'))

# Exact GP
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
#### post summary HSGP
x_names <- sprintf('x[%s]', seq(1:N))
out_hsgp <- as_draws_matrix(hsgp_fit)
subset_hsgp <- subset_draws(out_hsgp, variable = x_names)
summary_hsgp <- summarise_draws(subset_hsgp)
plot(hist(summary_hsgp$rhat), main = 'Rhats with s = 0.3')
sd_hsgp <- rep(NA, length(x_names))
abs_bias_hsgp <- rep(NA, length(x_names))
rmse_hsgp <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd_hsgp[i] <- sd(subset_hsgp[,i])
  abs_bias_hsgp[i] <- abs_bias_draws(subset_hsgp[,i], x_true[i])
  rmse_hsgp[i] <- rmse_draws(subset_hsgp[,i], x_true[i])
}
#### post summary GP
out_gp <- as_draws_matrix(gp_fit)
subset_gp <- subset_draws(out_gp, variable = x_names)
summary_gp <- summarise_draws(subset_gp)
sd_gp <- rep(NA, length(x_names))
abs_bias_gp <- rep(NA, length(x_names))
rmse_gp <- rep(NA, length(x_names))
for(i in 1:length(x_names)) {
  sd_gp[i] <- sd(subset_gp[,i])
  abs_bias_gp[i] <- abs_bias_draws(subset_gp[,i], x_true[i])
  rmse_gp[i] <- rmse_draws(subset_gp[,i], x_true[i])
}
mean(rmse_gp)
mean(abs_bias_gp)
mean(sd_gp)
mean(rmse_hsgp)
mean(abs_bias_hsgp)
mean(sd_hsgp)

# compare with naive estimate on the basis of the x_obs prior
n_obs <- 100000
rmse_x_naive <- vector(length = n_obs)
mae_x_naive <- vector(length = n_obs)
s_x <- 0.3
for (i in 1:n_obs) {
  # as if we only used x_obs to infer x
  x_obs <- rnorm(1, 0, s_x)
  draws_x_prior <- rnorm(2000, x_obs, s_x)
  rmse_x_naive[i] <- rmse_draws(draws_x_prior, 0)
  mae_x_naive[i] <- mae_draws(draws_x_prior, 0)
}
naive_rmse <- mean(rmse_x_naive)
naive_mae <- mean(mae_x_naive)
naive_rmse



##########
# N <- 50  # higher N will take very long. Better to stick to <=30
# dims <- 10
# # Set true input x
# x_min <- 0
# x_max <- 10
# #x_true <- seq(0.5, by = 0.5, length.out = N) #length.out = N)
# s <- 0.3 # latent sd # needed in case of latent GPs
# # Generate parameters
# marginal_sd_params_f <- c(3, 0.25)
# error_sd_params_f <- c(1, 0.25)
# ls_params_f <- c(1, 0.05)
# marginal_sd_params_g <- c(2, 0.25)
# error_sd_params_g <- c(0.7, 0.25)
# ls_params_g <- c(0.7, 0.05)
# intc_params_y1 <- c(0, 1)
# intc_params_y2 <- c(0, 10)
# delta <- 1e-12
# covfn <- 'se' # only needed for gp data
# 
# # Model specifications
# rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# # Specify priors for the hyperparameters (check if they differ from the simulation setting)
# marginal_sd_params_model_f <- marginal_sd_params_f
# error_sd_params_model_f <- error_sd_params_f 
# ls_params_model_f <- ls_params_f 
# marginal_sd_params_model_g <- marginal_sd_params_g
# error_sd_params_model_g <- error_sd_params_g
# ls_params_model_g <- ls_params_g
# intc_params_model_y1 <- intc_params_y1
# intc_params_model_y2 <- intc_params_y2
# latent_model <- 1 # 1 for latent inputs, 0 for manifest
# covfn_model <- 0  # 0 - SE; 1 = Matern3/2; 2 = Matern 5/2
# adapt_delta_model <- 0.95
# 
# input_dat <- sim_input(n_obs = N, x_min = x_min, x_max = x_max, s_x = s)
# x_true <- input_dat$x_true
# x_obs <- input_dat$x_obs
# params1 <- true_params_vary(n_obs = N, 
#                                dims = dims, 
#                                msd_pars = marginal_sd_params_f, 
#                                esd_pars = error_sd_params_f, 
#                                ls_pars = ls_params_f, 
#                                intc_pars = intc_params_y1) 
# params2 <- true_params_vary(n_obs = N, 
#                              dims = dims, 
#                              msd_pars = marginal_sd_params_g, 
#                              esd_pars = error_sd_params_g, 
#                              ls_pars = ls_params_g, 
#                              intc_pars = intc_params_y2) 
# 
# # Simulate data
# dat1 <- gp_sim_data(n_obs = N, 
#                      dims = dims,
#                      true_x = x_true,
#                      covfn = covfn,
#                      delta = delta,
#                      alpha = params1$alpha, 
#                      sigma = params1$sigma,
#                      rho = params1$rho, 
#                    intc = params1$intc)
# dat2 <- gp_sim_data(n_obs = N, 
#                     dims = dims, 
#                     true_x = x_true,
#                     covfn = covfn,
#                     delta = delta,
#                     alpha = params2$alpha, 
#                     sigma = params2$sigma,
#                     rho = params2$rho, 
#                     intc = params2$intc)
# y1 <- dat1[1:N,1:dims]
# f <- dat1[1:N,(dims+1):(2*dims)]
# y2 <- dat2[1:N,1:dims]
# g <- dat2[1:N,(dims+1):(2*dims)]
