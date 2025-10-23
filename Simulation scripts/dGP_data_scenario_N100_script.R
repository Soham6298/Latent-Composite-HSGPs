# Here we simulate the GP data scenario for DGP-LVM and other GP model specifications
# We run multiple trials to verify the effects of the modifications in covariance functions
# Ground truth data is simulated from derivative GP with all our proposed modifications

#libraries
library(rstan)
library(parallel)
library(doParallel)
library(foreach)
library(posterior)
library(data.table)
# Source fns file
source('indcompgpfns.R')
idhsgpmodel <- stan_model('pdHSGP.stan')
ihsgpmodel <- stan_model('pcHSGP.stan')
singlehsgpmodel <- stan_model('sHSGP.stan')
# tempdir set (set temp dir so that it doesn't overload disk: primarily for cluster computing)
Sys.setenv(TMPDIR = "/mnt/volume")
unlink(tempdir(), recursive = TRUE)
tempdir(check = TRUE)
#functions

## Specify data conditions
N <- 100  # higher N will take very long. Better to stick to <=30
dims <- c(5, 10, 20)
# Set true input x
x_min <- 0
x_max <- 10
s_x <- 0.3 # latent sd # needed in case of latent GPs
# Generate parameters
marginal_sd_params <- c(3, 0.25)
error_sd_params <- c(1, 0.25)
ls_params <- c(1, 0.05)
intc_params_y <- c(0, 5)
intc_params_yderiv <- c(0, 5)
delta <- 1e-12
covfn <- 'se'

# Model specifications
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
lambda <- 10  # derivative scale factor
marginal_sd_params_grad_model <- marginal_sd_params # alpha_grad ~ N( , )
error_sd_params_grad_model <- error_sd_params    # sigma_grad ~ N( , )
marginal_sd_params_model <- marginal_sd_params_grad_model * lambda  # alpha_obs ~ N( , )
error_sd_params_model <- error_sd_params_grad_model * lambda  # sigma_obs ~ N( , )
ls_params_model <- ls_params # rho ~ N( , ) or InvGamma( , )
intc_params_model_y <- intc_params_y
intc_params_model_yderiv <- intc_params_yderiv
covfn_model <- 0
# Set min basis functions for hsgp
min_m <- 30 
adapt_delta_model <- 0.95

# Generate sim data
n_sim <- 50
params <- list()
simdata <- list()
for (i in 1:n_sim) {
  set.seed(i)
  params[[i]] <- list()
  simdata[[i]] <- list() 
  for (j in 1:length(dims)) {
    x_true <- runif(N, x_min, x_max)
    params[[i]][[j]] <- true_params_deriv_vary(n_obs = N, 
                                               dims = dims[j], 
                                               deriv_scale = lambda,
                                               msd_pars = marginal_sd_params,
                                               esd_pars = error_sd_params,
                                               ls_pars = ls_params,
                                               intc_pars = intc_params_y,
                                               intc_deriv_pars = intc_params_yderiv)
    simdata[[i]][[j]] <- deriv_gp_sim_data(n_obs = N, 
                                           dims = dims[j],
                                           true_x = x_true, 
                                           s_x = s_x, 
                                           rho = params[[i]][[j]]$rho,
                                           alpha_obs = params[[i]][[j]]$alpha_obs, 
                                           alpha_grad = params[[i]][[j]]$alpha_grad,
                                           sigma_obs = params[[i]][[j]]$sigma_obs,
                                           sigma_grad = params[[i]][[j]]$sigma_grad,
                                           intc = params[[i]][[j]]$intc,
                                           intc_deriv = params[[i]][[j]]$intc_deriv,
                                           covfn = covfn,
                                           delta = delta)
  }
}
saveRDS(params, 'deriv_Params_hsgp_n100.rds')
saveRDS(simdata, 'deriv_simdata_hsgp_n100.rds')

# Sort variable names for output
x_names <- sprintf('x[%s]', seq(1:N))
rho_names <- list()
alpha_names <- list()
sigma_names <- list()
rho_f_names <- list()
alpha_f_names <- list()
sigma_f_names <- list()
rho_g_names <- list()
alpha_g_names <- list()
sigma_g_names <- list()
for( j in 1:length(dims)){
  rho_f_names[[j]] <- sprintf('rho_f[%s]', seq(1:dims[j]))
  alpha_f_names[[j]] <- sprintf('alpha_f[%s]', seq(1:dims[j]))
  sigma_f_names[[j]] <- sprintf('sigma_f[%s]', seq(1:dims[j]))
  rho_g_names[[j]] <- sprintf('rho_g[%s]', seq(1:dims[j]))
  alpha_g_names[[j]] <- sprintf('alpha_g[%s]', seq(1:dims[j]))
  sigma_g_names[[j]] <- sprintf('sigma_g[%s]', seq(1:dims[j]))
  rho_names[[j]] <- sprintf('rho[%s]', seq(1:dims[j]))
  alpha_names[[j]] <- sprintf('alpha[%s]', seq(1:dims[j]))
  sigma_names[[j]] <- sprintf('sigma[%s]', seq(1:dims[j]))
}

#start cluster
# parallel for deriv model fits
out_obsonly <- list()
out_obsonly1 <- list()
cores <- 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_obshsgp = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  set.seed(i)
  out_obsonly1[[i]] <- list()
  for (j in 1:length(dims)) {
    y <- simdata[[i]][[j]][1:N,1:dims[j]]
    x_obs <- simdata[[i]][[j]]$x_obs[1:N]
    x_true <- simdata[[i]][[j]]$x_true[1:N]
    hsgp_fit <- singlehsgp_model(hsgpmodel = singlehsgpmodel,
                                 n_obs = N, 
                                 m_obs = min_m,
                                 dims = dims[j], 
                                 y = y, 
                                 x_min = x_min,
                                 x_max = x_max,
                                 c = 5/4,
                                 inputs = x_obs, 
                                 latent_sd = s_x, 
                                 latent_inputs = 1, 
                                 covfn = 0, 
                                 rho_prior = rho_prior_model,
                                 ls_param = ls_params_model, 
                                 msd_param = marginal_sd_params_model, 
                                 esd_param = error_sd_params_model, 
                                 intc_y = intc_params_model_y,
                                 is_vary = 1,
                                 is_corr = 1,
                                 iter = 2000, 
                                 warmup = 1000, 
                                 chains = 1, 
                                 cores = 1, 
                                 init = 0, 
                                 adapt_delta = adapt_delta_model)
    out_hsgp_x <- compare_summary(model = hsgp_fit, 
                                variable = x_names, 
                                true_variable = simdata[[i]][[j]]$x_true[1:N], 
                                n_obs = N, 
                                m_approx = 'obs_hsgp',
                                dims = dims[j], 
                                variable_class = 'x',
                                model_name = 'obs_hsgp',
                                sim_id = i)
    out_hsgp_rho <- compare_summary(model = hsgp_fit, 
                                  variable = rho_names[[j]], 
                                  true_variable = params[[i]][[j]]$rho, 
                                  n_obs = N, 
                                  m_approx = 'obs_hsgp',
                                  dims = dims[j],
                                  variable_class = 'rho', 
                                  model_name = 'obs_hsgp',
                                  sim_id = i)
    out_hsgp_alpha <- compare_summary(model = hsgp_fit, 
                                        variable = alpha_names[[j]], 
                                        true_variable = params[[i]][[j]]$alpha_obs, 
                                        n_obs = N,
                                        m_approx = 'obs_hsgp',
                                        dims = dims[j], 
                                        variable_class = 'alpha',
                                        model_name = 'obs_hsgp',
                                        sim_id = i)
    out_hsgp_sigma <- compare_summary(model = hsgp_fit, 
                                        variable = sigma_names[[j]], 
                                        true_variable = params[[i]][[j]]$sigma_obs,
                                        n_obs = N, 
                                        m_approx = 'obs_hsgp',
                                        dims = dims[j], 
                                        variable_class = 'sigma',
                                        model_name = 'obs_hsgp',
                                        sim_id = i)
    out_obsonly1[[i]][[j]] <- rbind(out_hsgp_x, out_hsgp_rho, out_hsgp_alpha, out_hsgp_sigma)
  }
  out_obsonly[[i]] <- rbindlist(out_obsonly1[[i]])
}
compare_obshsgp <- rbindlist(model_fit_obshsgp)

# parallel for HSGP
out_derivonly <- list()
out_derivonly1 <- list()
cores = 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_derivhsgp = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  set.seed(i)
  out_derivonly1[[i]] <- list()
  for (j in 1:length(dims)) {
    y <- simdata[[i]][[j]][(N+1):(2*N),1:dims[j]]
    x_obs <- simdata[[i]][[j]]$x_obs[1:N]
    x_true <- simdata[[i]][[j]]$x_true[1:N]
    hsgp_fit <- singlehsgp_model(hsgpmodel = singlehsgpmodel,
                           n_obs = N, 
                           m_obs = min_m,
                           dims = dims[j], 
                           y = y, 
                           x_min = x_min,
                           x_max = x_max,
                           c = 5/4,
                           inputs = x_obs, 
                           latent_sd = s_x, 
                           latent_inputs = 1, 
                           covfn = 3, 
                           rho_prior = rho_prior_model,
                           ls_param = ls_params_model, 
                           msd_param = marginal_sd_params_grad_model, 
                           esd_param = error_sd_params_grad_model, 
                           intc_y = intc_params_model_yderiv,
                           is_vary = 1,
                           is_corr = 1,
                           iter = 2000, 
                           warmup = 1000, 
                           chains = 1, 
                           cores = 1, 
                           init = 0, 
                           adapt_delta = adapt_delta_model)
    out_hsgp_x <- compare_summary(model = hsgp_fit, 
                                  variable = x_names, 
                                  true_variable = simdata[[i]][[j]]$x_true[1:N], 
                                  n_obs = N,
                                  m_approx = 'deriv_hsgp', 
                                  dims = dims[j],
                                  variable_class = 'x', 
                                  model_name = 'deriv_hsgp', 
                                  sim_id = i)
    out_hsgp_rho <- compare_summary(model = hsgp_fit,
                                      variable = rho_names[[j]], 
                                      true_variable = params[[i]][[j]]$rho,
                                      n_obs = N, 
                                      m_approx = 'deriv_hsgp', 
                                      dims = dims[j], 
                                      variable_class = 'rho', 
                                      model_name = 'deriv_hsgp',
                                      sim_id = i)
    out_hsgp_alpha <- compare_summary(model = hsgp_fit, 
                                        variable = alpha_names[[j]], 
                                        true_variable = params[[i]][[j]]$alpha_grad, 
                                        n_obs = N, 
                                        m_approx = 'deriv_hsgp', 
                                        dims = dims[j], 
                                        variable_class = 'alpha', 
                                        model_name = 'deriv_hsgp', 
                                        sim_id = i)
    out_hsgp_sigma <- compare_summary(model = hsgp_fit,
                                        variable = sigma_names[[j]], 
                                        true_variable = params[[i]][[j]]$sigma_grad, 
                                        n_obs = N, 
                                        m_approx = 'deriv_hsgp', 
                                        dims = dims[j], 
                                        variable_class = 'sigma',
                                        model_name = 'deriv_hsgp', 
                                        sim_id = i)
    out_derivonly1[[i]][[j]] <- rbind(out_hsgp_x, out_hsgp_rho, out_hsgp_alpha, out_hsgp_sigma)
  }
  out_derivonly[[i]] <- rbindlist(out_derivonly1[[i]])
}
compare_derivhsgp <- rbindlist(model_fit_derivhsgp)

out_ihsgp <- list()
out_ihsgp1 <- list()
cores = 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_ihsgp = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  set.seed(i)
  out_ihsgp1[[i]] <- list()
  for (j in 1:length(dims)) {
    y1 <- simdata[[i]][[j]][1:N,1:dims[j]]
    y2 <- simdata[[i]][[j]][(N+1):(2*N),1:dims[j]]
    x_obs <- simdata[[i]][[j]]$x_obs[1:N]
    x_true <- simdata[[i]][[j]]$x_true[1:N]
    hsgp_fit <- ihsgp_model(hsgpmodel = ihsgpmodel,
                            n_obs = N, 
                            dims = dims[j],
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
                            iter = 2000, 
                            warmup = 1000, 
                            chains = 1, 
                            cores = 1, 
                            init = 0, 
                            adapt_delta = adapt_delta_model)
    out_ihsgp_x <- compare_summary(model = hsgp_fit, 
                                   variable = x_names, 
                                   true_variable = simdata[[i]][[j]]$x_true[1:N], 
                                   n_obs = N,
                                   m_approx = 'ihsgp', 
                                   dims = dims[j],
                                   variable_class = 'x', 
                                   model_name = 'ihsgp', 
                                   sim_id = i)
    out_ihsgp_rho_f <- compare_summary(model = hsgp_fit,
                                       variable = rho_f_names[[j]], 
                                       true_variable = params[[i]][[j]]$rho,
                                       n_obs = N, 
                                       m_approx = 'ihsgp', 
                                       dims = dims[j], 
                                       variable_class = 'rho_f', 
                                       model_name = 'ihsgp',
                                       sim_id = i)
    out_ihsgp_alpha_f <- compare_summary(model = hsgp_fit, 
                                         variable = alpha_f_names[[j]], 
                                         true_variable = params[[i]][[j]]$alpha_obs, 
                                         n_obs = N, 
                                         m_approx = 'ihsgp', 
                                         dims = dims[j], 
                                         variable_class = 'alpha_f', 
                                         model_name = 'ihsgp', 
                                         sim_id = i)
    out_ihsgp_sigma_f <- compare_summary(model = hsgp_fit,
                                         variable = sigma_f_names[[j]], 
                                         true_variable = params[[i]][[j]]$sigma_obs, 
                                         n_obs = N, 
                                         m_approx = 'ihsgp', 
                                         dims = dims[j], 
                                         variable_class = 'sigma_f',
                                         model_name = 'ihsgp', 
                                         sim_id = i)
    out_ihsgp_rho_g <- compare_summary(model = hsgp_fit,
                                       variable = rho_g_names[[j]], 
                                       true_variable = params[[i]][[j]]$rho,
                                       n_obs = N, 
                                       m_approx = 'ihsgp', 
                                       dims = dims[j], 
                                       variable_class = 'rho_g', 
                                       model_name = 'ihsgp',
                                       sim_id = i)
    out_ihsgp_alpha_g <- compare_summary(model = hsgp_fit, 
                                         variable = alpha_g_names[[j]], 
                                         true_variable = params[[i]][[j]]$alpha_grad, 
                                         n_obs = N, 
                                         m_approx = 'ihsgp', 
                                         dims = dims[j], 
                                         variable_class = 'alpha_g', 
                                         model_name = 'ihsgp', 
                                         sim_id = i)
    out_ihsgp_sigma_g <- compare_summary(model = hsgp_fit,
                                         variable = sigma_g_names[[j]], 
                                         true_variable = params[[i]][[j]]$sigma_grad, 
                                         n_obs = N, 
                                         m_approx = 'ihsgp', 
                                         dims = dims[j], 
                                         variable_class = 'sigma_g',
                                         model_name = 'ihsgp', 
                                         sim_id = i)
    out_ihsgp1[[i]][[j]] <- rbind(out_ihsgp_x, out_ihsgp_rho_f, out_ihsgp_alpha_f, out_ihsgp_sigma_f,
                                  out_ihsgp_rho_g, out_ihsgp_alpha_g, out_ihsgp_sigma_g)
  }
  out_ihsgp[[i]] <- rbindlist(out_ihsgp1[[i]])
}
compare_ihsgp <- rbindlist(model_fit_ihsgp)

# parallel for HSGP deriv
out_idhsgp <- list()
out_idhsgp1 <- list()
cores = 50
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_idhsgp = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  set.seed(i)
  out_idhsgp1[[i]] <- list()
  for (j in 1:length(dims)) {
    y1 <- simdata[[i]][[j]][1:N,1:dims[j]]
    y2 <- simdata[[i]][[j]][(N+1):(2*N),1:dims[j]]
    x_obs <- simdata[[i]][[j]]$x_obs[1:N]
    x_true <- simdata[[i]][[j]]$x_true[1:N]
    hsgp_fit <- idhsgp_model(hsgpmodel = idhsgpmodel,
                             n_obs = N, 
                             dims = dims[j],
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
                             rho_prior = rho_prior_model,
                             ls_param = ls_params_model,
                             msd_param_f = marginal_sd_params_model,
                             msd_param_g = marginal_sd_params_grad_model,
                             esd_param_f = error_sd_params_model,
                             esd_param_g = error_sd_params_grad_model,
                             intc_yf = intc_params_model_y,
                             intc_yg = intc_params_model_yderiv,
                             is_vary = 1, 
                             is_corr = 1,
                             iter = 2000, 
                             warmup = 1000, 
                             chains = 1, 
                             cores = 1, 
                             init = 0, 
                             adapt_delta = adapt_delta_model)
    out_idhsgp_x <- compare_summary(model = hsgp_fit, 
                                    variable = x_names, 
                                    true_variable = simdata[[i]][[j]]$x_true[1:N], 
                                    n_obs = N,
                                    m_approx = 'idhsgp', 
                                    dims = dims[j],
                                    variable_class = 'x', 
                                    model_name = 'idhsgp', 
                                    sim_id = i)
    out_idhsgp_rho <- compare_summary(model = hsgp_fit,
                                      variable = rho_names[[j]], 
                                      true_variable = params[[i]][[j]]$rho,
                                      n_obs = N, 
                                      m_approx = 'idhsgp', 
                                      dims = dims[j], 
                                      variable_class = 'rho', 
                                      model_name = 'idhsgp',
                                      sim_id = i)
    out_idhsgp_alpha_f <- compare_summary(model = hsgp_fit, 
                                          variable = alpha_f_names[[j]], 
                                          true_variable = params[[i]][[j]]$alpha_obs, 
                                          n_obs = N, 
                                          m_approx = 'idhsgp', 
                                          dims = dims[j], 
                                          variable_class = 'alpha_f', 
                                          model_name = 'idhsgp', 
                                          sim_id = i)
    out_idhsgp_sigma_f <- compare_summary(model = hsgp_fit,
                                          variable = sigma_f_names[[j]], 
                                          true_variable = params[[i]][[j]]$sigma_obs, 
                                          n_obs = N, 
                                          m_approx = 'idhsgp', 
                                          dims = dims[j], 
                                          variable_class = 'sigma_f',
                                          model_name = 'idhsgp', 
                                          sim_id = i)
    out_idhsgp_alpha_g <- compare_summary(model = hsgp_fit, 
                                          variable = alpha_g_names[[j]], 
                                          true_variable = params[[i]][[j]]$alpha_grad, 
                                          n_obs = N, 
                                          m_approx = 'idhsgp', 
                                          dims = dims[j], 
                                          variable_class = 'alpha_g', 
                                          model_name = 'idhsgp', 
                                          sim_id = i)
    out_idhsgp_sigma_g <- compare_summary(model = hsgp_fit,
                                          variable = sigma_g_names[[j]], 
                                          true_variable = params[[i]][[j]]$sigma_grad, 
                                          n_obs = N, 
                                          m_approx = 'idhsgp', 
                                          dims = dims[j], 
                                          variable_class = 'sigma_g',
                                          model_name = 'idhsgp', 
                                          sim_id = i)
    out_idhsgp1[[i]][[j]] <- rbind(out_idhsgp_x, out_idhsgp_rho, out_idhsgp_alpha_f, out_idhsgp_sigma_f, out_idhsgp_alpha_g, out_idhsgp_sigma_g)
  }
  out_idhsgp[[i]] <- rbindlist(out_idhsgp1[[i]])
}
compare_idhsgp <- rbindlist(model_fit_idhsgp)

compare_table <- rbind(compare_obshsgp, compare_derivhsgp, compare_ihsgp, compare_idhsgp)
## comparison table edit
compare_table$data_id = paste0(compare_table$sim_id, '_', 
                               compare_table$d, '_',
                               compare_table$m)
saveRDS(compare_table, 'dGP_data_scenario_n100.rds')
