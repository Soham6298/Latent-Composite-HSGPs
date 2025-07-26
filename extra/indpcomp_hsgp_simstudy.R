## Full simulation study for HSGPs and exact GPs. We generate data and fit models on cluster computers for 50 simulation trials.

#libraries
library(rstan)
library(posterior)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
# Import functions
source('indcompgpfns.R')
# Compile stan models
hsgpmodel <- stan_model('indpcomphsgp_maternclass.stan')
gpmodel <- stan_model('indpcompexactGP_maternclass.stan')
# Set temporary folder (only to prevent cluster memory overflow)
#Sys.setenv(TMPDIR = "/mnt/volume")
#unlink(tempdir(), recursive = TRUE)
#tempdir(check = TRUE)

# Specify data generating conditions
N <- 20
dims <- c(2, 3)
x_min <- 0
x_max <- 10
s_x <- 0.3
marginal_sd_params_f <- c(3, 0.25)
error_sd_params_f <- c(1, 0.25)
ls_params_f <- c(1, 0.05)
marginal_sd_params_g <- c(2, 0.25)
error_sd_params_g <- c(1, 0.25)
ls_params_g <- c(1.5, 0.05)
intc_params_y1 <- c(0, 1)
intc_params_y2 <- c(0, 5)
delta <- 1e-12
covfn <- 'se'

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
# Set covariance function based scaler for HSGPs
if(covfn_model==1){
  const <- 3.42
}else if(covfn_model==2){
  const <- 2.65
}else {
  const <- 1.75
}
# Compute minimum number of basis functions
min_m_f <- ceiling(const * (5/4) * (x_max - x_min)/(2*ls_params_model_f[1])) # mean(rho) = 0.3, S = 0.49 (half range of x), c = 5/4
min_m_g <- ceiling(const * (5/4) * (x_max - x_min)/(2*ls_params_model_g[1]))
# Set acceptance probability for MCMC 
adapt_delta_model <- 0.95

# Generate simulated data
n_sim <- 2
params1 <- list()
simdata1 <- list()
params2 <- list()
simdata2 <- list()
input_data <- list()
for (i in 1:n_sim) {
  set.seed(i)
  params1[[i]] <- list()
  params2[[i]] <- list()
  simdata1[[i]] <- list() 
  simdata2[[i]] <- list()
  input_data[[i]] <- list()
  for (j in 1:length(dims)) {
    input_data[[i]][[j]] <- sim_input(n_obs = N, x_min = x_min, x_max = x_max, s_x = s_x)
    params1[[i]][[j]] <- true_params_vary(n_obs = N, 
                                          dims = dims[j], 
                                          msd_pars = marginal_sd_params_f, 
                                          esd_pars = error_sd_params_f, 
                                          ls_pars = ls_params_f, 
                                          intc_pars = intc_params_y1) 
    params2[[i]][[j]] <- true_params_vary(n_obs = N, 
                                               dims = dims[j], 
                                               msd_pars = marginal_sd_params_g, 
                                               esd_pars = error_sd_params_g, 
                                               ls_pars = ls_params_g, 
                                               intc_pars = intc_params_y2)
    simdata1[[i]][[j]] <- gp_sim_data(n_obs = N, 
                                      dims = dims[j], 
                                      covfn = covfn,
                                      true_x = input_data[[i]][[j]]$x_true, 
                                      rho = params1[[i]][[j]]$rho, 
                                      alpha = params1[[i]][[j]]$alpha, 
                                      sigma = params1[[i]][[j]]$sigma,
                                      intc = params1[[i]][[j]]$intc)
    simdata2[[i]][[j]] <- gp_sim_data(n_obs = N, 
                                      dims = dims[j], 
                                      covfn = covfn,
                                      true_x = input_data[[i]][[j]]$x_true,
                                      rho = params2[[i]][[j]]$rho, 
                                      alpha = params2[[i]][[j]]$alpha, 
                                      sigma = params2[[i]][[j]]$sigma,
                                      intc = params2[[i]][[j]]$intc)
  }
}
# Save simulated data
saveRDS(params1, 'indpcomp_trueparams1_se.rds')
saveRDS(simdata1, 'indpcomp_simdata1_se.rds')
saveRDS(params2, 'indpcomp_trueparams2_se.rds')
saveRDS(simdata2, 'indpcomp_simdata2_se.rds')
saveRDS(input_data, 'indpcomp_input_data_se.rds')

# Model fit and summary
## declare variable names
x_names <- sprintf('x[%s]', seq(1:N))
rho_f_names <- list()
alpha_f_names <- list()
sigma_f_names <- list()
rho_g_names <- list()
alpha_g_names <- list()
sigma_g_names <- list()
for (i in 1:length(dims)) {
  rho_f_names[[i]] <- sprintf('rho_f[%s]', seq(1:dims[i]))
  alpha_f_names[[i]] <- sprintf('alpha_f[%s]', seq(1:dims[i]))
  sigma_f_names[[i]] <- sprintf('sigma_f[%s]', seq(1:dims[i]))
  rho_g_names[[i]] <- sprintf('rho_g[%s]', seq(1:dims[i]))
  alpha_g_names[[i]] <- sprintf('alpha_g[%s]', seq(1:dims[i]))
  sigma_g_names[[i]] <- sprintf('sigma_g[%s]', seq(1:dims[i]))
}
# parallel for exact GP
out_gp <- list()
out_gp1 <- list()
cores = 2
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_exact = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  set.seed(i)
  out_gp1[[i]] <- list()
  for (j in 1:length(dims)) {
    y1 <- simdata1[[i]][[j]][,1:dims[j]]
    y2 <- simdata2[[i]][[j]][,1:dims[j]]
    x_obs <- input_data[[i]][[j]]$x_obs
    x_true <- input_data[[i]][[j]]$x_true
    gp_fit <- gp_model(gpmodel = gpmodel, 
                       n_obs = N, 
                       dims = dims[j],
                       x_min = x_min,
                       x_max = x_max,
                       y_f = y1, 
                       y_g = y2, 
                       inputs = x_obs, 
                       latent_sd = s_x, 
                       latent_inputs = latent_model, 
                       covfn = covfn_model, 
                       rho_prior = rho_prior_model,
                       ls_param_f = ls_params_model_f,
                       ls_param_g = ls_params_model_g, 
                       msd_param_f = marginal_sd_params_model_f,
                       msd_param_g = marginal_sd_params_model_g,
                       esd_param_f = error_sd_params_model_f,
                       esd_param_g = error_sd_params_model_g,
                       intc_yf = intc_params_model_y1,
                       intc_yg = intc_params_model_y2,
                       is_vary = 1, 
                       is_corr = 1, 
                       iter = 2000, 
                       warmup = 1000, 
                       chains = 1, 
                       cores = 2, 
                       init = 0, 
                       adapt_delta = adapt_delta_model)
    out_gp_x <- compare_summary(model = gp_fit, 
                                variable = x_names, 
                                true_variable = input_data[[i]][[j]]$x_true, 
                                n_obs = N, 
                                m_approx = 'exact', 
                                dims = dims[j], 
                                variable_class = 'x', 
                                model_name = 'gp', 
                                sim_id = i)
    out_gp_rho_f <- compare_summary(model = gp_fit, 
                                  variable = rho_f_names[[j]],
                                  true_variable = params1[[i]][[j]]$rho, 
                                  n_obs = N, 
                                  m_approx = 'exact',
                                  dims = dims[j], 
                                  variable_class = 'rho_f', 
                                  model_name = 'gp', 
                                  sim_id = i)
    out_gp_alpha_f <- compare_summary(model = gp_fit, 
                                    variable = alpha_f_names[[j]],
                                    true_variable = params1[[i]][[j]]$alpha,
                                    n_obs = N, 
                                    m_approx = 'exact',
                                    dims = dims[j], 
                                    variable_class = 'alpha_f', 
                                    model_name = 'gp',
                                    sim_id = i)
    out_gp_sigma_f <- compare_summary(model = gp_fit, 
                                    variable = sigma_f_names[[j]],
                                    true_variable = params1[[i]][[j]]$sigma,
                                    n_obs = N,
                                    m_approx = 'exact',
                                    dims = dims[j], 
                                    variable_class = 'sigma_f', 
                                    model_name = 'gp',
                                    sim_id = i)
    out_gp_rho_g <- compare_summary(model = gp_fit, 
                                    variable = rho_g_names[[j]],
                                    true_variable = params2[[i]][[j]]$rho, 
                                    n_obs = N, 
                                    m_approx = 'exact',
                                    dims = dims[j], 
                                    variable_class = 'rho_g', 
                                    model_name = 'gp', 
                                    sim_id = i)
    out_gp_alpha_g <- compare_summary(model = gp_fit, 
                                      variable = alpha_g_names[[j]],
                                      true_variable = params2[[i]][[j]]$alpha,
                                      n_obs = N, 
                                      m_approx = 'exact',
                                      dims = dims[j], 
                                      variable_class = 'alpha_g', 
                                      model_name = 'gp',
                                      sim_id = i)
    out_gp_sigma_g <- compare_summary(model = gp_fit, 
                                      variable = sigma_g_names[[j]],
                                      true_variable = params2[[i]][[j]]$sigma,
                                      n_obs = N,
                                      m_approx = 'exact',
                                      dims = dims[j], 
                                      variable_class = 'sigma_g', 
                                      model_name = 'gp',
                                      sim_id = i)
    out_gp1[[i]][[j]] <- rbind(out_gp_x, out_gp_rho_f, out_gp_alpha_f, out_gp_sigma_f,
                               out_gp_rho_g, out_gp_alpha_g, out_gp_sigma_g)
  }
  out_gp[[i]] <- rbindlist(out_gp1[[i]])
}
compare_exact <- rbindlist(model_fit_exact)

# parallel for HSGP
out_hsgp <- list()
out_hsgp1 <- list()
cores = 2
cl <- parallel::makeCluster(cores, type="PSOCK")
doParallel::registerDoParallel(cl)
# model fit list
model_fit_hsgp = foreach(i = 1:n_sim) %dopar% {
  library(rstan)
  library(posterior)
  library(data.table)
  set.seed(i)
  out_hsgp1[[i]] <- list()
    for (j in 1:length(dims)) {
      y1 <- simdata1[[i]][[j]][,1:dims[j]]
      y2 <- simdata1[[i]][[j]][,1:dims[j]]
      x_obs <- input_data[[i]][[j]]$x_obs
      x_true <- input_data[[i]][[j]]$x_true
      hsgp_fit <- hsgp_model(hsgpmodel = hsgpmodel,
                             n_obs = N, 
                             dims = dims[j],
                             m_obs_f = min_m_f,
                             m_obs_g = min_m_g,
                             x_min = x_min,
                             x_max = x_max,
                             c = 5/4,
                             y_f = y1, 
                             y_g = y2, 
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
                             intc_yf = intc_params_model_y1,
                             intc_yg = intc_params_model_y2,
                             is_vary = 1, 
                             is_corr = 1, 
                             iter = 2000, 
                             warmup = 1000, 
                             chains = 2, 
                             cores = 1, 
                             init = 0, 
                             adapt_delta = adapt_delta_model)
      out_hsgp_x <- compare_summary(model = hsgp_fit, 
                                    variable = x_names, 
                                    true_variable = input_data[[i]][[j]]$x_true, 
                                    n_obs = N,
                                    m_approx = 'hsgp', 
                                    dims = dims[j],
                                    variable_class = 'x', 
                                    model_name = 'hsgp', 
                                    sim_id = i)
      out_hsgp_rho_f <- compare_summary(model = hsgp_fit,
                                      variable = rho_f_names[[j]], 
                                      true_variable = params1[[i]][[j]]$rho,
                                      n_obs = N, 
                                      m_approx = 'hsgp', 
                                      dims = dims[j], 
                                      variable_class = 'rho_f', 
                                      model_name = 'hsgp',
                                      sim_id = i)
      out_hsgp_alpha_f <- compare_summary(model = hsgp_fit, 
                                        variable = alpha_f_names[[j]], 
                                        true_variable = params1[[i]][[j]]$alpha, 
                                        n_obs = N, 
                                        m_approx = 'hsgp', 
                                        dims = dims[j], 
                                        variable_class = 'alpha_f', 
                                        model_name = 'hsgp', 
                                        sim_id = i)
      out_hsgp_sigma_f <- compare_summary(model = hsgp_fit,
                                        variable = sigma_f_names[[j]], 
                                        true_variable = params1[[i]][[j]]$sigma, 
                                        n_obs = N, 
                                        m_approx = 'hsgp', 
                                        dims = dims[j], 
                                        variable_class = 'sigma_f',
                                        model_name = 'hsgp', 
                                        sim_id = i)
      out_hsgp_rho_g <- compare_summary(model = hsgp_fit,
                                        variable = rho_g_names[[j]], 
                                        true_variable = params2[[i]][[j]]$rho,
                                        n_obs = N, 
                                        m_approx = 'hsgp', 
                                        dims = dims[j], 
                                        variable_class = 'rho_g', 
                                        model_name = 'hsgp',
                                        sim_id = i)
      out_hsgp_alpha_g <- compare_summary(model = hsgp_fit, 
                                          variable = alpha_g_names[[j]], 
                                          true_variable = params2[[i]][[j]]$alpha, 
                                          n_obs = N, 
                                          m_approx = 'hsgp', 
                                          dims = dims[j], 
                                          variable_class = 'alpha_g', 
                                          model_name = 'hsgp', 
                                          sim_id = i)
      out_hsgp_sigma_g <- compare_summary(model = hsgp_fit,
                                          variable = sigma_g_names[[j]], 
                                          true_variable = params2[[i]][[j]]$sigma, 
                                          n_obs = N, 
                                          m_approx = 'hsgp', 
                                          dims = dims[j], 
                                          variable_class = 'sigma_g',
                                          model_name = 'hsgp', 
                                          sim_id = i)
      out_hsgp1[[i]][[j]] <- rbind(out_hsgp_x, out_hsgp_rho_f, out_hsgp_alpha_f, out_hsgp_sigma_f,
                                   out_hsgp_rho_g, out_hsgp_alpha_g, out_hsgp_sigma_g)
  }
  out_hsgp[[i]] <- rbindlist(out_hsgp1[[i]])
}
compare_hsgp <- rbindlist(model_fit_hsgp)

## Bind to long table
compare_table <- rbind(compare_exact, compare_hsgp)

## comparison table
compare_table$data_id = paste0(compare_table$sim_id, '_', 
                               compare_table$m, '_',
                               compare_table$n, '_',
                               compare_table$d)
# Save simulation study results
saveRDS(compare_table, 'indpcompgp_simout_se.rds')
