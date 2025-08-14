#libraries
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
# Source functions
source('indcompgpfns.R')
idhsgpmodel <- stan_model('Stan models/id_hsgp.stan')
ihsgpmodel <- stan_model('Stan models/i_hsgp.stan')
# Read data
spl <- read.csv('case study data/moment_spliced_old_genes.csv')
obs <- read.csv('case study data/observation_old_genes.csv')
velo <- read.csv('case study data/velocity_old_genes.csv')
# Data cleanup
spl$type <- as.factor(obs$celltype)
spl$time <- as.factor(obs$numerical_prior)
velo$type <- as.factor(obs$celltype)
velo$time <- as.factor(obs$numerical_prior)
str(spl)
# Sub-sample data # seed 24551 worked with genes 1,4,7 # seed 651 all genes
set.seed(5445)
subsample_spl <- spl %>%
  group_by(time) %>%
  slice_sample(n = 10) %>%  # sample n rows per group
  ungroup()
subsample_velo <- subsample_spl %>%
  select(X) %>% 
  inner_join(velo, by = "X")
# Arrange in concatenated format
# data for model
y_spl <- subsample_spl[,-c(1, 13, 14)] #(1, 9, 10)
y_velo <- subsample_velo[,-c(1, 13, 14)] #(1, 9, 10)
t <- as.numeric(as.character(subsample_spl$time))
## Standardize data
y_spl <- as.data.frame(scale(y_spl))
y_velo <- as.data.frame(scale(y_velo))
colMeans(y_spl)
colMeans(y_velo)
apply(y_spl, 2, sd)
apply(y_velo, 2, sd)
# Convert df to long format
spl_long <- y_spl %>%
  mutate(time = t) %>%
  pivot_longer(cols = -time, names_to = "D", values_to = "value")
velo_long <- y_velo %>%
  mutate(time = t) %>%
  pivot_longer(cols = -time, names_to = "D", values_to = "value")
# Plot all in one figure with facets
ggplot(spl_long, aes(x = time, y = value)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~ D, scales = "free_y") +
  labs(x = "Exp time", y = "counts") + ggtitle('Spliced data')
ggplot(velo_long, aes(x = time, y = value)) +
  theme_bw() +
  geom_point() +
  facet_wrap(~ D, scales = "free_y") +
  labs(x = "Exp time", y = "counts") + ggtitle('Velocity data')
#y_spl <- y_spl[,-c(2,3,5,6)]
#y_velo <- y_velo[,-c(2,3,5,6)]
# set for prior measurement SD for latent x
s_x <- 0.05
# Number of cells (sample size)
N <- nrow(y_spl)
# Number of genes (output dimensions)
D <- ncol(y_spl)
# Model specifications
true_x_max <- 1
true_x_min <- 0
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)
#lambda <-  abs(mean(colMeans(y_spl / y_velo, na.rm = TRUE))) # derivative scale factor
col_sd_spl <- apply(y_spl, 2, sd, na.rm = TRUE)
col_sd_velo <- apply(y_velo, 2, sd, na.rm = TRUE)
marginal_sd_params_grad_model <-  c(0.5, 0.1)#c(mean(col_sd_velo), sd(col_sd_velo)) # alpha_grad ~ N( , )
error_sd_params_grad_model <- c(0.5, 0.1)#c(mean(col_sd_velo), sd(col_sd_velo))    # sigma_grad ~ N( , )
marginal_sd_params_model <- c(0.5, 0.1)#c(mean(col_sd_spl), sd(col_sd_spl))  # alpha_obs ~ N( , )
error_sd_params_model <- c(0.5, 0.1)#c(mean(col_sd_spl), sd(col_sd_spl)) # sigma_obs ~ N( , )
ls_params_model <- c(0.3, 0.1) # rho ~ N( , ) or InvGamma( , )
spl_means <- colMeans(y_spl)
velo_means <- colMeans(y_velo)
intc_params_model_y <- c(0, 0.1)#c(mean(spl_means), sd(spl_means))
intc_params_model_yderiv <- c(0, 0.1)
adapt_delta_model <- 0.95
# Model fit
# SE
init_per_chain <- list(
  #x_temp = rep(0, N),
  rho_temp = rep(0.3, D),
  alpha_f_temp = rep(0.5, D),
  alpha_g_temp = rep(0.5, D),
  sigma_f_temp = rep(0.5, D),
  sigma_g_temp = rep(0.5, D)
)
init_list <- list(init_per_chain, init_per_chain)
min_m <- 10
hsgp_fit <- idhsgp_model(hsgpmodel = idhsgpmodel,
                         n_obs = N, 
                         dims = D,
                         m_obs_f = min_m,
                         m_obs_g = min_m,
                         x_min = true_x_min,
                         x_max = true_x_max,
                         c = 5/4,
                         y_f = y_spl, 
                         y_g = y_velo, 
                         inputs = t, 
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
                         chains = 2, 
                         cores = 2, 
                         init = init_list, 
                         adapt_delta = adapt_delta_model)
#saveRDS(hsgp_fit,'id_hsgp_test.R')
print(hsgp_fit, pars = c('rho', 'alpha_f', 'alpha_g', 'sigma_f', 'sigma_g'))
traceplot(hsgp_fit, pars = c('x'))#c('rho','alpha_f', 'alpha_g', 'sigma_f', 'sigma_g'))#



#######
ihsgp_fit <- ihsgp_model(hsgpmodel = ihsgpmodel,
                         n_obs = N, 
                         dims = D,
                         m_obs_f = min_m,
                         m_obs_g = min_m,
                         x_min = true_x_min,
                         x_max = true_x_max,
                         c = 5/4,
                         y_f = y_spl, 
                         y_g = y_velo, 
                         inputs = t, 
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
                         chains = 2, 
                         cores = 2, 
                         init = 0, 
                         adapt_delta = adapt_delta_model)
#saveRDS(ihsgp_fit,'i_hsgp_test.R')
print(ihsgp_fit, pars = c('rho_f', 'rho_g', 'alpha_f', 'alpha_g', 'sigma_f', 'sigma_g'))
traceplot(ihsgp_fit, pars = c('x[53]'))#c('rho','alpha_f', 'alpha_g', 'sigma_f', 'sigma_g'))#

