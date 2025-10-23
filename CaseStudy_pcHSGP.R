#libraries
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)
# Source functions
source('indcompgpfns.R')
ihsgpmodel <- stan_model('pcHSGP.stan')
# Read data
spl <- read.csv('casestudydata/spliced_layer.csv')
uspl <- read.csv('casestudyata/unspliced_layer.csv')
obs <- read.csv('casestudydata/obs.csv')
# Data cleanup
spl$type <- as.factor(obs$celltype)
spl$time <- as.factor(obs$numerical_prior)
uspl$type <- as.factor(obs$celltype)
uspl$time <- as.factor(obs$numerical_prior)
# Arrange in concatenated format
# data for model
y_spl <- spl[,-c(1, 16, 17)] 
y_uspl <- uspl[,-c(1, 16, 17)]
t <- as.numeric(as.character(spl$time))
## Standardize data
y_spl <- as.data.frame(scale(y_spl))
y_uspl <- as.data.frame(scale(y_uspl))
colMeans(y_spl)
colMeans(y_uspl)
apply(y_spl, 2, sd)
apply(y_uspl, 2, sd)
# Use only for the short case study
#gene_list <- c("Blvrb","Rbms2","Skap1","Smim1","Yipf5")
#y_spl <- subset(y_spl, select = gene_list)
#y_uspl <- subset(y_uspl, select = gene_list)
# set for prior measurement SD for latent x
s_x <- 0.1
# Number of cells (sample size)
N <- nrow(y_spl)
# Number of genes (output dimensions)
D <- ncol(y_spl)
# Model specifications
true_x_max <- 1
true_x_min <- 0
rho_prior_model <- 0 # 0 = normal; 1 = invgamma;
# Specify priors for the hyperparameters (check if they differ from the simulation setting)

marginal_sd_params_grad_model <-  c(0.5, 0.1)
error_sd_params_grad_model <- c(0.5, 0.1)
marginal_sd_params_model <- c(0.5, 0.1)
error_sd_params_model <- c(0.5, 0.1)
ls_params_model <- c(0.3, 0.1) 
intc_params_model_y <- c(0, 0.1)#c(mean(spl_means), sd(spl_means))
intc_params_model_yderiv <- c(0, 0.1)
adapt_delta_model <- 0.95
# Model fit
# SE
min_m <- 10
ihsgp_fit <- ihsgp_model(hsgpmodel = ihsgpmodel,
                         n_obs = N, 
                         dims = D,
                         m_obs_f = min_m,
                         m_obs_g = min_m,
                         x_min = true_x_min,
                         x_max = true_x_max,
                         c = 5/4,
                         y_f = y_spl, 
                         y_g = y_uspl, 
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
                         chains = 4, 
                         cores = 4, 
                         init = 0, 
                         adapt_delta = adapt_delta_model)

# Save draws output
ihsgp_draws <- as_draws_df(ihsgp_fit)
saveRDS(ihsgp_draws,'spl_uspl_draws_fullgene.R')
