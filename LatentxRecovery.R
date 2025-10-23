## We perform model evaluation for exact and HSGPs based on latent input estimation from simulation studies
## Need to load the simulation result by running the simulation study code or the available results in the folder
## Take care with the ordering of model labels from the sim results. It is crucial to reproduce the figures. 

library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(brms)
library(lemon)
library(grid)
library(gtable)
library(SBC)
library(data.table)
# Import functions
source('indcompgpfns.R')

# Import results from exact, deriv and HSGPs
compare_table <- readRDS('simulation results/pcGP_data_scenario.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
compare_table$model_name <- as.factor(compare_table$model_name)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)
compare_x <- subset(compare_table, class == 'x')
levels(compare_x$m)

# Fit summary models
formula_rmse <- bf(rmse ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
                   sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

m_x_rmse <- brm(formula_rmse, data = compare_x, chains = 2, cores = 2, file_refit = 'on_change')

## Extract summary results as conditional eff data
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

m_rmse_eff <- conditional_effects(m_x_rmse, effects = 'm', 
                                    conditions = make_conditions(m_x_rmse, 'd'),
                                    resolution = 300)
m_rmse_eff_s <- conditional_effects(m_x_rmse, effects = 'true_value:m', 
                                      conditions = make_conditions(m_x_rmse, 'd'),
                                      resolution = 300)

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

# Prepare plots
label_outdims <- c('D = 5','D = 10','D = 20')
# Check and change labels according to the comparison models 
label_models <- c('pcGP', 'pcHSGP') #c('sHSGP', 'sdHSGP', 'pcHSGP', 'pdHSGP')#c('dGP','pcGP', 'pdGP', 'pcHSGP', 'pdHSGP')#
# Posterior bias plots
df_rmse_eff <- as.data.frame(m_rmse_eff$`m`)
# Reorder lables to match the label_models order (comment out for pcGP data scenario)
#df_rmse_eff$effect1__ <- factor(df_rmse_eff$effect1__, levels = c('obs_hsgp', 'deriv_hsgp','ihsgp','idhsgp'))#c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))#
levels(df_rmse_eff$cond__) <- label_outdims
levels(df_rmse_eff$effect1__) <- label_models
p_rmse_eff <- ggplot(df_rmse_eff, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size = 50,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'RMSE') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('(a)')

df_rmse_eff_s <- as.data.frame(m_rmse_eff_s$`true_value:m`)
# Reorder lables to match the label_models order (comment out for pcGP data scenario)
#df_rmse_eff_s$effect2__ <- factor(df_rmse_eff_s$effect2__, levels = c('obs_hsgp', 'deriv_hsgp', 'ihsgp', 'idhsgp'))#c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))#
levels(df_rmse_eff_s$cond__) <- label_outdims
levels(df_rmse_eff_s$effect2__) <- label_models
p_rmse_eff_s <- ggplot(df_rmse_eff_s, aes(x = effect1__, y = estimate__, 
                                          colour = effect2__, fill = effect2__)) +
  theme_bw(base_size = 50,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_rmse_eff_s$lower__, ymax = df_rmse_eff_s$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = 'RMSE', colour = 'Models', fill = 'Models') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), legend.position = 'bottom') +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + 
  scale_fill_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + ggtitle('(b)')

# Combine the plots
p_latentx_eff <- (p_rmse_eff + p_rmse_eff_s) + plot_layout(axis_titles = 'collect')

ggsave('igp_data_compare_latentx.pdf',
       p_latentx_eff,
       dpi = 300,
       width = 80,
       height = 30,
       units = 'cm')

