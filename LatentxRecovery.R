## We perform model evaluation for exact and HSGPs based on latent input estimation from simulation studies
## Need to load the output dataframe from the simulation study code before generating figures here

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
compare_table <- readRDS('simulation results/comp_deriv_se_simout.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)
compare_x <- subset(compare_table, class == 'x')
#compare_x <- subset(compare_x, m != 'exact')
compare_rho <- subset(compare_table, class == 'rho')
compare_alpha <- subset(compare_table, class == 'alpha')
compare_sigma <- subset(compare_table, class == 'sigma')
#levels(compare_x$m)
# Change factor labels according to different simulation studies
#compare_x$m <- factor(compare_x$m, levels = c('exact', '22', '26', '30', 'pyroVI'))

# Fit summary models
formula_rmse <- bf(rmse ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
                   sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

#formula_sd <- bf(sd ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m),
#                 sigma ~ (1 + m) * d + (1 + m | data_id) + s(true_value, by = m))

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
# Check and change labels according to the number of basis functions for HSGPs
label_models <- c('Deriv GP', 'Exact GP', 'HSGP')
# Posterior bias plots
df_rmse_eff <- as.data.frame(m_rmse_eff$`m`)
levels(df_rmse_eff$cond__) <- label_outdims
levels(df_rmse_eff$effect1__) <- label_models
p_rmse_eff <- ggplot(df_rmse_eff, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  annotate('point', x=0.6, y=naive_rmse, colour = '#D55E00', size = 5.5) +
  annotate('text', x = 0.8, y = naive_rmse + 0.02, label = 'Prior', size = 7, colour = '#D55E00') +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'RMSE') +
  #guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000" )) + ggtitle('(a)')

df_rmse_eff_s <- as.data.frame(m_rmse_eff_s$`true_value:m`)
levels(df_rmse_eff_s$cond__) <- label_outdims
levels(df_rmse_eff_s$effect2__) <- label_models
p_rmse_eff_s <- ggplot(df_rmse_eff_s, aes(x = effect1__, y = estimate__, 
                                          colour = effect2__, fill = effect2__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_rmse_eff_s$lower__, ymax = df_rmse_eff_s$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = 'RMSE', colour = 'Models', fill = 'Models') +
  #guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#E69F00", "#009E73", "#CC79A7")) + 
  scale_fill_manual(values = c("#E69F00", "#009E73", "#CC79A7")) + ggtitle('(b)')

# Combine the plots
p_latentx_eff <- (p_rmse_eff + p_rmse_eff_s) + plot_layout(axis_titles = 'collect')

ggsave('comp_se_latentx.pdf',
       p_latentx_eff,
       dpi = 300,
       width = 60,
       height = 20,
       units = 'cm')

