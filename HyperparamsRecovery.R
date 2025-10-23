## We check model performance for all fitted models based on class of model parameter recovery across simulation studies
## Need to load the simulation result by running the simulation study code or the available results in the folder
## Take care with the ordering of model labels from the sim results. It is crucial to reproduce the figures. 
## We only show hyperparameter recovery for the most complicated case of D = 20. The other cases are qualitatively similar

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
source('indcompgpfns.R')
# Summary model
## load simulation study output
compare_table <- readRDS('simulation results/cgp_data_compare.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
# Designate levels properly to set figure legends and guides
levels(compare_table$m)
# Check and change labels according to the number of basis functions for HSGPs
compare_table$m <- ordered(compare_table$m)
compare_table <- compare_table %>% 
  mutate(m = recode(m,'igp' = 'pcGP', 'ihsgp' = 'pcHSGP'))
# Use for dGP data scenario
#compare_table <- compare_table %>% 
#  mutate(m = recode(m, "derivgp" = 'dGP', 'igp' = 'pcGP', 'idgp' = 'pdGP', 'ihsgp' = 'pcHSGP', 'idhsgp' = 'pdHSGP'))
# Use for dGP data scenario N = 100
#compare_table <- compare_table %>% 
#  mutate(m = recode(m, "obs_hsgp" = 'sHSGP', 'deriv_hsgp' = 'sdHSGP', 'ihsgp' = 'pcHSGP', 'idhsgp' = 'pdHSGP'))
compare_table <- compare_table %>% 
  mutate(d = recode(d, "5" = 'D = 5', '10' = 'D = 10', '20' = 'D = 20'))
compare_table <- compare_table %>% 
  mutate(label = case_when(
    class %in% c('rho', 'rho_f', 'rho_g') ~ "Length scale",
    class %in% c('alpha_obs', 'alpha_grad', 'alpha_f', 'alpha_g', 'alpha') ~ "Marginal SD",
    class %in% c('sigma_obs', 'sigma_grad', 'sigma_f', 'sigma_g', 'sigma') ~ "Error SD",
    class == 'x' ~ 'x'
    ))
# Set data for generating figure
compare_table$data_id <- as.factor(compare_table$data_id)
compare_rho <- subset(compare_table, label == 'Length scale' & d == 'D = 20')
compare_alpha <- subset(compare_table, label == 'Marginal SD' & d == 'D = 20')
compare_sigma <- subset(compare_table, label == 'Error SD' & d == 'D = 20')
# Generate figure
p_rho_rmse <- ggplot(compare_rho, aes(x = label, y = rmse, colour = m)) +
  theme_bw(base_size=15,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  facet_wrap(~label) +
  labs(x = '', y = 'RMSE', colour = 'Models') +
  guides(fill = 'none') +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank()) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + ggtitle('')

p_alpha_rmse <- ggplot(compare_sigma, aes(x = label, y = rmse, colour = m)) +
  theme_bw(base_size=15,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  facet_wrap(~label) +
  labs(x = '', y = 'RMSE', colour = 'Models') +
  guides(fill = 'none') +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank()) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + ggtitle('')

p_sigma_rmse <- ggplot(compare_alpha, aes(x = label, y = rmse, colour = m)) +
  theme_bw(base_size=15,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  facet_wrap(~ label) +
  labs(x = '', y = 'RMSE', colour = 'Models') +
  guides(fill = 'none') +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank()) +
  scale_colour_manual(values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) + ggtitle('')

hyperparams_plot <- (p_rho_rmse + p_alpha_rmse + p_sigma_rmse) + plot_layout(axis_titles = 'collect', guides = 'collect')

ggsave('pcGP_data_hyperparams.pdf',
       hyperparams_plot,
       dpi = 300,
       width = 25,
       height = 8,
       units = 'cm')
