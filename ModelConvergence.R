# We show and check MCMC convergence for all fitted models across the simulated data scenarios
# We check Rhats, Bulk-ESS and Tail-ESS for fitted models
## Need to load the simulation result by running the simulation study code or the available results in the folder
## Take care with the ordering of model labels from the sim results. It is crucial to reproduce the figures. 

# Libraries
library(dplyr)
library(tidyverse)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(lemon)
library(grid)
library(gtable)
library(gridExtra)
## Set variable names, sample points and number of trials
# Source functions
source('indcompgpfns.R')
simdata_out <- readRDS('simulation results/pcGP_data_scenario.rds')  #change according to output file name from SimStudy
# Change according to simulation scenario
fix_plot_labels <- c('pcGP', 'pcHSGP') #c('sHSGP', 'sdHSGP','pcHSGP', 'pdHSGP')#c('dGP','pcGP', 'pdGP', 'pcHSGP', 'pdHSGP') 
simdata_out$d <- as.factor(simdata_out$d)
simdata_out$m <- as.factor(simdata_out$m)
simdata_out$m <- ordered(simdata_out$m)
levels(simdata_out$m)
# Use to reorder for dGP data scenarios
#simdata_out$m <- ordered(simdata_out$m, levels = c('obs_hsgp', 'deriv_hsgp','ihsgp','idhsgp'))#c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))#
simdata_out$n <- as.factor(simdata_out$n)
simdata_out$sim_id <- as.factor(simdata_out$sim_id)
levels(simdata_out$m)
# Convergence for latent x
simout_x <- subset(simdata_out, class == 'x')
convsummary <- simout_x %>%
  group_by(data_id) %>%
  summarise(sim_id = first(sim_id),
            n = first(n),
            m = first(m),
            d = first(d),
            mrhat = mean(rhat),
            mbess = mean(bess),
            mtess = mean(tess),
            model_name = first(model_name))
convsummary$rhat_name <- 'Rhat'
convsummary$bess_name <- 'Bulk-ESS'
convsummary$tess_name <- 'Tail-ESS'
# Rhat plot
x_rhat_summary_plot <- ggplot(convsummary, aes(x = d, y = mrhat, colour = m)) + 
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_boxplot(linewidth = 1) +
  facet_wrap(~rhat_name) +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) +
  ggtitle('(a)')
# Bulk-ESS plot
x_ess_bulk_summary_plot <- ggplot(convsummary, aes(x = d, y = mbess, colour = m)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_boxplot(linewidth = 1) +
  facet_wrap(~bess_name) +
  scale_y_log10() +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))
# Tail-ESS plot
x_ess_tail_summary_plot <- ggplot(convsummary, aes(x = d, y = mtess, colour = m)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_boxplot(linewidth = 1) +
  facet_wrap(~tess_name) +
  scale_y_log10() +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))
# combine diags for latent x
p_latentx <- x_rhat_summary_plot + x_ess_bulk_summary_plot + x_ess_tail_summary_plot + 
  plot_layout(axis_titles = 'collect', guides = 'collect') & theme(axis.title = element_blank())

## Convergence for GP hyperparameters
simout_x <- subset(simdata_out, class != 'x')
convsummary <- simout_x %>%
  group_by(data_id) %>%
  summarise(sim_id = first(sim_id),
            n = first(n),
            m = first(m),
            d = first(d),
            mrhat = mean(rhat),
            mbess = mean(bess),
            mtess = mean(tess),
            model_name = first(model_name))
levels(convsummary$m)
convsummary$rhat_name <- 'Rhat'
convsummary$bess_name <- 'Bulk-ESS'
convsummary$tess_name <- 'Tail-ESS'
# Rhat plot
pars_rhat_summary_plot <- ggplot(convsummary, aes(x = d, y = mrhat, colour = m)) + 
  theme_bw(base_size = 35, base_family = 'Times') + 
  geom_boxplot(linewidth = 1) +
  facet_wrap(~rhat_name) +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) +
  ggtitle('(b)')
# Bulk-ESS plot
pars_ess_bulk_summary_plot <- ggplot(convsummary, aes(x = d, y = mbess, colour = m)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_boxplot(linewidth = 1) +
  facet_wrap(~bess_name) +
  scale_y_log10() +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))
# Tail-ESS plot
pars_ess_tail_summary_plot <- ggplot(convsummary, aes(x = d, y = mtess, colour = m)) + 
  theme_bw(base_size = 35, base_family = 'Times') +
  geom_boxplot(linewidth = 1) +
  facet_wrap(~tess_name) +
  scale_y_log10() +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))
# Combine diags for hyperparameters
p_hyperpars <- pars_rhat_summary_plot + pars_ess_bulk_summary_plot + pars_ess_tail_summary_plot + 
  plot_layout(axis_titles = 'collect', guides = 'collect')

# Combine all figures
p_convdiag <- p_latentx / p_hyperpars + plot_layout(axis_titles = 'collect', guides = 'collect') & 
  theme(axis.title.y = element_blank())
ggsave('pcGP_data_valid.pdf',
       p_convdiag,
       dpi = 300,
       width = 50,
       height = 30,
       units = 'cm')

