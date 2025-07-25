# We show and check MCMC convergence for HSGP and the exact GP fitted models for the simulated data scenarios
# We check Rhats, Bulk-ESS and Tail-ESS for fitted models

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
simdata_out <- readRDS('simulation results/id_hsgps_simout.rds')  #change according to output file name from SimStudy
fix_plot_labels <- c('Deriv-GP','i-GP', 'id-GP', 'i-HSGP', 'id-HSGP') # Change according to simulation scenario
simdata_out$d <- as.factor(simdata_out$d)
simdata_out$m <- as.factor(simdata_out$m)
simdata_out$n <- as.factor(simdata_out$n)
#simdata_out$data_id <- as.factor(simdata_out$data_id)
simdata_out$sim_id <- as.factor(simdata_out$sim_id)

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
levels(simout_x$m)
#simout_x$m <- ordered(simout_x$m)
#levels(simout_x$m) <- c("derivgp", "idgp", "idhsgp", "igp", "ihsgp")
simout_x$m <- factor(simout_x$m, levels = c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))
# Rhat plot
x_rhat_summary_plot <- ggplot(convsummary, aes(x = d, y = mrhat, colour = model_name)) + 
  theme_bw(base_size = 30, base_family = 'Times') + 
  geom_boxplot(linewidth = 1) +
  facet_wrap(~rhat_name) +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) +
  ggtitle('(a)')
# Bulk-ESS plot
x_ess_bulk_summary_plot <- ggplot(convsummary, aes(x = d, y = mbess, colour = model_name)) + 
  theme_bw(base_size = 30, base_family = 'Times') +
  geom_boxplot(linewidth = 1) +
  facet_wrap(~bess_name) +
  scale_y_log10() +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))
# Tail-ESS plot
x_ess_tail_summary_plot <- ggplot(convsummary, aes(x = d, y = mtess, colour = model_name)) + 
  theme_bw(base_size = 30, base_family = 'Times') +
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
convsummary$rhat_name <- 'Rhat'
convsummary$bess_name <- 'Bulk-ESS'
convsummary$tess_name <- 'Tail-ESS'
levels(simout_x$m)
#simout_x$m <- ordered(simout_x$m)
#levels(simout_x$m) <- c("derivgp", "idgp", "idhsgp", "igp", "ihsgp")
simout_x$m <- factor(simout_x$m, levels = c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))
# Rhat plot
pars_rhat_summary_plot <- ggplot(convsummary, aes(x = d, y = mrhat, colour = model_name)) + 
  theme_bw(base_size = 30, base_family = 'Times') + 
  geom_boxplot(linewidth = 1) +
  facet_wrap(~rhat_name) +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2")) +
  ggtitle('(b)')
# Bulk-ESS plot
pars_ess_bulk_summary_plot <- ggplot(convsummary, aes(x = d, y = mbess, colour = model_name)) + 
  theme_bw(base_size = 30, base_family = 'Times') +
  geom_boxplot(linewidth = 1) +
  facet_wrap(~bess_name) +
  scale_y_log10() +
  labs(x = 'Output dimensions', y = 'Values', colour = 'Models') + 
  theme(axis.ticks = element_line(linewidth = 3)) +
  scale_colour_manual(labels = fix_plot_labels, values = c("#56B4E9", "#E69F00", "#009E73", "#CC79A7", "#0072B2"))
# Tail-ESS plot
pars_ess_tail_summary_plot <- ggplot(convsummary, aes(x = d, y = mtess, colour = model_name)) + 
  theme_bw(base_size = 30, base_family = 'Times') +
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
ggsave('id_hsgps_valid.pdf',
       p_convdiag,
       dpi = 300,
       width = 50,
       height = 30,
       units = 'cm')

