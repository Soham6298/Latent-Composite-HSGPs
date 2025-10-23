## This script generates the case study results and figures post data analysis
## The objects needs to be obtained by running the case study analysis script

#libraries
library(rstan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)
library(patchwork)
# Source functions
source('indcompgpfns.R')
# Read data for case study
obs_data <- read.csv('case study data/obs.csv')
# Read draws object from case study analysis
us5 <- readRDS('IndptGPs_casestudy/spl_uspl_draws_fullgene.R')
vs5 <- readRDS('IndptGPs_casestudy/spl_velo_draws_fullgene.R')
N <- nrow(obs_data)
# Summarise draws object
x_names <- sprintf('x[%s]', seq(1:N))
# Unspliced+Spliced
us5_x <- subset_draws(us5, variable = x_names)
us5_summary <- summarise_draws(us5_x)
us5_summary$obs_t <- obs_data$numerical_prior 
us5_summary$type <- as.factor(obs_data$celltype)
us5_summary$prior_name <- 'Prior time'
us5_summary$posterior_name <- 'Unspliced + Spliced'
us5_summary <- us5_summary %>% 
  mutate(type = recode(type,'Blood progenitors 1' = 'BP1', 'Blood progenitors 2' = 'BP2',
                       "Erythroid1" = 'E1', "Erythroid2" = 'E2', "Erythroid3" = 'E3'))
# Velocity+Spliced
vs5_x <- subset_draws(vs5, variable = x_names)
vs5_summary <- summarise_draws(vs5_x)
vs5_summary$obs_t <- obs_data$numerical_prior 
vs5_summary$type <- as.factor(obs_data$celltype)
vs5_summary$posterior_name <- 'Spliced + Velocity'
vs5_summary <- vs5_summary %>% 
  mutate(type = recode(type,'Blood progenitors 1' = 'BP1', 'Blood progenitors 2' = 'BP2',
                       "Erythroid1" = 'E1', "Erythroid2" = 'E2', "Erythroid3" = 'E3'))
# Plots
p_prior <- ggplot(us5_summary, aes(x = type, y = obs_t)) +
  theme_bw(base_size=40,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  geom_jitter(height = 0, width = 0.1) +
  facet_wrap(~prior_name) +
  labs(x = 'Cell type', y = 'Ordering') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)) + ggtitle('(a)')

p_latentx <- ggplot(us5_summary, aes(x = type, y = mean)) +
  theme_bw(base_size=40,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  geom_jitter(height = 0, width = 0.1) +
  facet_wrap(~posterior_name) +
  labs(x = 'Cell type', y = 'Ordering') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)) + ggtitle('(b)')

p_latentx_velo <- ggplot(vs5_summary, aes(x = type, y = mean)) +
  theme_bw(base_size=40,
           base_family = 'Times') +
  geom_violin(linewidth = 1, position = position_dodge(width = 0.7)) +
  geom_jitter(height = 0, width = 0.1) +
  facet_wrap(~posterior_name) +
  labs(x = 'Cell type', y = 'Ordering') +
  guides(fill = 'none') +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1)) + ggtitle('(c)')

casestudyplot <- p_prior + p_latentx + p_latentx_velo + plot_layout(axis_titles = 'collect')
ggsave('case_study_14gene.pdf',
       casestudyplot,
       dpi = 300,
       width = 60,
       height = 20,
       units = 'cm')

# Scatterplots between the prior and posterior times
p_scatter_us <- ggplot(us5_summary, aes(x = obs_t, y = mean)) +
  theme_bw(base_size=40,
           base_family = 'Times') +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = '#D55E00', linewidth = 1.5, linetype = 'dashed') +
  facet_wrap(~posterior_name) +
  labs(x = 'Prior time', y = 'Posterior time') +
  guides(fill = 'none') +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = 'none') + ggtitle('(a)')

p_scatter_vs <- ggplot(vs5_summary, aes(x = obs_t, y = mean)) +
  theme_bw(base_size=40,
           base_family = 'Times') +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = '#D55E00', linewidth = 1.5, linetype = 'dashed') +
  facet_wrap(~posterior_name) +
  labs(x = 'Prior time', y = 'Posterior time') +
  guides(fill = 'none') +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = 'none') + ggtitle('(b)')

us_vs <- data.frame(posterior_us = us5_summary$mean, posterior_vs = vs5_summary$mean, type = us5_summary$type, scattername = 'Comparing both posteriors')

p_scatter_us_vs <- ggplot(us_vs, aes(x = posterior_us, y = posterior_vs)) +
  theme_bw(base_size=40,
           base_family = 'Times') +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = '#D55E00', linewidth = 1.5, linetype = 'dashed') +
  facet_wrap(~scattername) +
  labs(x = 'Unspliced+Spliced', y = 'Spliced+Velocity') +
  guides(fill = 'none') +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) + ggtitle('(c)')

casestudyscatter <- p_scatter_us + p_scatter_vs + p_scatter_us_vs + plot_layout(axis_titles = 'collect', guides = 'collect')
ggsave('case_study_14gene_scatter.pdf',
       casestudyscatter,
       dpi = 300,
       width = 60,
       height = 20,
       units = 'cm')