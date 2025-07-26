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

# Import results from exact and HSGPs
compare_table <- readRDS('simulation results/id_hsgps_simout_sameintc.rds')
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)
compare_x <- subset(compare_table, class == 'x')
levels(compare_x$m)
# Change factor labels according to different simulation studies
compare_x$m <- ordered(compare_x$m)
levels(compare_x$m) <- c("derivgp", "idgp", "idhsgp", "igp", "ihsgp")
compare_x$m <- factor(compare_x$m, levels = c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))
# Create indicator for each simulation conditions
compare_x$subset_id <- paste0(compare_x$n,'_',compare_x$m,'_',compare_x$d)
# Separate results accordingly
list_df <- split(compare_x, compare_x$subset_id)
# Compute log gamma scores for each condition
gamma_stat <- list()
for(i in 1:15){
  gamma_stat[[i]] <- rep(NA, 20)
    for(j in 1:20){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('x','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
}
# Create long table formate for the log gamma scores
gamma_stat_long <- list()
for(i in 1:15){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
}
# Combine the log gamma scores
log_gamma_stat <- rbindlist(gamma_stat_long)
# Clean the results (if required)
log_gamma_stat$log_gamma = ifelse(is.infinite(log_gamma_stat$log_gamma),
                                  -1000, log_gamma_stat$log_gamma)
# Compute log gamma - critical value
log_gamma_stat$log_gamma_diff <- log_gamma_stat$log_gamma - log(SBC:::adjust_gamma(N = 50, K = 1000, L = 1))
# Subset the results for unique values
log_gamma_data <- subset(log_gamma_stat, sim_id==1) 

# Fit summary models to log gamma scores
m_log_gamma <- brm(bf(log_gamma_diff ~ (1 + m) * d + s(true_value, by = m),
                        sigma ~ (1 + m) * d + s(true_value, by = m)),
                     data = log_gamma_data, chains = 2, cores = 2, file_refit = 'on_change')

# Extract results as conditional effects
m_log_gamma_eff <- conditional_effects(m_log_gamma, effects = 'm', 
                                         conditions = make_conditions(m_log_gamma, 'd'),
                                         resolution = 300)

m_log_gamma_eff_s <- conditional_effects(m_log_gamma, effects = 'true_value:m', 
                                           conditions = make_conditions(m_log_gamma, 'd'),
                                           resolution = 300)#, method = 'posterior_predict')

cols <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# Generate plots
# errorbar plots
label_outdims <- c('D = 5','D = 10','D = 20')
# Check and change labels according to the number of basis functions for HSGPs
label_models <- c("derivgp", "igp", "idgp", "ihsgp", "idhsgp")
df_log_gamma_eff <- as.data.frame(m_log_gamma_eff$m)
levels(df_log_gamma_eff$cond__) <- label_outdims
levels(df_log_gamma_eff$effect1__) <- label_models
p_log_gamma_eff <- ggplot(df_log_gamma_eff, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~(differenced))) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(a)')

# Spline plots
df_log_gamma_eff_s <- as.data.frame(m_log_gamma_eff_s$`true_value:m`)
levels(df_log_gamma_eff_s$cond__) <- label_outdims
levels(df_log_gamma_eff_s$effect2__) <- label_models

p_log_gamma_eff_s <- ggplot(data = df_log_gamma_eff_s, aes(x = effect1__, y = estimate__, 
                                                               colour = effect2__, fill = effect2__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_ribbon(aes(ymin = df_log_gamma_eff_s$lower__, ymax = df_log_gamma_eff_s$upper__), alpha = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  facet_wrap(~cond__) +
  labs(x = 'True value', y = expression(log~gamma~(differenced)), colour = 'Models', fill = 'Models') +
  #guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6),
        axis.title.x = element_blank(), axis.title.y = element_blank()) +# legend.position = 'none') +
  scale_colour_manual(values = c("#CC79A7", "#E69F00", "#009E73", "#F0E442", "#D55E00")) + 
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#009E73", "#F0E442", "#D55E00")) #+ ggtitle('(a)')


# Combine plots
p_log_gamma_eff <- (p_log_gamma_eff + p_log_gamma_eff_s) + plot_layout(axis_titles = 'collect')

ggsave('id_hsgps_log_gamma_eff.pdf',
       p_log_gamma_eff,
       dpi = 300,
       width = 80,
       height = 20,
       units = 'cm')

