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
compare_table <- readRDS('simulation results/idhsgp_single_simout_n100_extra.rds')
n_samples <- compare_table$n[1]
compare_table$sim_id <- as.factor(compare_table$sim_id)
compare_table$n <- as.factor(compare_table$n)
compare_table$m <- as.factor(compare_table$m)
str(compare_table$m)
compare_table$d <- as.factor(compare_table$d)
compare_table$data_id <- as.factor(compare_table$data_id)
compare_x <- subset(compare_table, class == 'rho')
levels(compare_x$m)
n_models <- length(levels(compare_x$m))-1
n_d <- length(levels(compare_table$d))
#n_samples <- compare_table$n[1]
# Change factor labels according to different simulation studies
compare_x$m <- ordered(compare_x$m)
#levels(compare_x$m) <- c("derivgp", "idgp", "idhsgp", "igp", "ihsgp")
compare_x$m <- ordered(compare_x$m, levels = c('obs_hsgp', 'deriv_hsgp', 'ihsgp','idhsgp'))#factor(compare_x$m, levels = c("derivgp", "igp", "idgp", "ihsgp", "idhsgp"))
# Create indicator for each simulation conditions
compare_x$subset_id <- paste0(compare_x$n,'_',compare_x$m,'_',compare_x$d)
max_rank = max(compare_x$ranks)
# Separate results accordingly
list_df <- split(compare_x, compare_x$subset_id)
# Compute log gamma scores for each condition
gamma_stat <- list()
for(i in 1:(n_models*n_d)){
  if(list_df[[i]]$d[1] == '5'){
    gamma_stat[[i]] <- rep(NA, 5)
    for(j in 1:5){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('rho','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
  }else if(list_df[[i]]$d[1] == '10'){
    gamma_stat[[i]] <- rep(NA, 10)
    for(j in 1:10){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('rho','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
  }else if(list_df[[i]]$d[1] == '20'){
    gamma_stat[[i]] <- rep(NA, 20)
    for(j in 1:20){
      testrank_subset <- subset(list_df[[i]],  pars == paste0('rho','[',j,']'))
      gamma_stat[[i]][j] <- log_gamma_statistic(testrank_subset$rank, max_rank = 1000)
    }
  }
}
# Create long table formate for the log gamma scores
gamma_stat_long <- list()
for(i in 1:(n_models*n_d)){
  if(list_df[[i]]$d[1] == '5'){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
  }else if(list_df[[i]]$d[1] == '10'){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
  }else if(list_df[[i]]$d[1] == '20'){
    gamma_stat_long[[i]] <- data.frame(log_gamma = gamma_stat[[i]], 
                                       n = list_df[[i]]$n, 
                                       m = list_df[[i]]$m, 
                                       d = list_df[[i]]$d, 
                                       sim_id = list_df[[i]]$sim_id,
                                       subset_id = list_df[[i]]$subset_id,
                                       true_value = list_df[[i]]$true_value)
  }
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
# Separate log gamma scores per sample size
log_gamma_data_5 <- subset(log_gamma_data, d == 5)
log_gamma_data_10 <- subset(log_gamma_data, d == 10)
log_gamma_data_20 <- subset(log_gamma_data, d == 20)

# Fit summary models to log gamma scores
m_log_gamma5 <- brm(bf(log_gamma_diff ~ (1 + m), #+ s(true_value, by = m),
                        sigma ~ (1 + m)), #+ s(true_value, by = m)),
                     data = log_gamma_data_5, chains = 2, cores = 2, file_refit = 'on_change')

m_log_gamma10 <- brm(bf(log_gamma_diff ~ (1 + m), #+ s(true_value, by = m),
                        sigma ~ (1 + m)), #+ s(true_value, by = m)),
                     data = log_gamma_data_10, chains = 2, cores = 2, file_refit = 'on_change')

m_log_gamma20 <- brm(bf(log_gamma_diff ~ (1 + m), #+ s(true_value, by = m),
                         sigma ~ (1 + m)), #+ s(true_value, by = m)),
                      data = log_gamma_data_20, chains = 2, cores = 2, file_refit = 'on_change')

# Extract results as conditional effects
m_log_gamma_eff5 <- conditional_effects(m_log_gamma5, effects = 'm',
                                         resolution = 300)
m_log_gamma_eff10 <- conditional_effects(m_log_gamma10, effects = 'm',
                                         resolution = 300)
m_log_gamma_eff20 <- conditional_effects(m_log_gamma20, effects = 'm',
                                          resolution = 300)

cols <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# Generate plots

# Check and change labels according to the number of basis functions for HSGPs
label_models <- c('HSGP(SE)', 'HSGP(dSE)','id-HSGP')
df_log_gamma_eff5 <- as.data.frame(m_log_gamma_eff5$m)
levels(df_log_gamma_eff5$effect1__) <- label_models
# Plot for D = 5
p_log_gamma_eff5 <- ggplot(df_log_gamma_eff5, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  #facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~- "critical value")) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6)) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(a) D = 5')

# Plot for D = 10
df_log_gamma_eff10 <- as.data.frame(m_log_gamma_eff10$m)
levels(df_log_gamma_eff10$effect1__) <- label_models
p_log_gamma_eff10 <- ggplot(df_log_gamma_eff10, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  #facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~- "critical value")) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6), axis.title.y = element_blank()) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(b) D = 10')

# Plot for D = 20
df_log_gamma_eff20 <- as.data.frame(m_log_gamma_eff20$m)
levels(df_log_gamma_eff20$effect1__) <- label_models
p_log_gamma_eff20 <- ggplot(df_log_gamma_eff20, aes(x = effect1__, y = estimate__)) +
  theme_bw(base_size=35,
           base_family = 'Times') +
  geom_point(size = 3.5 ,
             position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.5,
                linewidth = 1.0,
                position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0, colour = "#0072B2", linetype = 'dashed', linewidth = 1.5) +
  #facet_wrap(~cond__) +
  labs(x = 'Models', y = expression(log~gamma~(differenced))) +
  guides(fill = 'none') + 
  theme(axis.ticks = element_line(linewidth = 3), 
        axis.text.x = element_text(angle = 35, vjust = 0.7, hjust = 0.6), axis.title.y = element_blank()) +
  scale_colour_manual(values = c("#000000")) + ggtitle('(c) D = 20')


# Combine plots
p_log_gamma_eff <- (p_log_gamma_eff5 + p_log_gamma_eff10 + p_log_gamma_eff20) + plot_layout(axis_titles = 'collect')
ggsave('idhsgps_n100_log_gamma_eff_rho.pdf',
       p_log_gamma_eff,
       dpi = 300,
       width = 50,
       height = 20,
       units = 'cm')

