library(tidyverse)
library(ggpubr)

#source('utils_plot_inference.R')
source('utils_cluster_alloc.R')

## MERS input files
output_file_likelihood_CI_mers <- '../results/mers/df_inference_10000_central.rds'
file_cluster_alloc_mers <- '../data/mers/cluster_alloc_mers.rds'

## Measles input files
output_file_likelihood_CI_measles <- '../results/measles/df_inference_pacenti_10000_central.rds'
file_cluster_alloc_measles <- '../data/measles/cluster_alloc_measles_pacenti.rds'

## Loading files
cluster_alloc_mers <- readRDS(file_cluster_alloc_mers)
df_CI_mers <- readRDS(output_file_likelihood_CI_mers)

cluster_alloc_measles <- readRDS(file_cluster_alloc_measles)
df_CI_measles <- readRDS(output_file_likelihood_CI_measles)

col_MERS <- 'firebrick'
col_measles <- '#669d1c'

## MERS reproduction number estimate
plt_R_mers <- df_CI_mers %>%  
  filter(param == 'R') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim)) +
  theme_bw() +
  geom_hline(yintercept = 0.47, col = 'gray22', alpha = 0.7) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.29, ymax = 0.8),
            alpha = 0.2, fill = 'gray22') +
  geom_point(col = col_MERS) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), col = col_MERS) +
  scale_y_continuous(limits = c(0., 1.0),
                     name = 'R estimate',
                     expand = expansion(mult = c(0., 0.1))) +
  scale_x_discrete(name = '')  +
  theme(panel.grid.minor = element_blank())

## MERS dispersion parameter estimate
plt_k_mers <- df_CI_mers %>% 
  filter(param == 'k') %>%
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim)) +
  geom_hline(yintercept = 0.26, col = 'gray22', alpha = 0.7) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.09, ymax = 1.24),
            alpha = 0.2, fill = 'gray22') +
  theme_bw() +
  geom_point(col = col_MERS) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), col = col_MERS) +
  scale_y_continuous(trans = 'log',
                     name = 'k estimate',
                     breaks = c(seq(0.001, 0.009, 0.001), seq(0.01, 0.09, 0.01),
                                seq(0.1, 0.9, 0.1), seq(1.0, 10.0, 1.0)),
                     labels = c(0.001, 0.002, rep('', 2), 0.005, rep('', 4),
                                0.01, 0.02, rep('', 2), 0.05, rep('', 4),
                                0.1, 0.2, rep('', 2), 0.5, rep('', 4),
                                1.0, 2.0, rep('', 2), 5.0, rep('', 4), 10.0),
                     expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(ylim = c(0.02, 10)) +
  scale_x_discrete(name = '') +
  theme(panel.grid.minor = element_blank())

## Measles R estimate
plt_R_measles <- df_CI_measles %>% 
  filter(param == 'R') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim))+
  theme_bw() +
  geom_point(col = col_measles) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), col = col_measles) +
  scale_y_continuous(limits = c(0., NA),
                     breaks = seq(0.0, 1.6, 0.2),
                     name = 'R estimate',
                     expand = expansion(mult = c(0., 0.1))) +
  scale_x_discrete(name = '')  +
  theme(panel.grid.minor = element_blank())

## Measles dispersion parameter estimate
plt_k_measles <- df_CI_measles %>% 
  filter(param == 'k') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim)) +
  theme_bw() +
  geom_point(col = col_measles) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), col = col_measles) +
  scale_y_continuous(trans = 'log',
                     name = 'k estimate',
                     breaks = c(seq(0.001, 0.009, 0.001),
                                seq(0.01, 0.09, 0.01),
                                seq(0.1, 0.9, 0.1),
                                seq(1.0, 10.0, 1.0)),
                     labels = c(0.001, 0.002, rep('', 2), 0.005, rep('', 4),
                                0.01, 0.02, rep('', 2), 0.05, rep('', 4),
                                0.1, 0.2, rep('', 2), 0.5, rep('', 4),
                                1.0, 2.0, rep('', 2), 5.0, rep('', 4), 10.0),
                     expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(ylim = c(0.001, 10)) +
  scale_x_discrete(name = '') +
  theme(panel.grid.minor = element_blank())

## Distribution of the size of clusters of identical sequences
plt_size_cluster_measles <- plt_cluster_size_distrib(cluster_alloc_measles, col_measles)
plt_size_cluster_mers <- plt_cluster_size_distrib(cluster_alloc_mers, col_MERS)

plt_mers_measles <- ggarrange(plt_size_cluster_mers + ggtitle('MERS'),
                              plt_R_mers + ggtitle(''),
                              plt_k_mers + ggtitle(''), 
                              plt_size_cluster_measles + ggtitle('Measles'), 
                              plt_R_measles + ggtitle(''),
                              plt_k_measles + ggtitle(''), 
                              ncol = 3, nrow = 2,
                              labels = 'AUTO')

plot(plt_mers_measles)


########## Check the sensitivity of our estimates to the maximum cluster size threshold c_max used in the inference
##### -> No impact here. 
vec_max_cluster_size_inference <- c('10000', '50000', '1e+05')
df_CI_mers <- Reduce('bind_rows', lapply(vec_max_cluster_size_inference, FUN = function(max_cluster_size_inference){
  output_file_likelihood_CI_mers <- paste0('../results/mers/df_inference_', max_cluster_size_inference, '_central.rds')
  readRDS(output_file_likelihood_CI_mers) %>% 
    mutate(max_cluster_size_inference = max_cluster_size_inference)
}))
df_CI_measles <- Reduce('bind_rows', lapply(vec_max_cluster_size_inference, FUN = function(max_cluster_size_inference){
  output_file_likelihood_CI_measles <- paste0('../results/measles/df_inference_pacenti_', max_cluster_size_inference, '_central.rds')
  readRDS(output_file_likelihood_CI_measles) %>% 
    mutate(max_cluster_size_inference = max_cluster_size_inference)
}))

df_CI_mers %>% 
  arrange(param, p_detect) %>% 
  View()

df_CI_measles %>% 
  arrange(param, p_detect) %>% 
  View()
