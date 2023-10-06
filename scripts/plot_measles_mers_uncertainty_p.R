library(tidyverse)
library(ggpubr)
library(colorspace)

#source('utils_plot_inference.R')
source('utils_cluster_alloc.R')

## Loading files
output_file_likelihood_CI_mers_central <- '../results/mers/df_inference_10000_central.rds'
output_file_likelihood_CI_mers_lower <- '../results/mers/df_inference_10000_lower.rds'
output_file_likelihood_CI_mers_upper <- '../results/mers/df_inference_10000_upper.rds'
output_file_likelihood_CI_measles_central <- '../results/measles/df_inference_pacenti_10000_central.rds'
output_file_likelihood_CI_measles_lower <- '../results/measles/df_inference_pacenti_10000_lower.rds'
output_file_likelihood_CI_measles_upper <- '../results/measles/df_inference_pacenti_10000_upper.rds'

df_CI_mers_central <- readRDS(output_file_likelihood_CI_mers_central)
df_CI_mers_lower <- readRDS(output_file_likelihood_CI_mers_lower)
df_CI_mers_upper <- readRDS(output_file_likelihood_CI_mers_upper)

df_CI_measles_central <- readRDS(output_file_likelihood_CI_measles_central)
df_CI_measles_lower <- readRDS(output_file_likelihood_CI_measles_lower)
df_CI_measles_upper <- readRDS(output_file_likelihood_CI_measles_upper)

df_CI_measles_uncertainty <- bind_rows(df_CI_measles_central %>% mutate(p_assumption = 'Central'),
                                      df_CI_measles_lower %>% mutate(p_assumption = 'Lower bound'),
                                      df_CI_measles_upper %>% mutate(p_assumption = 'Upper bound')) %>% 
  mutate(p_assumption = factor(p_assumption, levels = c('Lower bound', 'Central', 'Upper bound')))

df_CI_mers_uncertainty <- bind_rows(df_CI_mers_central %>% mutate(p_assumption = 'Central'),
                                    df_CI_mers_lower %>% mutate(p_assumption = 'Lower bound'),
                                    df_CI_mers_upper %>% mutate(p_assumption = 'Upper bound')) %>% 
  mutate(p_assumption = factor(p_assumption, levels = c('Lower bound', 'Central', 'Upper bound')))


## Define colours for plot
col_mers <- 'firebrick'
col_mers_lower <- colorspace::lighten(col_mers, 0.7)
col_mers_upper <- colorspace::darken(col_mers, 0.7)
col_measles <- '#669d1c'
col_measles_lower <- colorspace::lighten(col_measles, 0.7)
col_measles_upper <- colorspace::darken(col_measles, 0.7)

plt_R_measles_uncertainty <- df_CI_measles_uncertainty %>% 
  filter(param == 'R') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim, colour = as.factor(p_assumption)))+
  theme_bw() +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), position = position_dodge(0.5)) +
  scale_y_continuous(limits = c(0., NA),
                     breaks = seq(0.0, 1.6, 0.2),
                     name = 'R estimate',
                     expand = expansion(mult = c(0., 0.1))) +
  scale_x_discrete(name = '')  +
  scale_color_manual(values = c(col_measles_lower, col_measles, col_measles_upper),
                     name = 'Assumption for p') +
  theme(panel.grid.minor = element_blank())

plt_k_measles_uncertainty <- df_CI_measles_uncertainty %>% 
  filter(param == 'k') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim, colour = as.factor(p_assumption))) +
  theme_bw() +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), position = position_dodge(0.5)) +
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
  scale_color_manual(values = c(col_measles_lower, col_measles, col_measles_upper),
                     name = 'Assumption for p') +
  theme(panel.grid.minor = element_blank())

plt_R_mers_uncertainty <- df_CI_mers_uncertainty %>% 
  filter(param == 'R') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim, colour = as.factor(p_assumption)))+
  theme_bw() +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), position = position_dodge(0.5)) +
  scale_y_continuous(limits = c(0., 1.),
                     breaks = seq(0.0, 1.6, 0.2),
                     name = 'R estimate',
                     expand = expansion(mult = c(0., 0.1))) +
  scale_x_discrete(name = '')  +
  scale_color_manual(values = c(col_mers_lower, col_mers, col_mers_upper),
                     name = 'Assumption for p') +
  theme(panel.grid.minor = element_blank())

plt_k_mers_uncertainty <- df_CI_mers_uncertainty %>% 
  filter(param == 'k') %>% 
  mutate(obs_process = case_when(p_detect == 1.0 ~  'All infections\ndetected',
                                 p_detect == 0.5 ~ 'Half of infections\ndetected',
                                 TRUE ~ as.character(p_detect))) %>% 
  ggplot(aes(x = as.factor(obs_process), y = mle_estim, colour = as.factor(p_assumption))) +
  theme_bw() +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95), position = position_dodge(0.5)) +
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
  scale_color_manual(values = c(col_mers_lower, col_mers, col_mers_upper),
                     name = 'Assumption for p') +
  theme(panel.grid.minor = element_blank())

panel_uncertainty_mers_measles <- ggarrange(
  ggarrange(plt_R_mers_uncertainty + ggtitle('MERS'), 
            plt_k_mers_uncertainty + ggtitle(''),
            common.legend = T, legend = 'bottom',
            labels = c('A', 'B')),
  ggarrange(plt_R_measles_uncertainty + ggtitle('Measles'), 
            plt_k_measles_uncertainty + ggtitle(''),
            common.legend = T, legend = 'bottom',
            labels = c('C', 'D')),
  nrow = 2
)

plot(panel_uncertainty_mers_measles)
