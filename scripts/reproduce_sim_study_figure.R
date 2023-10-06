library(tidyverse)
library(ggpubr)
library(RColorBrewer)

## Load inference results
df_inference <- readRDS('../results/sim_study/df_inference_all.rds')

## Figure 2A
plt_R_inference <- df_inference %>%
  filter(k == 0.1, param == 'R',
         max_cluster_size_sim == 10000,
         max_cluster_size_inference == 10000,
         size_dataset == 1000,
         seed == 100) %>%
  rename(true_val = R0) %>%
  ggplot(aes(x = true_val, y = mle_estim,
             group = as.factor(i_scenario),
             shape = as.factor(p_detect),
             colour = as.factor(ratio_delay))) +
  geom_hline(yintercept = 1.6/0.6, colour = '#008b8b', linetype = 'dotted') +
  geom_hline(yintercept = 2.0, colour = '#ff8749', linetype = 'dotted') +
  geom_hline(yintercept = 6./5., colour = '#956f5d', linetype = 'dotted') +
  geom_segment(aes(x = true_val - 0.12, xend = true_val + 0.12,
                   y = true_val, yend = true_val),
               linetype = 'dashed',
               color = 'black')  +
  geom_point(position = position_dodge(0.15)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.15)) +
  scale_x_continuous(breaks = seq(0., 3.0, 0.5),
                     limits = c(0.35, 3.15),
                     name = 'True R') +
  scale_y_continuous(limits = c(0, 3.0),
                     expand = expansion(mult = c(0., 0.05)),
                     name = 'Estimated R',
                     breaks = seq(0., 3.0, 0.5)) +
  scale_shape_manual(name = 'Proportion of infections\nsequenced',
                     values = c(17, 16)) +
  scale_colour_manual(name = expression(d[mut]/d[gen]),
                      values = c('#008b8b', '#ff8749', '#956f5d'),
                      breaks = c(0.6, 1, 5),
                      labels = c('0.6 (p = 38%)',
                                 '1 (p = 50%)',
                                 '5 (p = 83%)')) +
  theme_classic() +
  theme(panel.grid.minor = element_blank())

## Figure 2B
plt_k_inference <- df_inference %>%
  filter(R0 == 1.0, param == 'k',
         max_cluster_size_sim == 10000,
         max_cluster_size_inference == 10000,
         size_dataset == 5000,
         seed == 100) %>%
  rename(true_val = k) %>%
  ggplot(aes(x = true_val, y = mle_estim,
             group = as.factor(i_scenario),
             shape = as.factor(p_detect),
             colour = as.factor(ratio_delay))) +
  geom_segment(aes(x = true_val*exp(0.35),
                   xend = true_val*exp(-0.35),
                   y = true_val, yend = true_val),
               linetype = 'dashed',
               color = 'black') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.4)) +
  scale_x_continuous(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0),
                     trans = 'log',
                     name = 'True k',
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
                     labels = scales::comma,
                     trans = 'log',
                     name = 'Estimated k',
                     expand = expansion(mult = c(0.05, 0.05))) +
  scale_shape_manual(name = 'Proportion of infections\nsequenced',
                     values = c(17, 16)) +
  scale_colour_manual(name = expression(d[mut]/d[gen]),
                      values = c('#008b8b', '#ff8749', '#956f5d'),
                      breaks = c(0.6, 1, 5),
                      labels = c('0.6 (p = 38%)',
                                 '1 (p = 50%)',
                                 '5 (p = 83%)')) +
  theme_classic() +
  theme(panel.grid.minor = element_blank())


plot(plt_R_inference)
plot(plt_k_inference)
