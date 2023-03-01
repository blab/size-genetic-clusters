library(tidyverse)
library(ggpubr)
library(RColorBrewer)

## Load inference results
df_inference <- readRDS('../results/sim_study_conditional_extinction/df_inference_cond_extinction.rds')

## Figure 2C
vec_k <- unique(df_inference$k)

plt_inference_k_cond_extinction <- df_inference %>% 
  filter(size_dataset == 1000,
         param == 'k',
         seed == 100) %>% 
  rename(true_val = k) %>% 
  ggplot(aes(x = as.factor(R0),
             group = as.factor(i_scenario),
             #group = interaction(true_val, R0, ratio_delay),
             shape = as.factor(ratio_delay),
             colour = as.factor(true_val))) +
  geom_hline(data = tibble(true_val = vec_k),
             aes(yintercept = true_val,
                 colour = as.factor(true_val)),
             linetype = 'dashed') +
  geom_point(aes(y = mle_estim),
             position = position_dodge(0.6)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.6)) +
  scale_x_discrete(name = 'True reproduction number R',
                   breaks = seq(0., 3.0, 0.25),
                   expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(trans = 'log',
                     name = 'Dispersion parameter k estimate',
                     expand = expansion(mult = c(0.0, 0.0)),
                     breaks = c(0.001, 0.002, 0.005,
                                0.01, 0.02, 0.05,
                                0.1, 0.2, 0.5,
                                1.0, 2.0, 5.0,
                                10.),
                     labels = c('0.001', '0.002', '0.005',
                                '0.01', '0.02', '0.05',
                                '0.1', '0.2', '0.5',
                                '1', '2', '5',
                                '10'),
                     limits = c(0.001, 10.)) +
  scale_shape_manual(name = expression(d[mut]/d[gen]),
                     values = 15:17,
                     breaks = c(0.6, 1, 5),
                     labels = c('0.6 (p = 38%)',
                                '1 (p = 50%)',
                                '5 (p = 83%)')) +
  scale_colour_manual(breaks = vec_k, 
                      values = brewer.pal(9, 'RdBu'),
                      name = 'True k') +
  theme_classic() +
  ggtitle('Inference conditional on cluster exctinction')

plot(plt_inference_k_cond_extinction)
