library(tidyverse)
library(ggpubr)

## Load inference results
df_inference <- readRDS('../results/ncov_NZ/df_inference_10000_central.rds')
df_inference_lower <- readRDS('../results/ncov_NZ/df_inference_10000_lower.rds')
df_inference_upper <- readRDS('../results/ncov_NZ/df_inference_10000_upper.rds')

df_uncertainty <- bind_rows(
  df_inference %>% mutate(p_assumption = 'Central'),
  df_inference_lower %>% mutate(p_assumption = 'Lower bound'),
  df_inference_upper %>% mutate(p_assumption = 'Upper bound')
)%>% 
  mutate(p_assumption = factor(p_assumption, levels = c('Lower bound', 'Central', 'Upper bound')))

## Plot results
plt_k_uncertainty <- df_uncertainty %>% 
  filter(param == 'k', p_detect >= 0.5) %>% 
  ggplot(aes(x = period, y = mle_estim, 
             group = interaction(p_assumption, p_detect),
             shape = as.factor(p_assumption),
             colour = as.factor(p_detect))) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.5)) +
  scale_x_discrete(name = '', breaks = 'All', labels = '') +
  scale_y_continuous(trans = 'log',
                     limits = c(0.1, 2.0),
                     breaks = c(seq(0.1, 2.0, 0.1)),
                     labels = c(seq(0.1, 0.5, 0.1), rep('', 4), 
                                1.0, rep('', 4), 1.5, rep('', 4), 2.0),
                     expand = expansion(mult = c(0.05, 0.01)),
                     name = 'k estimate') +
  scale_colour_manual(values = c('#002628', '#008b8b','#95b1b0'),
                      name = 'Proportion of\ninfections\ndetected') +
  scale_shape_manual(name = 'Assumption for p', values = c(15, 16, 17)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks.x = element_blank())


plt_R_uncertainty <- df_uncertainty %>% 
  filter(param != 'k', p_detect >= 0.5) %>% 
  ggplot(aes(x = period, y = mle_estim, 
             group = interaction(p_assumption, p_detect),
             shape = as.factor(p_assumption),
             colour = as.factor(p_detect))) +
  geom_hline(yintercept = 1.0, linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.5)) +
  scale_x_discrete(name = '',
                   breaks = c('02) Apr-May 2020', '03) Jun-Dec 2020', 
                              '04) Jan-Apr 2021', '05) May-Jul 2021'),
                   labels = c('Apr-\nMay 20', 'Jun-\nDec 20', 'Jan-\nApr 21', 'May-\nJul 21')) +
  scale_y_continuous(limits = c(0.0, NA),
                     breaks = seq(0., 1.25, 0.25),
                     expand = expansion(mult = c(0.0, 0.05)),
                     name = 'R estimate') +
  scale_shape_manual(name = 'Assumption for p', values = c(15, 16, 17)) +
  scale_colour_manual(values = c('#002628', '#008b8b','#95b1b0'),
                      name = 'Proportion of infections\ndetected') +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plt_uncertainty <- ggarrange(plt_R_uncertainty, plt_k_uncertainty, common.legend = T,
                             widths = c(1., 0.5), legend = 'right', labels = 'AUTO')

plot(plt_uncertainty)
