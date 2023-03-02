library(tidyverse)
library(ggpubr)

## Load inference results
df_inference <- readRDS('../results/ncov_NZ/df_inference_NZ_10000.rds') 

## Dispersion parameter estimates
plt_inference_k <- df_inference %>% 
  filter(param == 'k') %>% 
  ggplot(aes(x = period, y = mle_estim, 
             group = as.factor(p_detect),
             colour = as.factor(p_detect))) +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.4)) +
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
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks.x = element_blank())


## Reproduction number estimates
plt_inference_R <- df_inference %>% 
  filter(param != 'k') %>% 
  ggplot(aes(x = period, y = mle_estim, 
             group = as.factor(p_detect),
             colour = as.factor(p_detect))) +
  geom_hline(yintercept = 1.0, linetype = 'dashed', color = 'darkgrey') +
  geom_point(position = position_dodge(0.4)) +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95),
                 position = position_dodge(0.4)) +
  scale_x_discrete(name = '',
                   breaks = c('02) Apr-May 2020', '03) Jun-Dec 2020', 
                              '04) Jan-Apr 2021', '05) May-Jul 2021'),
                   labels = c('Apr-\nMay 20', 'Jun-\nDec 20', 'Jan-\nApr 21', 'May-\nJul 21')) +
  scale_y_continuous(limits = c(0.0, NA),
                     expand = expansion(mult = c(0.0, 0.05)),
                     name = 'R estimate') +
  scale_colour_manual(values = c('#002628', '#008b8b','#95b1b0'),
                      name = 'Proportion of infections\ndetected') +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

## Load and plot cluster size distribution
file_cluster_by_period <- '../data/ncov_NZ/df_cluster_by_period_NZ.rds'
df_name_period <- tibble(i_period = 1:4,
                         name_period = readRDS(file_cluster_by_period) %>% 
                           select(period) %>% unlist() %>% as.vector() %>% unique())

cluster_by_period <- readRDS(file_cluster_by_period) %>% 
  mutate(cluster_size_crop = case_when(cluster_size >= 5 ~ '5+',
                                       TRUE ~ as.character(cluster_size))) %>% 
  group_by(period, cluster_size_crop) %>% 
  summarise(n_clusters = sum(n_clusters)) %>% 
  group_by(period) %>% 
  mutate(prop = n_clusters/sum(n_clusters))

plt_period <- cluster_by_period %>%
  mutate(period = case_when(period == '02) Apr-May 2020' ~ 'Apr-May 20',
                            period == '03) Jun-Dec 2020' ~ 'Jun-Dec 20',
                            period == '04) Jan-Apr 2021' ~ 'Jan-Apr 21',
                            period == '05) May-Jul 2021' ~ 'May-Jul 21',
                            TRUE ~ 'Other'
  )) %>% 
  mutate(period = factor(period, levels = c('Apr-May 20', 'Jun-Dec 20', 'Jan-Apr 21', 'May-Jul 21'))) %>% 
  ggplot(aes(x = cluster_size_crop, y = n_clusters)) +
  geom_bar(stat = 'identity', position = position_dodge(0.8),
           fill = '#95b1b0') + 
  facet_wrap(. ~ period, nrow = 2, scales = 'free_y') +
  scale_x_discrete(name = 'Size of clusters\nof identical sequences') +
  scale_y_continuous(name = 'Number of clusters',
                     expand = expansion(mult = c(0., 0.05))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = '#002628'),
        strip.text = element_text(colour = 'white'))


## Panel with all the results
panel_NZ_zero_covid <- ggarrange(plt_period +
                                   ggtitle('SARS-CoV-2 in New Zealand'), 
                                 plt_inference_R + theme(legend.position = 'none'),
                                 plt_inference_k, ncol = 3,
                                 widths = c(1.25, 1.05, 1.))

plot(panel_NZ_zero_covid)


## Check the sensitivity of the estimates to max_cluster_size_inference (10,000 vs 50,000)
df_inference_2 <- readRDS('../results/ncov_NZ/df_inference_NZ_50000.rds')
df_inference %>% 
  left_join(df_inference_2, by = c('param', 'p_detect')) %>% 
  mutate(is_same = (lower_95.x == lower_95.y)) %>% 
  filter(! is_same) 
