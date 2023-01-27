source('utils_distrib.R')

## This script performs simulations to quantify how much the largest transmission event
## contributes to the size of a cluster of identical sequences

### Reproduction number used for the simulation
R <- 1.0

### Probability that transmission occurs before mutation used in the analysis
p <- 0.7

### Dispersion parameters used for the simulations
vec_k <- c(0.05, 0.1, 1.0)

### Simulating cluster sizes
set.seed(1996)

df_sim <- Reduce('bind_rows', lapply(vec_k, FUN = function(k){
  simulate_branching_process_until_max_cluster_size_save_largest_gen(n_sim = 1e5,
                                                                     R = R,
                                                                     k = k, p,
                                                                     max_cluster_size = 1e5) %>% 
    mutate(k = k)
})) 


### Displaying the results of the simulation

## Proportion of the cluster due to the largest transmission event
plt_prop_contrib <- df_sim %>% 
  mutate(tot_cluster_size = tot_cluster_size - 1) %>% 
  filter(tot_cluster_size >= 1) %>% 
  mutate(cat_tot_cluster_size = factor(case_when(tot_cluster_size < 20 ~ '<20',
                                                 tot_cluster_size < 40 ~ '[20-40)',
                                                 tot_cluster_size < 60 ~ '[40-60)',
                                                 tot_cluster_size < 80 ~ '[60-80)',
                                                 tot_cluster_size < 100 ~ '[80-100)',
                                                 tot_cluster_size < 120 ~ '[100-120)',
                                                 tot_cluster_size < 140 ~ '[120-140)',
                                                 TRUE ~ '140+'),
                                       levels = c('<20', '[20-40)', '[40-60)', '[60-80)',
                                                  '[80-100)', '[100-120)', '[120-140)', 
                                                  '140+'))) %>% 
  group_by(k, cat_tot_cluster_size) %>% 
  mutate(quantity_of_interest = contrib_largest_transmission_event/tot_cluster_size) %>% 
  summarise(y0 = min(quantity_of_interest),
            y05 = quantile(quantity_of_interest, 0.95),
            y25 = quantile(quantity_of_interest, 0.25),
            y50 = median(quantity_of_interest),
            y75 = quantile(quantity_of_interest, 0.75),
            y95 = quantile(quantity_of_interest, 0.05),
            y100 = max(quantity_of_interest)) %>% 
  ungroup() %>% 
  ggplot(aes(cat_tot_cluster_size)) +
  geom_boxplot(aes(ymin = y05, lower = y25, middle = y50, upper = y75, ymax = y95,
                   fill = as.factor(k)),
               stat = "identity",
               width = 0.6,
               position = position_dodge(0.8)) +
  scale_fill_manual(name = 'Dispersion\nparameter k',
                    breaks = c(0.05, 0.1, 1),
                    values = c('#56423d', '#bea69f', '#ffe6d8')) +
  theme_bw() +
  scale_x_discrete(name = '\nTotal number of offsprings\nwith identical sequences') +
  scale_y_continuous(name = 'Contribution of the largest\ntransmission event to cluster size',
                     limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.05))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0))



## Size of the largest transmission event as a function of the size of the cluster
plt_size_largest_cluster <- df_sim %>% 
  mutate(tot_cluster_size = tot_cluster_size - 1) %>% 
  filter(tot_cluster_size >= 1) %>% 
  mutate(cat_tot_cluster_size = factor(case_when(tot_cluster_size < 20 ~ '<20',
                                                 tot_cluster_size < 40 ~ '[20-40)',
                                                 tot_cluster_size < 60 ~ '[40-60)',
                                                 tot_cluster_size < 80 ~ '[60-80)',
                                                 tot_cluster_size < 100 ~ '[80-100)',
                                                 tot_cluster_size < 120 ~ '[100-120)',
                                                 tot_cluster_size < 140 ~ '[120-140)',
                                                 TRUE ~ '140+'),
                                       levels = c('<20', '[20-40)', '[40-60)', '[60-80)',
                                                  '[80-100)', '[100-120)', '[120-140)', 
                                                  '140+'))) %>% 
  group_by(k, cat_tot_cluster_size) %>% 
  mutate(quantity_of_interest = contrib_largest_transmission_event) %>% 
  summarise(y0 = min(quantity_of_interest),
            y05 = quantile(quantity_of_interest, 0.95),
            y25 = quantile(quantity_of_interest, 0.25),
            y50 = median(quantity_of_interest),
            y75 = quantile(quantity_of_interest, 0.75),
            y95 = quantile(quantity_of_interest, 0.05),
            y100 = max(quantity_of_interest)) %>% 
  ungroup() %>% 
  ggplot(aes(cat_tot_cluster_size)) +
  geom_boxplot(aes(ymin = y05, lower = y25, middle = y50, upper = y75, ymax = y95,
                   fill = as.factor(k)),
               stat = "identity",
               width = 0.6,
               position = position_dodge(0.8)) +
  scale_fill_manual(name = 'Dispersion\nparameter k',
                    breaks = c(0.05, 0.1, 1),
                    values = c('#56423d', '#bea69f', '#ffe6d8')) +
  theme_bw() +
  scale_x_discrete(name = '\nTotal number of offsprings\nwith identical sequences') +
  scale_y_continuous(name = 'Number of offsprings with\nidentical sequences generated\nduring largest transmission event',
                     trans = 'log',
                     breaks = c(1, 2, 5, 10, 20, 50, 
                                100, 200, 500,  1000, 2000, 5000,
                                10000)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        panel.grid.minor = element_blank())


panel_1.0 <- ggpubr::ggarrange(plt_size_largest_cluster, plt_prop_contrib,
                               ncol = 2,
                               common.legend = T)

plot(panel_1.0)
