library(tidyverse)

# Load inference results
df_inference <- readRDS('../results/ncov_WA/df_inference.rds')

# Display p-values
df_inference %>% 
  mutate(date_beginning_window = as.Date(date_beginning_window)) %>% 
  select(variant, date_beginning_window, p_val_split_vs_combined, p_detect_cases) %>% 
  arrange(date_beginning_window) %>% 
  mutate(p_val_split_vs_combined_char = case_when(p_val_split_vs_combined == 0. ~ '<1e-5',
                                                  TRUE ~ formatC(p_val_split_vs_combined, format = "e", digits = 1)
  )) %>% 
  mutate(signif_level = case_when(p_val_split_vs_combined < 0.01 ~ '< 1%',
                                  p_val_split_vs_combined < 0.05 ~ '< 5%',
                                  TRUE ~ 'Non significant')) %>% 
  mutate(variant = case_when(variant == 'Omicron_BA.1' ~ 'Omicron (BA.1)',
                             variant == 'Omicron_BA.2' ~ 'Omicron (BA.2)',
                             variant == 'Omicron_BA.45' ~ 'Omicron (BA.4, BA.5)',
                             TRUE ~ variant),
         variant = factor(variant, levels = c('D614G', '20G', 'Epsilon', 'Alpha', 'Delta', 
                                              'Omicron (BA.1)',
                                              'Omicron (BA.2)', 'Omicron (BA.4, BA.5)'))
  ) %>% 
  ggplot(aes(y = as.factor(variant),
             x = as.factor(paste0(p_detect_cases*100, '%')),
             label = p_val_split_vs_combined_char,
             colour = signif_level)) +
  geom_point(alpha = 0.0) +
  geom_text(size = 3) +
  scale_x_discrete(name = '\nProportion of infections detected') +
  scale_color_manual(name = 'Signficance level',
                     values = c('darkred', 'orange2', 'darkgrey'),
                     breaks = c('< 1%', '< 5%', 'Non significant')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())


# Plot cluster size distribution
seq_name_variants <- c('D614G', 'Epsilon', 'Alpha', 'Delta',
                       'Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.45')
seq_label_variants <- c('non D614G', 'Epsilon', 'Alpha', 'Delta',
                        'Omicron (BA.1)', 'Omicron (BA.2)', 'Omicron (BA.4, BA.5)')

cluster_size_threshold <- 9
df_cluster_distrib <- Reduce('bind_rows', lapply(1:length(seq_name_variants), FUN = function(i_var){
  name_variant <- seq_name_variants[i_var]
  label_variant <- seq_label_variants[i_var]
  readRDS(paste0('../data/ncov_WA/cluster_distrib_', name_variant, '.rds'))$df_cluster_distrib %>% 
    mutate(cluster_size = case_when(n_in_same_clique >= cluster_size_threshold ~ paste0(as.character(cluster_size_threshold), '+'),
                                    TRUE ~ as.character(n_in_same_clique))) %>% 
    group_by(is_variant, cluster_size) %>% 
    summarise(n_clusters = sum(n_clusters)) %>% 
    mutate(variant = name_variant, label_var = label_variant)
}))

plt_cluster_distrib <- df_cluster_distrib %>% ungroup() %>% 
  mutate(label_var = factor(label_var, levels = seq_label_variants)) %>% 
  ggplot() +
  geom_text(x = Inf, y = Inf,
            aes(label = label_var,
                color = (cluster_size == paste0(cluster_size_threshold, '+'))),
            vjust = 1.0, hjust = 1.0) +
  geom_bar(aes(x = cluster_size, y = n_clusters, group = is_variant, fill = is_variant),
           stat = 'identity', position = position_dodge(0.8)) +
  facet_wrap(label_var ~ is_variant, scales = 'free_y', ncol = 2) +
  scale_x_discrete(name = 'Cluster size') +
  scale_y_continuous(name = 'Number of clusters') +
  scale_colour_manual(values = c('white', 'black'), guide = 'none') +
  scale_fill_manual(name = 'Is variant?', values = c('darkgrey', 'darkred'),
                    breaks = c(F, T), labels = c('No', 'Yes')) +
  theme_classic() +
  theme(strip.text.x = element_text(size = 0),
        strip.background = element_blank())
plot(plt_cluster_distrib)
