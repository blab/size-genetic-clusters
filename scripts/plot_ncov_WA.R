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
  geom_text(show_guide  = FALSE, size = 3) +
  scale_x_discrete(name = '\nProportion of infections detected') +
  scale_color_manual(name = 'Signficance level',
                     values = c('darkred', 'orange2', 'darkgrey'),
                     breaks = c('< 1%', '< 5%', 'Non significant')) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(color = ))
