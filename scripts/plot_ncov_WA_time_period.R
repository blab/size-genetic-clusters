library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Load inference results
df_inference <- readRDS('../results/ncov_WA/df_inference_time_period.rds') %>% 
  mutate(R_var_greater_R_non_var = (R_var_MLE_split > R_non_var_MLE_split))

# Plot the p-values for the different variants as a function of the time window used for the analysis
palette_for_variants <-  c('darkslateblue', '#FB9A99', '#A6CEE3', '#33A02C', '#FED976', '#FEB24C', '#FC4E2A')

plt_for_fig_5C <- df_inference %>% 
  filter(p_detect_cases == 0.5) %>%
  ggplot(aes(x = i_day, y = log10(p_val_split_vs_combined),  
             group = variant_char, colour = as.factor(variant_char))) +
  geom_point(aes(shape = as.factor(R_var_greater_R_non_var)), 
             alpha = 0.8) +
  scale_color_manual(name = 'Variant',
                     breaks = c(
                       'D614G', 'Epsilon', 'Alpha', 'Delta', 
                       'Omicron (BA.1)', 'Omicron (BA.2)', 'Omicron (BA.4/BA.5)'
                     ),
                     values = palette_for_variants) +
  scale_fill_manual(name = 'Variant',
                    # breaks = c('D614G', 'Epsilon', 'Alpha', 'Delta', 
                    #            'Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.45'),
                    breaks = c(
                      'D614G', 'Epsilon', 'Alpha', 'Delta', 
                      'Omicron (BA.1)', 'Omicron (BA.2)', 'Omicron (BA.4/BA.5)'
                    ),
                    values = c('darkslateblue', '#FB9A99', '#A6CEE3','#33A02C', 
                               '#FED976', '#FEB24C', '#FC4E2A')) +
  scale_x_continuous(name = 'Days after 10 variant sequences sampled',
                     expand = expansion(mult = c(0.05, 0.05)),
                     breaks = c(0, seq(10, 60, 10))) +
  scale_y_continuous(
    breaks = log10(c(seq(1e-6, 9e-6, 1e-6), seq(1e-5, 9e-5, 1e-5), seq(1e-4, 9e-4, 1e-4),
                     seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 1., 1e-1))),
    labels = c(expression("\u{2264} "*10^{-6}), '', '', '', '', '', '', '', '', 
               expression(10^{-5}), '', '', '', '', '', '', '', '', 
               expression(10^{-4}), '', '', '', '', '', '', '', '', 
               expression(10^{-3}), '', '', '', '', '', '', '', '',  
               expression(10^{-2}), '', '', '', '', '', '', '', '',  
               expression(10^{-1}), '', '', '', '', '', '', '', '',  
               expression(10^{0})),
    name = 'p-value',
    expand = expansion(mult = c(0.03, 0.03))) +
  scale_shape_manual(name = 'R MLE estimates consistent with\nvariant transmission advantage', 
                     breaks = c(T, F),
                     labels = c('Yes', 'No'),
                     values = c(16, 17)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.spacing.x = unit(0.8, "lines"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 12)) +
  facet_wrap(. ~ variant_char, ncol = 4)


# Sensitivity analysis exploring the impact of the proportion of infections detected as cases
plt_sensitivity_p_detect_cases <- df_inference %>% 
  ggplot(aes(x = i_day, y = log10(p_val_split_vs_combined),  
             group = variant_char, colour = as.factor(variant_char))) +
  geom_point(aes(shape = as.factor(R_var_greater_R_non_var)), 
             alpha = 0.8) +
  scale_color_manual(name = 'Variant',
                     breaks = c(
                       'D614G', 'Epsilon', 'Alpha', 'Delta', 
                       'Omicron (BA.1)', 'Omicron (BA.2)', 'Omicron (BA.4/BA.5)'
                     ),
                     values = palette_for_variants) +
  scale_fill_manual(name = 'Variant',
                    # breaks = c('D614G', 'Epsilon', 'Alpha', 'Delta', 
                    #            'Omicron_BA.1', 'Omicron_BA.2', 'Omicron_BA.45'),
                    breaks = c(
                      'D614G', 'Epsilon', 'Alpha', 'Delta', 
                      'Omicron (BA.1)', 'Omicron (BA.2)', 'Omicron (BA.4/BA.5)'
                    ),
                    values = c('darkslateblue', '#FB9A99', '#A6CEE3','#33A02C', 
                               '#FED976', '#FEB24C', '#FC4E2A')) +
  scale_x_continuous(name = 'Days after 10 variant sequences sampled',
                     expand = expansion(mult = c(0.05, 0.05)),
                     breaks = c(0, seq(10, 60, 10))) +
  scale_y_continuous(
    breaks = log10(c(seq(1e-6, 9e-6, 1e-6), seq(1e-5, 9e-5, 1e-5), seq(1e-4, 9e-4, 1e-4),
                     seq(1e-3, 9e-3, 1e-3), seq(1e-2, 9e-2, 1e-2), seq(1e-1, 1., 1e-1))),
    labels = c(expression("\u{2264} "*10^{-6}), '', '', '', '', '', '', '', '', 
               expression(10^{-5}), '', '', '', '', '', '', '', '', 
               expression(10^{-4}), '', '', '', '', '', '', '', '', 
               expression(10^{-3}), '', '', '', '', '', '', '', '',  
               expression(10^{-2}), '', '', '', '', '', '', '', '',  
               expression(10^{-1}), '', '', '', '', '', '', '', '',  
               expression(10^{0})),
    name = 'p-value',
    expand = expansion(mult = c(0.03, 0.03))) +
  scale_shape_manual(name = 'R MLE estimates consistent with\nvariant transmission advantage', 
                     breaks = c(T, F),
                     labels = c('Yes', 'No'),
                     values = c(16, 17)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.spacing.x = unit(0.8, "lines"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 12)) +
  facet_grid(variant_char ~ p_detect_cases)

plot(plt_sensitivity_p_detect_cases)
