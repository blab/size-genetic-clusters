library(tidyverse)
library(RColorBrewer)

## Load inference results
df_inference <- readRDS('../results/sim_study_trans_adv/df_res_sim_study_trans_adv.rds')

vec_increase_R <- sort(unique(df_inference$increase_R))


## Functions to plot the results
get_plot_sensitivity <- function(curr_baseline_R = 0.75, curr_ratio_delay = 1.0, curr_k = 0.1){
  
  curr_p <- curr_ratio_delay/ (1. + curr_ratio_delay)
  curr_trans_adv_threshold <- 1. / (curr_baseline_R * curr_p) - 1.
  x_line_to_draw <- which(vec_increase_R[-1] > curr_trans_adv_threshold)[1] - 0.5
  
  plt_sensitivity_ratio <- df_inference %>% 
    filter(k == curr_k, 
           increase_R > 0,
           baseline_R == curr_baseline_R,
           abs(ratio_delay_mutation_transmission - curr_ratio_delay) < 1e-6) %>% 
    mutate(is_significant = (p_val < 0.05)) %>% 
    group_by(baseline_R, increase_R, k, ratio_delay_mutation_transmission, size_dataset) %>%
    summarise(n_signif = sum(is_significant),
              n_tot = n()) %>% ungroup() %>%
    mutate(sensitivity = n_signif/n_tot,
           rank_increase_R = dense_rank(increase_R)) %>%
    ggplot(aes(x = rank_increase_R,
               y = sensitivity,
               group = interaction(size_dataset, ratio_delay_mutation_transmission, rank_increase_R),
               fill = as.factor(size_dataset))) +
    geom_bar(position = position_dodge(0.8), width = 0.8,
             stat = 'identity') +
    scale_x_continuous(name = 'True transmission advantage',
                       breaks = 1:(length(vec_increase_R) - 1),
                       label = paste0(vec_increase_R[-1]*100, '%')) +
    scale_y_continuous(name = 'Sensitivity', limits = c(0., 1.),
                       expand = expansion(mult = c(0., 0.0))) +
    scale_fill_manual(breaks = c(50, 100, 1000, 5000),
                      name = 'Number of clusters',
                      values = rev(brewer.pal(8, 'Spectral')[1:4])) +
    theme_classic() +
    theme(plot.margin = unit(c(2, 1, 1, 1), "lines")) +
    coord_cartesian(clip = "off") +
    annotation_custom(grob = linesGrob(gp = gpar(lty = 'dashed')),
                      xmin = x_line_to_draw, xmax = x_line_to_draw,
                      ymin = 0., ymax = 1.1) +
    annotation_custom(grob = textGrob(expression(R[V] < 1/p),
                                      hjust = 1.0), 
                      xmin = x_line_to_draw - 0.2, xmax = x_line_to_draw - 0.2,
                      ymin = 1.05, ymax = 1.05) +
    annotation_custom(grob = textGrob(expression(R[V] > 1/p),
                                      hjust = 0.0), 
                      xmin = x_line_to_draw + 0.2, xmax = x_line_to_draw + 0.2,
                      ymin = 1.05, ymax = 1.05)
  
  return(plt_sensitivity_ratio)
}
get_plot_bias_trans_adv <- function(curr_baseline_R = 0.75, curr_ratio_delay = 1.0, curr_k = 0.1){
  
  curr_p <- curr_ratio_delay/ (1. + curr_ratio_delay)
  curr_trans_adv_threshold <- 1. / (curr_baseline_R * curr_p) - 1.
  x_line_to_draw <- which(vec_increase_R[-1] > curr_trans_adv_threshold)[1] - 0.5
  
  
  plt_bias_trans_adv <- df_inference %>% 
    filter(baseline_R == curr_baseline_R, k == curr_k, increase_R > 0, abs(ratio_delay_mutation_transmission - curr_ratio_delay) < 1e-6) %>% 
    mutate(is_significant = (p_val < 0.05)) %>% 
    mutate(transm_adv = R_var_MLE_split / R_non_var_MLE_split) %>%
    mutate(quantity_of_interest = (transm_adv - 1.0 - increase_R)) %>% 
    group_by(increase_R, size_dataset, ratio_delay_mutation_transmission) %>% 
    summarise(y0 = min(quantity_of_interest),
              y05 = quantile(quantity_of_interest, 0.95),
              y25 = quantile(quantity_of_interest, 0.25),
              y50 = median(quantity_of_interest),
              y75 = quantile(quantity_of_interest, 0.75),
              y95 = quantile(quantity_of_interest, 0.05),
              y100 = max(quantity_of_interest)) %>% 
    ggplot(aes(x = as.factor(increase_R),
               fill = as.factor(size_dataset))) +
    geom_hline(yintercept = 0., linetype = 'dashed', color = 'darkgrey') +
    geom_boxplot(aes(ymin = y05, lower = y25, middle = y50, upper = y75, ymax = y95),
                 position = position_dodge(0.7), width = 0.6,
                 stat = "identity") +
    scale_x_discrete(name = 'True transmission advantage',
                     breaks = vec_increase_R,
                     label = paste0(vec_increase_R*100, '%')) +
    scale_y_continuous(name = 'Transmission advantage bias') +
    scale_fill_manual(breaks = c(50, 100, 1000, 5000),
                      name = 'Number of clusters',
                      values = rev(brewer.pal(8, 'Spectral')[1:4])) +
    theme_classic() +
    theme(plot.margin = unit(c(2, 1, 1, 1), "lines")) +
    coord_cartesian(clip = "off") +
    annotation_custom(grob = linesGrob(gp = gpar(lty = 'dashed')),
                      xmin = x_line_to_draw, xmax = x_line_to_draw,
                      ymin = -3.2, ymax = 3.1) +
    annotation_custom(grob = textGrob(expression(R[V] < 1/p),
                                      hjust = 1.0), 
                      xmin = x_line_to_draw - 0.2, xmax = x_line_to_draw - 0.3,
                      ymin = 3., ymax = 3.) +
    annotation_custom(grob = textGrob(expression(R[V] > 1/p),
                                      hjust = 0.0), 
                      xmin = x_line_to_draw + 0.2, xmax = x_line_to_draw + 0.2,
                      ymin = 3.0, ymax = 3.0)
  
  return(plt_bias_trans_adv)
  
}
get_boxplot_bias_k <- function(curr_baseline_R = 0.75, curr_ratio_delay = 1.0, curr_k = 0.1,
                               curr_size_dataset = 5000){
  plt <- df_inference %>% 
    filter(abs(ratio_delay_mutation_transmission - curr_ratio_delay) < 1e-6,
           baseline_R == curr_baseline_R,
           increase_R > 0.,
           k == curr_k,
           size_dataset == curr_size_dataset) %>% 
    group_by(size_dataset, ratio_delay_mutation_transmission, increase_R) %>% 
    mutate(quantity_of_interest_1 = k_MLE_combined,
           quantity_of_interest_2 = k_MLE_split) %>% 
    summarise(y0_1 = min(quantity_of_interest_1),
              y05_1 = quantile(quantity_of_interest_1, 0.975),
              y25_1 = quantile(quantity_of_interest_1, 0.25),
              y50_1 = median(quantity_of_interest_1),
              y75_1 = quantile(quantity_of_interest_1, 0.75),
              y95_1 = quantile(quantity_of_interest_1, 0.025),
              y100_1 = max(quantity_of_interest_1),
              y0_2 = min(quantity_of_interest_2),
              y05_2 = quantile(quantity_of_interest_2, 0.975),
              y25_2 = quantile(quantity_of_interest_2, 0.25),
              y50_2 = median(quantity_of_interest_2),
              y75_2 = quantile(quantity_of_interest_2, 0.75),
              y95_2 = quantile(quantity_of_interest_2, 0.025),
              y100_2 = max(quantity_of_interest_2)) %>% 
    pivot_longer(cols = -c('size_dataset', 'ratio_delay_mutation_transmission', 'increase_R'),
                 values_to = 'value',
                 names_sep = '_',
                 names_to = c('quantile', 'id')) %>% 
    ungroup() %>% 
    pivot_wider(names_from = 'quantile', values_from = 'value') %>% 
    mutate(rank_increase_R = dense_rank(increase_R)) %>% 
    ggplot(aes(x = as.numeric(rank_increase_R))) +
    geom_hline(yintercept = 0.1, linetype = 'dashed', color = 'darkgrey') +
    geom_boxplot(aes(ymin = y05, lower = y25, middle = y50, upper = y75, ymax = y95,
                     group = interaction(increase_R, size_dataset, id),
                     fill = as.factor(id)),
                 stat = "identity",
                 position = position_dodge(0.6),
                 width = 0.5)  +
    scale_y_continuous(trans = 'log', name = 'k estimate',
                       expand = expansion(mult = c(0.05, 0.05)),
                       breaks = c(0.01, 0.02, 0.03, 0.04, 0.05,
                                  0.06, 0.07, 0.08, 0.09, 
                                  0.1, 0.11, 0.12, 0.13,
                                  0.2, 0.3, 0.4, 0.5, 1.0)) +
    scale_x_continuous(name = 'True transmission advantage',
                       breaks = 1:(length(vec_increase_R) - 1),
                       labels = paste0(vec_increase_R[-1]*100, '%')) +
    scale_fill_manual(name = 'Accounting for\ntwo genetic\nsubpopulations',
                      values = c('gray60', 'darkcyan'),
                      breaks = 1:2,
                      labels = c('No', 'Yes')) +
    theme_classic() +
    theme(plot.margin = unit(c(2, 1, 1, 1), "lines")) +
    coord_cartesian(ylim = c(0.06, 0.12), clip = 'off')
  
  return(plt)
}

## Sensitivity of our inference framework
get_plot_sensitivity()

## Bias in the transmission advantage estimation
get_plot_bias_trans_adv()

## Bias in the dispersion parameter when not accounting for 
## two genetic subpopulations
get_boxplot_bias_k()

