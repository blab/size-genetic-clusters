library(tidyverse)
library(ggpubr)

## Load simulation results
df_SNPs <- readRDS('../results/withinhost_sim/df_SNPs.rds')

## Generate proportion of pairs with same consensus genome
df_SNPs_summary <- df_SNPs %>% 
  group_by(id_pair, i_rep, i_scenario, rep) %>% 
  mutate(is_different_consensus = ((freq_donor >= 0.5 & freq_recepient < 0.5) | (freq_donor < 0.5 & freq_recepient > 0.5))) %>% 
  summarise(n_different_loci_consensus = sum(is_different_consensus)) %>% 
  group_by(i_scenario) %>% 
  summarise(n_different_consensus = sum(n_different_loci_consensus > 0.),
            n_pairs = n(), 
            prop_different_consensus = n_different_consensus/n_pairs,
            n_same_consensus = n_pairs - n_different_consensus,
            p_same_consensus = n_same_consensus/n_pairs) %>% 
  left_join(df_scenarios, by = 'i_scenario') %>% 
  mutate(ratio_d_mut_d_gen = av_delay_between_mutations/inf_duration,
         p_trans_before_mut = ratio_d_mut_d_gen/(ratio_d_mut_d_gen + 1))


df_SNPs_summary <- bind_cols(df_SNPs_summary, 
                             Reduce('bind_rows', lapply(1:nrow(df_SNPs_summary), FUN = function(i_row){
                               vec_CI <- prop.test(x = df_SNPs_summary$n_same_consensus[i_row],
                                                   n = df_SNPs_summary$n_pairs[i_row])$conf.int %>%
                                 as.numeric()
                               names(vec_CI) <- c('lower_prop', 'upper_prop')
                               vec_CI
                             })))

## Plot the proportion of pairs with identical consensus genomes across scenarios and compare it with the probability
## that a transmission event occurs before a mutation one
plt_estim_p_trans_before_mut <- df_SNPs_summary %>% 
  mutate(bottleneck_size_char = paste0('Bottleneck size = ', bottleneck_size),
         bottleneck_size_char = factor(bottleneck_size_char, levels = paste0('Bottleneck size = ', c(1, 2, 10))),
         inf_duration_char = paste0('Inf. duration = ', inf_duration),
         inf_duration_char = factor(inf_duration_char, levels = paste0('Inf. duration = ', c(6, 12, 100, 200)))) %>% 
  ggplot(aes(x = as.factor(bottleneck_size), y = p_same_consensus)) +
  geom_hline(aes(yintercept = p_trans_before_mut),
             linetype = 'dashed', color = 'darkgrey') +
  geom_point(size = 2) +
  geom_linerange(aes(ymin = lower_prop, ymax = upper_prop)) +
  theme_classic() +
  facet_wrap(~ inf_duration_char, ncol = 4, 
             scales = 'free_y') +
  scale_x_discrete(name = 'Transmission bottleneck size') +
  scale_y_continuous(limits = c(0., NA),
                     expand = expansion(mult = c(0., 0.05)),
                     name = 'Proportion of pairs with\nidentical consensus') +
  theme(strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

plot(plt_estim_p_trans_before_mut)


## Frequency of iSNV in donor and recipient (Sanity check that we are actually simulating the right transmission bottleneck scenario)
plt_freq_donor_recipient <- df_SNPs %>% 
  filter(rep == 1) %>% # Don't plot all the pairs (plot is a bit big otherwise)
  filter(freq_donor + freq_recepient > 0.) %>% 
  group_by(id_pair, i_rep, i_scenario) %>% 
  mutate(fixed_in_donor = freq_donor > 0.5,
         fixed_in_recipient = freq_recepient > 0.5,
         same_fixation_donor_recipient = (fixed_in_donor == fixed_in_recipient)) %>% 
  mutate(bottleneck_size_char = paste0('Bottleneck size = ', bottleneck_size),
         bottleneck_size_char = factor(bottleneck_size_char, levels = paste0('Bottleneck size = ', c(1, 2, 10))),
         inf_duration_char = paste0('Inf. duration = ', inf_duration),
         inf_duration_char = factor(inf_duration_char, levels = paste0('Inf. duration = ', c(6, 12, 100, 200)))) %>% 
  ggplot(aes(x = freq_donor, y = freq_recepient)) +
  geom_point(alpha = 0.005, colour = 'gray22') +
  facet_grid(inf_duration_char ~ bottleneck_size_char) +
  theme_classic() +
  scale_x_continuous(name = 'Frequency in donor', limits = c(0., 1.)) +
  scale_y_continuous(name = 'Frequencing in recipient', limits = c(0., 1.)) +
  coord_fixed() +
  theme(strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white', size = 12),
        axis.text = element_text(size = 12),
        panel.spacing = unit(1, 'lines'))

plot(plt_freq_donor_recipient)
