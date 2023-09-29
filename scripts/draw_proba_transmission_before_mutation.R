library(dplyr)
library(readr)
library(ggpubr)
library(ggrepel)
library(ggh4x)
library(RColorBrewer)
source('utils_proba_trans_before_mut.R')

## Load values for the mutation rate and the generation time
df_pathogens <- read_csv('../data/proba_trans_before_mut/gt_mutation_rate_with_uncertainty.csv')

## Convert mutation rate in mutation per day
df_pathogens$subs_per_day <- sapply(1:nrow(df_pathogens), FUN = function(i_row){
  if(! is.na(df_pathogens$subs_per_site_per_day[i_row])){
    df_pathogens$subs_per_site_per_day[i_row] * df_pathogens$genome_length[i_row]
  } else if(! is.na(df_pathogens$subs_per_year[i_row])){
    df_pathogens$subs_per_year[i_row] / 365.25
  } else if(! is.na(df_pathogens$subs_per_site_per_year[i_row])) {
    df_pathogens$subs_per_site_per_year[i_row] * df_pathogens$genome_length[i_row] / 365.25
  } 
})
df_pathogens$subs_per_day_lower <- sapply(1:nrow(df_pathogens), FUN = function(i_row){
  df_pathogens$subs_per_site_per_year_lower[i_row] * df_pathogens$genome_length[i_row] / 365.25
})
df_pathogens$subs_per_day_upper <- sapply(1:nrow(df_pathogens), FUN = function(i_row){
  df_pathogens$subs_per_site_per_year_upper[i_row] * df_pathogens$genome_length[i_row] / 365.25
})

## Run the simulations for all pathogens
n_sim <- 1e7
set.seed(9876)
list_sim_pathogens <-  lapply(1:nrow(df_pathogens), FUN = function(i_pathogen){
  sim_pathogen <- compute_proba_transmission_before_mutation_rate_per_day(n_sim = n_sim,
                                                                          mean_GT = df_pathogens$mean_GT[i_pathogen],
                                                                          sd_GT = df_pathogens$sd_GT[i_pathogen],
                                                                          rate_sub_per_day = df_pathogens$subs_per_day[i_pathogen])
})
list_sim_pathogens_lower_GT <-  lapply(1:nrow(df_pathogens), FUN = function(i_pathogen){
  sim_pathogen <- compute_proba_transmission_before_mutation_rate_per_day(n_sim = n_sim,
                                                                          mean_GT = df_pathogens$lower_mean_GT[i_pathogen],
                                                                          sd_GT = df_pathogens$sd_GT[i_pathogen],
                                                                          rate_sub_per_day = df_pathogens$subs_per_day[i_pathogen])
})
list_sim_pathogens_upper_GT <-  lapply(1:nrow(df_pathogens), FUN = function(i_pathogen){
  sim_pathogen <- compute_proba_transmission_before_mutation_rate_per_day(n_sim = n_sim,
                                                                          mean_GT = df_pathogens$upper_mean_GT[i_pathogen],
                                                                          sd_GT = df_pathogens$sd_GT[i_pathogen],
                                                                          rate_sub_per_day = df_pathogens$subs_per_day[i_pathogen])
})
list_sim_pathogens_lower_mut_rate <-  lapply(1:nrow(df_pathogens), FUN = function(i_pathogen){
  sim_pathogen <- compute_proba_transmission_before_mutation_rate_per_day(n_sim = n_sim,
                                                                          mean_GT = df_pathogens$mean_GT[i_pathogen],
                                                                          sd_GT = df_pathogens$sd_GT[i_pathogen],
                                                                          rate_sub_per_day = df_pathogens$subs_per_day_lower[i_pathogen])
})
list_sim_pathogens_upper_mut_rate <-  lapply(1:nrow(df_pathogens), FUN = function(i_pathogen){
  sim_pathogen <- compute_proba_transmission_before_mutation_rate_per_day(n_sim = n_sim,
                                                                          mean_GT = df_pathogens$mean_GT[i_pathogen],
                                                                          sd_GT = df_pathogens$sd_GT[i_pathogen],
                                                                          rate_sub_per_day = df_pathogens$subs_per_day_upper[i_pathogen])
})


## Reproduction number threshold for different pathogens
df_p_trans_before_mut <- Reduce('bind_rows', lapply(1:length(list_sim_pathogens), FUN = function(i_pathogen){
  c('pathogen' = df_pathogens$pathogen[i_pathogen],
    'p_trans_before_mut' = list_sim_pathogens[[i_pathogen]]$p_trans_before_mut,
    'p_trans_before_mut_lower_GT' = list_sim_pathogens_lower_GT[[i_pathogen]]$p_trans_before_mut,
    'p_trans_before_mut_upper_GT' = list_sim_pathogens_upper_GT[[i_pathogen]]$p_trans_before_mut,
    'p_trans_before_mut_lower_mut_rate' = list_sim_pathogens_lower_mut_rate[[i_pathogen]]$p_trans_before_mut,
    'p_trans_before_mut_upper_mut_rate' = list_sim_pathogens_upper_mut_rate[[i_pathogen]]$p_trans_before_mut
  )
})) %>% 
  mutate(p_trans_before_mut = as.numeric(p_trans_before_mut),
         p_trans_before_mut_lower_GT = as.numeric(p_trans_before_mut_lower_GT),
         p_trans_before_mut_upper_GT = as.numeric(p_trans_before_mut_upper_GT),
         p_trans_before_mut_lower_mut_rate = as.numeric(p_trans_before_mut_lower_mut_rate),
         p_trans_before_mut_upper_mut_rate = as.numeric(p_trans_before_mut_upper_mut_rate),
         R_max = 1./p_trans_before_mut) %>% 
  group_by(pathogen) %>% 
  mutate(lower_p_trans_before_mut = min(p_trans_before_mut_upper_GT, p_trans_before_mut_upper_mut_rate, na.rm = T),
         upper_p_trans_before_mut = max(p_trans_before_mut_lower_GT, p_trans_before_mut_lower_mut_rate, na.rm = T))

#saveRDS(df_p_trans_before_mut, '../results/df_p_trans_before_mut_with_uncertainty.rds')

## Display the different estimates
breaks_p_trans_before_mut <- seq(0., 1., 0.2)
breaks_R_to_plot <- c(1., 1.25, 1.5, 2.0, 2.5, 5.0)
corr_break_R <- 1./breaks_R_to_plot
                  
plt_p_uncertainty <- df_p_trans_before_mut %>% 
  left_join(df_pathogens, by = 'pathogen') %>% 
  mutate(expected_nb_mut_gen = mean_GT*subs_per_day) %>% 
  ggplot(aes(x = expected_nb_mut_gen,
             y = p_trans_before_mut,
             colour = as.factor(pathogen))) +
  geom_point(size = 1.0) +
  geom_linerange(aes(ymin = lower_p_trans_before_mut, ymax = upper_p_trans_before_mut)) +
  geom_text_repel(aes(label = pathogen)) +
  scale_x_continuous(name = 'Expected number of mutations during 1 generation',
                     #limits = c(0.2, NA),
                     breaks = c(0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 2.0)) +
  scale_y_continuous(limits = c(0., 1.), breaks = breaks_p_trans_before_mut,
                     expand = expansion(mult = c(0., 0.)),
                     name = 'Probability p that transmission\noccurs before mutation',
                     sec.axis = sec_axis(~ . * 1., name = "Reproduction number threshold", 
                                         breaks = corr_break_R,
                                         labels = breaks_R_to_plot)
  ) +
  geom_blank(data = tibble(expected_nb_mut_gen = c(0.1, 0.8, 1.0, 2.0),
                           y = c(0.5, 0.5, 0.5, 0.5),
                           pathogen = c(NA, NA, NA, NA),
                           p_trans_before_mut = c(0.5, 0.5, 0.5, 0.5))) +
  scale_colour_manual(values = c(brewer.pal(12, 'Paired')[-11], 'black')) +
  theme_classic() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_blank()) +
  facet_wrap(. ~ (expected_nb_mut_gen >= 1.0), scales = 'free_x') +
  force_panelsizes(cols = c(1.0, 0.15))


plot(plt_p_uncertainty)


### Display results
df_p_trans_before_mut %>% 
  left_join(df_pathogens, by = 'pathogen') %>% 
  mutate(pathogen = factor(pathogen, levels = c('MERS', 'Measles', 'Ebola', 'Zika', 'Mpox', 
                                                'Influenza A (H1N1pdm)', "Influenza A (H3N2)",
                                                'Mumps', 'RSV-A', 'SARS-CoV', 
                                                "SARS-CoV-2 (pre-Omicron)", "SARS-CoV-2 (Omicron)"))) %>% 
  ungroup() %>% 
  arrange(pathogen) %>% 
  mutate(p_trans_before_mut = round(p_trans_before_mut, 2),
         p_trans_before_mut_lower_GT = round(p_trans_before_mut_lower_GT, 2),
         p_trans_before_mut_upper_GT = round(p_trans_before_mut_upper_GT, 2),
         p_trans_before_mut_lower_mut_rate = round(p_trans_before_mut_lower_mut_rate, 2),
         p_trans_before_mut_upper_mut_rate = round(p_trans_before_mut_upper_mut_rate, 2))



## Illustration of the simulation study for the different pathogens
list_plots <- lapply(1:length(list_sim_pathogens), FUN = function(i_pathogen){
  plot_from_sim_proba(sim_proba = list_sim_pathogens[[i_pathogen]], 
                      plot_title = df_p_trans_before_mut$pathogen[i_pathogen],
                      plot_subtitle = paste0('p = ',
                                             round(df_p_trans_before_mut$p_trans_before_mut[i_pathogen], 2) * 100,
                                             '%'),
                      xmax = 101)
})
panel_proba_sim <- ggarrange(plotlist = list_plots,
                             ncol = 3, nrow = 4,
                             common.legend = T, legend = 'top')


plot(panel_proba_sim)

# png('../figures/supplementary_figures/figure_s15_panel_sim_proba_all_pathogens.png',
#     height = 9, width = 7.6, res = 350, units = 'in')
# plot(panel_proba_sim)
# dev.off()
