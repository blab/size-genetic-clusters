library(tidyverse)
library(ggpubr)
library(ggrepel)

source('utils_proba_trans_before_mut.R')

## Load values for the mutation rate and the generation time
df_pathogens <- read_csv('../data/proba_trans_before_mut/gt_mutation_rate.csv')

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

## Run the simulations for all pathogens
n_sim <- 1e7
set.seed(9876)
list_sim_pathogens <-  lapply(1:nrow(df_pathogens), FUN = function(i_pathogen){
  sim_pathogen <- compute_proba_transmission_before_mutation_rate_per_day(n_sim = n_sim,
                                                                          mean_GT = df_pathogens$mean_GT[i_pathogen],
                                                                          sd_GT = df_pathogens$sd_GT[i_pathogen],
                                                                          rate_sub_per_day = df_pathogens$subs_per_day[i_pathogen])
})

## Reproduction number threshold for different pathogens
df_p_trans_before_mut <- Reduce('bind_rows', lapply(1:length(list_sim_pathogens), FUN = function(i_pathogen){
  c('pathogen' = df_pathogens$pathogen[i_pathogen],
    'p_trans_before_mut' = list_sim_pathogens[[i_pathogen]]$p_trans_before_mut)
})) %>% 
  mutate(p_trans_before_mut = as.numeric(p_trans_before_mut),
         R_max = 1./p_trans_before_mut)

print(df_p_trans_before_mut)

df_curve_to_plot <- tibble(p_trans_before_mut = seq(0.01, 1.0, 0.01)) %>% 
  mutate(R_max = 1./p_trans_before_mut)

plt_trans_before_mut <- df_p_trans_before_mut %>% 
  filter(pathogen != 'SARS-CoV') %>% 
  ggplot(aes(x = p_trans_before_mut,
             y = R_max)) +
  geom_hline(colour = 'darkgrey', yintercept = 1.0, linetype = 'dashed') +
  geom_line(data = df_curve_to_plot) +
  geom_point(colour = 'darkcyan', size = 1.0) +
  geom_label_repel(aes(label = pathogen), box.padding = 1.7, max.overlaps = Inf) +
  scale_x_continuous(name = 'Probability that transmission\noccurs before mutation',
                     expand = expansion(mult = c(0., 0.)),
                     breaks = seq(0., 1.0, 0.1)) +
  scale_y_continuous(name = 'Reproduction number threshold',
                     expand = expansion(mult = c(0., 0.05)),
                     breaks = seq(0., 3.0, 0.5)) +
  coord_cartesian(xlim = c(0.3, 1.0), ylim = c(0.8, 3.0)) +
  theme_classic()

plot(plt_trans_before_mut)

## Illustration of the simulation study for the different pathogens
list_plots <- lapply(1:length(list_sim_pathogens), FUN = function(i_pathogen){
  plot_from_sim_proba(sim_proba = list_sim_pathogens[[i_pathogen]], 
                      plot_title = df_p_trans_before_mut$pathogen[i_pathogen],
                      plot_subtitle = paste0('p = ',
                                             round(df_p_trans_before_mut$p_trans_before_mut[i_pathogen], 2) * 100,
                                             '%'))
})
panel_proba_sim <- ggarrange(plotlist = list_plots,
                             ncol = 3, nrow = 4,
                             common.legend = T, legend = 'top')


plot(panel_proba_sim)
