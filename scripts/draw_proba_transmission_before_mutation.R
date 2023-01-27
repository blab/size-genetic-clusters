library(tidyverse)

source('utils_proba_trans_before_mut.R')

set.seed(765)
n_sim <- 1e7 

###############
### MERS
###############
# Mutation rate from clock regression in Dudas et al. eLife. 2018
rate_sub_per_site_per_year_mers <- 4.59 * 1e-4
length_genome_mers <- 30130

# Generation time distribution from Cauchemez et al. PNAS. 2016
mean_GT_mers <- 6.8 # (95% CrI: 6.0 - 7.8)
sd_GT_mers <- 6.3 # (95% CrI: 3.5 - 16.9)

sim_proba_mers <- compute_proba_transmission_before_mutation(n_sim,
                                                             mean_GT_mers, sd_GT_mers,
                                                             rate_sub_per_site_per_year_mers,
                                                             length_genome_mers)

print(sim_proba_mers$p_trans_before_mut)

plt_MERS <- plot_from_sim_proba(sim_proba_mers) + coord_cartesian(xlim = c(0, 50))
plt_MERS + theme(legend.position = 'bottom')

write.table(round(sim_proba_mers$p_trans_before_mut, 5), '../results/proba_trans_before_mut/p_mers.txt')

###############
### Measles
###############
# Mutation rate from 
rate_sub_per_site_per_year_measles <- 4.97e-4 # https://nextstrain.org/measles?l=clock 
length_genome_measles <- 15894

# Generation time distribution from https://reader.elsevier.com/reader/sd/pii/S0022519311003146?token=A6046F166054A3991F645E68FC6167652EB6BA1B02A9E20DDFB382FE39EE85903578B986777D7D0117CF44D0640FE02D&originRegion=us-east-1&originCreation=20221216223845
mean_GT_measles <- 11.2
sd_GT_measles <- 1.79

sim_proba_measles <- compute_proba_transmission_before_mutation(n_sim,
                                                             mean_GT_measles, sd_GT_measles,
                                                             rate_sub_per_site_per_year_measles,
                                                             length_genome_measles)

print(sim_proba_measles$p_trans_before_mut)

plt_measles <- plot_from_sim_proba(sim_proba_measles) +
  coord_cartesian(xlim = c(0, 50)) + 
  ggtitle('Measles')

plt_measles + theme(legend.position = 'bottom')

write.table(round(sim_proba_measles$p_trans_before_mut, 5), '../results/proba_trans_before_mut/p_measles.txt',
            row.names = F, col.names = F)


###############
### SARS-CoV-2 (pre-Omicron)
###############
# Mutation rate
rate_sub_per_year_covid <- 26.611 #https://nextstrain.org/ncov/gisaid/global/6m?l=clock

# Generation time
mean_GT_pre_omicron <- 5.9 # https://elifesciences.org/articles/70767.pdf
sd_GT_pre_omicron <- 4.8

sim_proba_pre_omicron <- compute_proba_transmission_before_mutation_rate_per_year(n_sim,
                                                                                  mean_GT_pre_omicron, sd_GT_pre_omicron,
                                                                                  rate_sub_per_year_covid)
print(sim_proba_pre_omicron$p_trans_before_mut)

plt_covid_pre_omicron <- plot_from_sim_proba(sim_proba_pre_omicron) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0., 0.17)) + 
  ggtitle('SARS-CoV-2 (pre-Omicron)')

plt_covid_pre_omicron + theme(legend.position = 'bottom')

write.table(round(sim_proba_pre_omicron$p_trans_before_mut, 5), '../results/proba_trans_before_mut/p_ncov_pre_omicron.txt',
            row.names = F, col.names = F)

###############
### SARS-CoV-2 (Omicron)
###############
# Mutation rate
rate_sub_per_year_covid <- 26.611 #https://nextstrain.org/ncov/gisaid/global/6m?l=clock

# Generation time
mean_GT_post_omicron <- mean_GT_pre_omicron - 1.0 # https://www.eurosurveillance.org/docserver/fulltext/eurosurveillance/27/6/eurosurv-27-6-1.pdf?expires=1671044970&id=id&accname=guest&checksum=FC76A4D5486FF82A96310E5A536CE70A + papier de Thomas
sd_GT_post_omicron <- sd_GT_pre_omicron

sim_proba_post_omicron <- compute_proba_transmission_before_mutation_rate_per_year(n_sim,
                                                                                   mean_GT_post_omicron, sd_GT_post_omicron,
                                                                                   rate_sub_per_year_covid)

print(sim_proba_post_omicron$p_trans_before_mut)

plt_covid_post_omicron <- plot_from_sim_proba(sim_proba_post_omicron) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0., 0.17)) +
  ggtitle('SARS-CoV-2 (post-Omicron)')

plt_covid_post_omicron + theme(legend.position = 'bottom')

write.table(round(sim_proba_post_omicron$p_trans_before_mut, 5), '../results/proba_trans_before_mut/p_ncov_post_omicron.txt',
            row.names = F, col.names = F)


##### All pathogens
panel_figure <- ggpubr::ggarrange(plt_MERS, plt_measles, 
                                  plt_covid_pre_omicron,
                                  plt_covid_post_omicron,
                                  common.legend = T, labels = 'AUTO')


plot(panel_figure)