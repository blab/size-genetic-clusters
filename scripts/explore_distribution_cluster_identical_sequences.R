source('utils_distrib.R')

library(tidyverse)
library(viridis)
library(metR)
library(RColorBrewer)

##### 1) Exploring the mean number of offsprings with identical genomes 
# for different assumptions regarding the reproduction number R, 
# and the probability p that 

### Defining the range of parameters to explore
vec_R0 <- seq(0.0, 2.2, 0.01)
vec_log_ratio <- seq(-3., 3., 0.001)

### Computing the mean number of offsprings with identical genomes (simply the product R * p)
df_mean_nb_offspring <- expand.grid(R0 = vec_R0, log_ratio = vec_log_ratio) %>% 
  mutate(ratio = exp(log_ratio),
         proba_trans_before_mut = get_proba_trans_before_mut(ratio)) %>% 
  mutate(mean_nb_offspring = R0 * proba_trans_before_mut)

### Displaying the results 
breaks_ratio <- c(0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0) # Breaks for the y-axis
log_breaks_ratio <- sapply(log(breaks_ratio), FUN = function(curr_log_break){
  vec_log_ratio[vec_log_ratio > curr_log_break][1]
})

plt_competition_mutation_transmission <- df_mean_nb_offspring %>% 
  ggplot(aes(x = R0, y = proba_trans_before_mut)) +
  geom_raster(aes(fill = mean_nb_offspring)) +
  geom_contour(aes(z = mean_nb_offspring), breaks = c(0.2, 0.6, 1.0, 1.4, 1.8), col = 'white') +
  geom_text_contour(aes(z = mean_nb_offspring),
                    nudge_x = 0.1, nudge_y = 0.1,
                    col = 'white', label.placer = label_placer_n(1)) +
  scale_fill_viridis(name = 'Mean number\nof offsprings\nwith identical\ngenomes', option = 'mako') +
  scale_y_continuous(name = expression(p),
                     trans = 'logit',
                     breaks = c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95),
                     expand = expansion(mult = c(0., 0.))) +
  scale_x_continuous(name = 'Reproduction number R',
                     expand = expansion(mult = c(0., 0.)))  +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = 'bottom') +
  guides(fill = guide_colorbar(barwidth = unit(1.7, 'in')))

# pdf('../figures/mean_nb_offspring_identical_sequences.pdf',
#     height = 3.5, width = 3.0)
# plot(plt_competition_mutation_transmission)
# dev.off()  

##### 2) Cumulative distribution function  
# for different assumptions regarding the reproduction number R, 
# and the dispersion parameter k

p <- 0.7
vec_R <- seq(0.5, 2.0, 0.5)
vec_k <- c(0.01, 0.1, 1.)

### Generating the distribution
df_theoretical_distrib <- Reduce('bind_rows', lapply(vec_R, FUN = function(R){
  Reduce('bind_rows', lapply(vec_k, FUN = function(k){
    theoretical_distrib_clusters(R, k, p, max_cluster_size = max_cluster_size) %>% 
      mutate(R = R, k = k, p = p)
  })) 
}))


### Plotting the distribution
plt_theoretical_distrib_cluster_size <- df_theoretical_distrib %>%
  mutate(k = paste0('k = ', k)) %>%
  ggplot(aes(x = cluster_size,
             y = cum_prob,
             group = interaction(k, R),
             colour = as.factor(R))) +
  geom_step() +
  facet_wrap(.~ k, nrow = 1) +
  scale_x_continuous(name = 'Size of cluster of identical sequences',
                     breaks = c(1, seq(10, 50, 10)),
                     expand = expansion(mult = c(0., 0.05))) +
  scale_y_continuous(name = 'Cumulative distribution\nfunction') +
  scale_colour_manual(name = 'R',
                      breaks = seq(0.5, 2.0, 0.5),
                      labels = c('0.5', '1.0', '1.5', '2.0'),
                      values = RColorBrewer::brewer.pal(5, 'BuPu')[-1]) +
  theme_bw() +
  coord_cartesian(xlim = c(1, 50), ylim = c(0.0, 1.0)) +
  theme(panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = 'gray22'),
        strip.text = element_text(colour = 'white'),
        panel.spacing = unit(1, 'lines'))

plot(plt_theoretical_distrib_cluster_size)

