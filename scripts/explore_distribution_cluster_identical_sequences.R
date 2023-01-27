source('utils_distrib.R')

library(tidyverse)
library(viridis)
library(metR)

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
  ggplot(aes(y = log_ratio, x = R0)) +
  geom_raster(aes(fill = mean_nb_offspring)) +
  geom_contour(aes(z = mean_nb_offspring), breaks = c(0.2, 0.6, 1.0, 1.4, 1.8), col = 'white') +
  geom_text_contour(aes(z = mean_nb_offspring),
                    nudge_x = 0.1, nudge_y = 0.1,
                    col = 'white', label.placer = label_placer_n(1)) +
  scale_fill_viridis(name = 'Mean number of offsprings\nwith identical genomes', option = 'mako') +
  scale_y_continuous(name = expression(d[mut]/d[gen]),
                     breaks = log_breaks_ratio, labels = breaks_ratio,
                     expand = expansion(mult = c(0., 0.))) +
  scale_x_continuous(name = 'Reproduction number R',
                     expand = expansion(mult = c(0., 0.)))  +
  theme_classic() +
  theme(axis.line = element_blank())
  
plot(plt_competition_mutation_transmission)

##### 2) Cumulative distribution function  
# for different assumptions regarding the reproduction number R, 
# and the probability p that 