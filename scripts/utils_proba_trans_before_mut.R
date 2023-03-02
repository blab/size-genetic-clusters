compute_proba_transmission_before_mutation_rate_per_day <- function(n_sim,
                                                                    mean_GT, sd_GT,
                                                                    rate_sub_per_day){
  
  ## Inputs: 
  ## n_sim: number of replicates
  ## mean_GT: mean generation time (in days)
  ## sd_GT: standard deviation of the generation time (in days)
  ## rate_sub_per_day: rate of substitution per day
  
  ## Draw vector of time to mutation
  vec_time_to_mut <- rexp(n = n_sim, rate = rate_sub_per_day)
  
  ## Define characteristics of the GT distribution (assumed Gamma)
  alpha_GT <- (mean_GT^2) / (sd_GT)^2
  beta_GT <- mean_GT / (sd_GT)^2
  
  ## Draw vector of transmission events
  vec_time_to_transmission <- rgamma(n = n_sim, shape = alpha_GT, rate = beta_GT)
  
  ## Empirical probability that transmission occurs before mutation
  p_trans_before_mut <- sum(vec_time_to_transmission < vec_time_to_mut) / n_sim
  
  return(list(p_trans_before_mut = p_trans_before_mut,
              vec_time_to_transmission = vec_time_to_transmission,
              vec_time_to_mut = vec_time_to_mut))
}

## Script to plot the simulated time distributions obtained with one of the above functions

plot_from_sim_proba <- function(sim_proba){
  col_mutation <- 'orange2'
  col_transmission <- 'darkslateblue'
  
  plt_sim <- tibble(time_trans = sim_proba$vec_time_to_transmission,
                    time_mut = sim_proba$vec_time_to_mut) %>% 
    ggplot() +
    geom_histogram(aes(x = time_trans,
                       y = after_stat(count / sum(count)),
                       fill = 'Time to transmission'),  
                   alpha= 0.5, binwidth = 1) +
    geom_histogram(aes(x = time_mut, 
                       y = after_stat(count / sum(count)),
                       fill = 'Time to mutation'),
                   alpha = 0.5, binwidth = 1) +
    theme_classic() +
    scale_fill_manual(name = '',
                      breaks = c('Time to transmission', 'Time to mutation'),
                      values = c(col_transmission, col_mutation)) +
    scale_y_continuous(name = 'Probability',
                       expand = expansion(mult = c(0., 0.02))) +
    scale_x_continuous(name = 'Time (in days)',
                       expand = expansion(mult = c(0.03, 0.05)))
  
  return(plt_sim)
}