library(tidyverse)
library(doParallel)
library(foreach)
library(ggpubr)

source('utils_distrib.R')
source('utils_inference.R')

# Define the characteristics used to generate synthetic clusters
curr_k <- 0.1 # Dispersion parameter k
curr_R <- 1.0 # Reproduction number R
curr_seed <- 100 # Seed used for the random number generation
curr_max_cluster_size_inference <- 10000 # Maximum cluster size used in the inference (c_max)
curr_ratio_delay <- 1.0 # Ratio between the average delay before mutation and transmission
curr_p <- curr_ratio_delay / (1. + curr_ratio_delay) # Probability that transmission occurs before mutation
curr_p_detect <- 0.1 # Proportion of infections sequenced
curr_max_cluster_size_sim <- 10000 # Maximum cluster size used to generate synthetic clusters

# Generate synthetic clusters
set.seed(curr_seed)
## Simulate branching process
df_res <- simulate_branching_process_until_max_cluster_size_save_largest_gen(50000, curr_R, curr_k, curr_p,
                                                                             max_cluster_size = curr_max_cluster_size_sim)

## Simulate observation process
df_res$cluster_size_obs <- rbinom(n = 50000, size = df_res$tot_cluster_size, prob = curr_p_detect)
df_res <- df_res %>% filter(cluster_size_obs > 0)

## Define the observation datasets
vec_size_dataset <- c(50, 100, 1000, 5000)
list_df_res <- lapply(vec_size_dataset, FUN = function(curr_size_obs){
  if(curr_size_obs > nrow(df_res)){
    NULL
  } else{
    df_res[1:curr_size_obs,] %>% 
      group_by(cluster_size_obs) %>% 
      summarise(n_clusters = n())
  }
})
if(Reduce('+', lapply(list_df_res, is.null)) > 0){
  print('NEED TO INCREASE THE NUMBER OF SIMULATIONS BEING PERFORMED')
}
list_max_cluster_size_inference <- lapply(list_df_res, FUN = function(curr_df_res){
  max(c(curr_max_cluster_size_inference, curr_df_res$cluster_size_obs + 1))
})


# Running our inference framework on our simulated data
n_cores <- 4
## Setting up the bounds for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Inference
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

t0 <- Sys.time()
df_inference <- foreach(i_size_dataset = 1:length(vec_size_dataset),
                        .packages = c('tidyverse'),
                        .combine = 'bind_rows') %dopar% {
                          
                          minus_log_lik <- function(param){
                            R <- param[1]
                            k <- param[2]
                            
                            curr_log_lik <- get_loglik_obs(vec_cluster_size = list_df_res[[i_size_dataset]]$cluster_size_obs,
                                                           vec_n_clusters = list_df_res[[i_size_dataset]]$n_clusters,
                                                           R = R, k = k,
                                                           p = curr_p, p_detect = curr_p_detect,
                                                           max_cluster_size = list_max_cluster_size_inference[[i_size_dataset]])
                            
                            
                            return(- (curr_log_lik))
                          }
                          
                          mle_estimate <- optim(par = c(1.0, 0.1),
                                                fn = minus_log_lik,
                                                method = 'L-BFGS-B',
                                                lower = c(R_min, k_min), upper = c(R_max, k_max))
                          
                          
                          ## Get CI
                          CI_R <- get_CI_from_vec_param(1, minus_log_lik, mle_estimate,
                                                        vec_param = seq(R_min, R_max, 0.01))
                          CI_k <- get_CI_from_vec_param(2, minus_log_lik, mle_estimate,
                                                        vec_param = c(seq(k_min, 0.11, 0.001),
                                                                      seq(0.12, 1.1, 0.01),
                                                                      seq(1.12, k_max, 0.1)))
                          
                          bind_cols(tibble(mle_estim = mle_estimate$par,
                                           param = c('R', 'k')),
                                    bind_rows(CI_R, CI_k)) %>% 
                            mutate(max_log_lik = - mle_estimate$value,
                                   i_size_dataset = i_size_dataset,
                                   size_dataset = vec_size_dataset[i_size_dataset],
                                   CV_status = mle_estimate$convergence)
                          
                        }
t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)

# Plot the results of the inference
plt_R <- df_inference %>% 
  filter(param == 'R') %>% 
  ggplot(aes(x = as.factor(size_dataset), y = mle_estim)) +
  geom_hline(yintercept = curr_R, linetype = 'dashed', color = 'darkgrey') +
  geom_point() +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95)) +
  scale_y_continuous(name = 'R estimate',
                     limits = c(0., NA), expand = expansion(mult = c(0., 0.05))) +
  scale_x_discrete(name = 'Number of clusters') +
  theme_classic()
  
plt_k <- df_inference %>% 
  filter(param == 'k') %>% 
  ggplot(aes(x = as.factor(size_dataset), y = mle_estim)) +
  geom_hline(yintercept = curr_k, linetype = 'dashed', color = 'darkgrey') +
  geom_point() +
  geom_linerange(aes(ymin = lower_95, ymax = upper_95)) +
  scale_y_continuous(name = 'k estimate',
                     trans = 'log',
                     breaks = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0),
                     limits = c(0.001, 10.),
                     expand = expansion(mult = c(0., 0.))) +
  scale_x_discrete(name = 'Number of clusters') +
  theme_classic()

ggarrange(plt_R, plt_k, ncol = 2)
