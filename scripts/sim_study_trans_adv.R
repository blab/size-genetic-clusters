library(tidyverse)
library(doParallel)
library(foreach)

source('utils_distrib.R')
source('utils_inference.R')

# Define the characteristics used to generate synthetic clusters
curr_baseline_R <- 0.75 # Reproduction number of the non-variant R_NV
curr_increase_trans <- 0.5 # Transmission advantage of the variant
curr_variant_R <- curr_baseline_R * (1.0 + curr_increase_trans) # Reproduction number of the variant R_V
curr_k <- 0.1 # Dispersion parameter k
curr_seed <- 1 # Seed used for the random number generation
curr_ratio_delay <- 1.0 # Ratio between the average delay before mutation and transmission
curr_p <- curr_ratio_delay / (1. + curr_ratio_delay) # Probability that transmission occurs before mutation
curr_p_detect <- 0.5 # Proportion of infections sequenced
max_cluster_size_inference <- 10000 # Maximum cluster size used in the inference (c_max)
curr_max_cluster_size_sim <- 100000  # Maximum cluster size used to generate synthetic clusters

# Generate synthetic clusters
set.seed(curr_seed)

## Simulate branching process
df_res_non_variant <- simulate_branching_process_until_max_cluster_size_save_largest_gen(50000, curr_baseline_R, curr_k, 
                                                                                         curr_p,
                                                                                         max_cluster_size = curr_max_cluster_size_sim)
df_res_variant <- simulate_branching_process_until_max_cluster_size_save_largest_gen(50000, curr_variant_R, curr_k, 
                                                                                     curr_p,
                                                                                     max_cluster_size = curr_max_cluster_size_sim)

## Simulate observation process
df_res_non_variant$cluster_size_obs <- rbinom(n = 50000, size = df_res_non_variant$tot_cluster_size, prob = curr_p_detect)
df_res_variant$cluster_size_obs <- rbinom(n = 50000, size = df_res_variant$tot_cluster_size, prob = curr_p_detect)

df_res_non_variant <- df_res_non_variant %>% filter(cluster_size_obs > 0)
df_res_variant <- df_res_variant %>% filter(cluster_size_obs > 0)

## Define the observation datasets
vec_size_dataset <- c(50, 100, 1000, 5000)
list_df_res_non_variant <- lapply(vec_size_dataset, FUN = function(curr_size_obs){
  df_res_non_variant[1:curr_size_obs,] %>% 
    group_by(cluster_size_obs) %>% 
    summarise(n_clusters = n())
})
list_df_res_variant <- lapply(vec_size_dataset, FUN = function(curr_size_obs){
  df_res_variant[1:curr_size_obs,]  %>% 
    group_by(cluster_size_obs) %>% 
    summarise(n_clusters = n())
})

if(Reduce('+', lapply(list_df_res_non_variant, is.null)) > 0){
  print('NEED TO INCREASE THE NUMBER OF SIMULATIONS BEING PERFORMED')
}
if(Reduce('+', lapply(list_df_res_variant, is.null)) > 0){
  print('NEED TO INCREASE THE NUMBER OF SIMULATIONS BEING PERFORMED')
}

# Running our inference framework on our simulated data
n_cores <- 4
## Setting up the bounds for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Quantile vector used to compute the p-value
vec_p_quantile <- seq(0.0, 1.0, 1e-5)
vec_quantile_function_chisq <- qchisq(vec_p_quantile, df = 1)

## Inference
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

t0 <- Sys.time()
df_LRT <- foreach(i_size_dataset = 1:length(vec_size_dataset),
                  .packages = c('tidyverse'),
                  .combine = 'bind_rows') %dopar% {
                    
                    minus_log_lik_combined <- function(param){
                      R <- param[1]
                      k <- param[2]
                      
                      curr_log_lik_non_var <- get_loglik_obs(vec_cluster_size = list_df_res_non_variant[[i_size_dataset]]$cluster_size_obs,
                                                             vec_n_clusters = list_df_res_non_variant[[i_size_dataset]]$n_clusters,
                                                             R = R, k = k,
                                                             p = curr_p, p_detect = curr_p_detect,
                                                             max_cluster_size = max(c(max_cluster_size_inference, list_df_res_non_variant[[i_size_dataset]]$cluster_size_obs)))
                      
                      curr_log_lik_var <- get_loglik_obs(list_df_res_variant[[i_size_dataset]]$cluster_size_obs,
                                                         list_df_res_variant[[i_size_dataset]]$n_clusters,
                                                         R, k,
                                                         p = curr_p, p_detect = curr_p_detect,
                                                         max_cluster_size = max(c(max_cluster_size_inference, list_df_res_variant[[i_size_dataset]]$cluster_size_obs)))
                      
                      return(- (curr_log_lik_non_var + curr_log_lik_var))
                    }
                    minus_log_lik_split <- function(param){
                      R_var <- param[1]
                      R_non_var <- param[2]
                      k <- param[3]
                      
                      curr_log_lik_non_var <- get_loglik_obs(list_df_res_non_variant[[i_size_dataset]]$cluster_size_obs,
                                                             list_df_res_non_variant[[i_size_dataset]]$n_clusters,
                                                             R_non_var, k,
                                                             p = curr_p, p_detect = curr_p_detect,
                                                             max_cluster_size = max(c(max_cluster_size_inference, list_df_res_non_variant[[i_size_dataset]]$cluster_size_obs)))
                      
                      curr_log_lik_var <- get_loglik_obs(list_df_res_variant[[i_size_dataset]]$cluster_size_obs,
                                                         list_df_res_variant[[i_size_dataset]]$n_clusters,
                                                         R_var, k,
                                                         p = curr_p, p_detect = curr_p_detect,
                                                         max_cluster_size = max(c(max_cluster_size_inference, list_df_res_variant[[i_size_dataset]]$cluster_size_obs)))
                      
                      return(- (curr_log_lik_non_var + curr_log_lik_var))
                    }
                    
                    mle_combined <- optim(par = c(1.0, 0.1),
                                          fn = minus_log_lik_combined,
                                          method = 'L-BFGS-B',
                                          lower = c(R_min, k_min), upper = c(R_max, k_max))
                    
                    mle_split <- optim(par = c(1.0, 1.0, 0.1),
                                       fn = minus_log_lik_split,
                                       method = 'L-BFGS-B',
                                       lower = c(R_min, R_min, k_min), upper = c(R_max, R_max, k_max))
                    
                    
                    ## Computate likelihood ratio test p-value
                    LRT <- as.numeric(2 * (mle_combined$value - mle_split$value))
                    curr_p_val <- 1.0 - vec_p_quantile[which(vec_quantile_function_chisq > LRT)[1]]
                    
                    c('R_MLE_combined' = mle_combined$par[1],
                      'k_MLE_combined' = mle_combined$par[2],
                      'max_log_lik_combined' = - mle_combined$value,
                      'R_var_MLE_split' = mle_split$par[1],
                      'R_non_var_MLE_split' = mle_split$par[2],
                      'k_MLE_split' = mle_split$par[3],
                      'max_log_lik_split' = - mle_split$value,
                      'LRT' = LRT,
                      'p_val' = curr_p_val,
                      'i_size_dataset' = i_size_dataset,
                      'size_dataset' = vec_size_dataset[i_size_dataset],
                      'CV_combined' = mle_combined$convergence,
                      'CV_split' = mle_split$convergence
                    )
                    
                  }
t1 <- Sys.time()
print(t1 - t0)

stopCluster(cl)

# Look at the results of the inference
df_LRT %>% 
  mutate(trans_adv_estim = R_var_MLE_split/R_non_var_MLE_split - 1.0) %>% 
  select(size_dataset, p_val, trans_adv_estim)



