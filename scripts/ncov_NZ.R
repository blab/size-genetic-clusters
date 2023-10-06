library(doParallel)
library(foreach)
library(tidyverse)

source('utils_inference.R')
source('utils_cluster_alloc.R')

## Input files
file_cluster_alloc <- '../data/ncov_NZ/df_cluster_by_period_NZ.rds'
file_prop_cases_sequenced_per_period <- '../data/ncov_NZ/df_prop_sequenced_per_period.rds'
file_proba_trans_before_mut <- '../results/proba_trans_before_mut/df_p_trans_before_mut_with_uncertainty.rds'

## Probability that transmission occurs before mutation
p_trans_before_mut_pre_omicron <- readRDS(file_proba_trans_before_mut) %>% 
  ungroup() %>% 
  filter(pathogen == 'SARS-CoV-2 (pre-Omicron)') %>% 
  select(p_trans_before_mut) %>%
  as.numeric()
p_trans_before_mut_post_omicron <- readRDS(file_proba_trans_before_mut) %>% 
  ungroup() %>% 
  filter(pathogen == 'SARS-CoV-2 (Omicron)') %>% 
  select(p_trans_before_mut) %>%
  as.numeric()

## Cluster size distribution
df_clusters_by_period <- readRDS(file_cluster_alloc)

vec_periods <- sort(unique(df_clusters_by_period$period))
n_periods <- length(vec_periods)

list_curr_df_clusters_pre_omicron <- lapply(vec_periods, FUN = function(curr_period){
  df_clusters_by_period %>% 
    filter(period == curr_period) %>% filter(! is_omicron)
})
list_curr_df_clusters_post_omicron <- lapply(vec_periods, FUN = function(curr_period){
  df_clusters_by_period %>% 
    filter(period == curr_period) %>% filter(is_omicron)
})

## Data describing the fraction of cases sequenced over time
df_prop_sequenced_per_period <- readRDS(file_prop_cases_sequenced_per_period)
vec_prop_cases_sequenced <- sapply(vec_periods, FUN = function(curr_period){
  prop_cases_sequenced <- as.numeric(unlist(df_prop_sequenced_per_period[df_prop_sequenced_per_period$period == curr_period, 'prop_cases_sequenced']))
})

## Scenario for the proportion of infection detected
vec_p_detect_cases <- c(0.5, 0.8, 1.0) # Proportion of infections detected (as cases)

## Setting up the bounds for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Setting up the inference
max_cluster_size_inference <- 10000
n_cores <- 3

## Running the inference
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
t0 <- Sys.time()
df_inference <- Reduce('bind_rows', foreach(i_p_detect = 1:length(vec_p_detect_cases), .packages = c('tidyverse')) %dopar% {
  curr_vec_p_detect_by_period <- vec_p_detect_cases[i_p_detect] * vec_prop_cases_sequenced
  
  ## Define functions to be minimized
  minus_log_lik <- function(param){
    vec_R <- param[1:n_periods]
    k <- param[n_periods + 1]
    
    vec_log_lik_by_period <- sapply(1:n_periods, FUN = function(i_period){
      curr_df_clusters_pre_omicron <- list_curr_df_clusters_pre_omicron[[i_period]]
      curr_df_clusters_post_omicron <- list_curr_df_clusters_post_omicron[[i_period]]
      
      
      if(nrow(curr_df_clusters_pre_omicron) == 0){
        curr_log_lik_pre_omicron <- 0.0
      } else{
        curr_log_lik_pre_omicron <- get_loglik_obs(vec_cluster_size = curr_df_clusters_pre_omicron$cluster_size,
                                                   vec_n_clusters = curr_df_clusters_pre_omicron$n_clusters,
                                                   R = vec_R[i_period], k = k,
                                                   p = p_trans_before_mut_pre_omicron, p_detect = curr_vec_p_detect_by_period[i_period],
                                                   max_cluster_size = max_cluster_size_inference)
      }
      
      if(nrow(curr_df_clusters_post_omicron) == 0){
        curr_log_lik_post_omicron <- 0.0
      } else{
        curr_log_lik_post_omicron <- get_loglik_obs(vec_cluster_size = curr_df_clusters_post_omicron$cluster_size,
                                                    vec_n_clusters = curr_df_clusters_post_omicron$n_clusters,
                                                    R = vec_R[i_period], k = k,
                                                    p = p_trans_before_mut_post_omicron, p_detect = curr_vec_p_detect_by_period[i_period],
                                                    max_cluster_size = max_cluster_size_inference)
      }
      
      (curr_log_lik_pre_omicron + curr_log_lik_post_omicron)
    })
    
    return(- sum(vec_log_lik_by_period))
    
  }
  
  mle_estim <- optim(par = c(rep(1.0, n_periods), 0.1),
                     fn = minus_log_lik,
                     method = 'L-BFGS-B',
                     lower = c(rep(R_min, n_periods), k_min), upper = c(rep(R_max, n_periods), k_max))
  
  list_CI_R <- Reduce('bind_rows', lapply(1:n_periods, FUN = function(i_period){
    CI_R <- get_CI_from_vec_param(i_period, minus_log_lik, mle_estim,
                   vec_param = seq(R_min, R_max, 0.01))
  }))
  
  CI_k <- get_CI_from_vec_param(n_periods + 1, minus_log_lik, mle_estim,
                 vec_param = seq(k_min, k_max, 0.001))
  
  
  
  bind_cols(tibble(mle_estim = mle_estim$par,
                   param = c(paste0('R_', 1:n_periods), 'k'),
                   period = c(vec_periods, 'All'),
                   p_detect = rep(vec_p_detect_cases[i_p_detect], n_periods + 1)),
            bind_rows(list_CI_R, CI_k))
})

t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)

## Plotting the results
print(df_inference)
