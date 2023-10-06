library(doParallel)
library(foreach)
library(tidyverse)

source('utils_inference.R')
source('utils_cluster_alloc.R')

## Input files
file_cluster_alloc <- '../data/measles/cluster_alloc_measles_pacenti.rds'
file_proba_trans_before_mut <- '../results/proba_trans_before_mut/df_p_trans_before_mut_with_uncertainty.rds'

cluster_alloc <- readRDS(file_cluster_alloc)

## Probability that transmission occurs before mutation
p_trans_before_mut <- readRDS(file_proba_trans_before_mut) %>% 
  ungroup() %>% 
  filter(pathogen == 'Measles') %>% 
  select(p_trans_before_mut) %>%
  as.numeric()

## Scenario for the proportion of infection sequenced
n_sequences <- sum(cluster_alloc$df_size_distrib$cluster_size * cluster_alloc$df_size_distrib$count)
n_suspected_cases <- 322 # See Pacenti et al for the 2017 outbreak
prop_cases_sequenced <- n_sequences/n_suspected_cases
vec_p_detect_cases <- c(0.5, 1.0) # Proportion of infections detected (as cases)

## Setting up the boundary for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Setting up the inference
max_cluster_size_inference <- 10000
n_cores <- 2

## Running the grid search
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

t0 <- Sys.time()
df_inference <- Reduce('bind_rows',
                       foreach(i_p_detect = 1:length(vec_p_detect_cases),
                               .packages = c('tidyverse')) %dopar% {
                                 
                                 ## Define the function to minimize
                                 minus_log_lik <- function(param){
                                   R <- param[1]
                                   k <- param[2]
                                   
                                   ll <- get_loglik_obs(vec_cluster_size = cluster_alloc$df_size_distrib$cluster_size,
                                                        vec_n_clusters = cluster_alloc$df_size_distrib$count,
                                                        R = R, k = k,
                                                        p = p_trans_before_mut,
                                                        p_detect = vec_p_detect_cases[i_p_detect] * prop_cases_sequenced, 
                                                        max_cluster_size = max_cluster_size_inference)
                                   
                                   
                                   return(-ll)
                                 }
                                 
                                 mle_estim <- optim(par = c(1.0, 0.1),
                                                    fn = minus_log_lik,
                                                    method = 'L-BFGS-B',
                                                    lower = c(R_min, k_min), upper = c(R_max, k_max))
                                 
                                 CI_R <- get_CI(1, minus_log_lik, mle_estim,
                                                param_min = R_min, param_max = R_max,
                                                increment_param = 0.01)
                                 CI_k <- get_CI(2, minus_log_lik, mle_estim,
                                                param_min = k_min, param_max = k_max,
                                                increment_param = 0.001)
                                 
                                 bind_cols(tibble(mle_estim = mle_estim$par,
                                                  param = c('R', 'k'),
                                                  p_detect = rep(vec_p_detect_cases[i_p_detect], 2)),
                                           bind_rows(CI_R, CI_k))
                                 
                               })
t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)

print(df_inference)

