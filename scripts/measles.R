library(doParallel)
library(foreach)

source('utils_inference.R')
source('utils_cluster_alloc.R')

## Input files
file_cluster_alloc <- '../data/measles/cluster_alloc_measles_pacenti.rds'
file_proba_trans_before_mut <- '../results/proba_trans_before_mut/p_measles.txt'

cluster_alloc <- readRDS(file_cluster_alloc)

## Probability that transmission occurs before mutation
p_trans_before_mut <- as.numeric(read.table(file_proba_trans_before_mut))

## Scenario for the proportion of infection sequenced
n_sequences <- sum(cluster_alloc$df_size_distrib$cluster_size * cluster_alloc$df_size_distrib$count)
n_suspected_cases <- 322 # See Pacenti et al for the 2017 outbreak
prop_cases_sequenced <- n_sequences/n_suspected_cases
vec_p_detect_cases <- c(0.5, 1.0) # Proportion of infections detected (as cases)

## Setting up the grid search
R0_min <- 0.01
R0_max <- 10.0
vec_R0 <- c(seq(R0_min, 1.0, 0.01), seq(1.1, R0_max, 0.1))
vec_k <- c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(1.0, 10.0, 0.1))

## Setting up the inference
max_cluster_size_inference <- 10000
n_cores <- 2

## Running the grid search
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)

t0 <- Sys.time()
df_grid_indep_obs <- Reduce('bind_rows', foreach(i_p_detect = 1:length(vec_p_detect_cases), .packages = c('tidyverse')) %dopar% {
  
  run_grid_search(vec_R = vec_R, vec_k = vec_k, p_detect = vec_p_detect_cases[i_p_detect] * prop_cases_sequenced, 
                  max_cluster_size_inference = max_cluster_size_inference,
                  cluster_alloc = cluster_alloc, p_trans_before_mut = p_trans_before_mut) %>% 
    mutate(p_detect = vec_p_detect_cases[i_p_detect])
  
})
t1 <- Sys.time()
print(t1 - t0)
stopCluster(cl)

## Getting CI from grid search
df_CI <- Reduce('bind_rows', lapply(vec_p_detect_cases, FUN = function(curr_p_detect){
  curr_df_grid <- df_grid_indep_obs %>%
    filter(p_detect == curr_p_detect)
  
  get_profile_likelihood_CI(curr_df_grid, col_loglik = 'log_lik') %>%
    mutate(p_detect = curr_p_detect)
}))

View(df_CI)

