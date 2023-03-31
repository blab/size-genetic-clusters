library(tidyverse)

source('utils_inference.R')
source('utils_cluster_alloc.R')

args <- commandArgs(trailingOnly = T)
input_file <- as.character(args[1]) # Input file for the cluster size distribution
p_trans_before_mut <- as.numeric(args[2]) # Probability that transmission occurs before mutation
p_seq <- as.numeric(args[3]) # Proportion of infections sequenced
compute_CI <- as.logical(as.numeric(args[4]))
output_file <- as.character(args[5]) # Output file for the results of the inference

print(input_file)
print(p_trans_before_mut)
print(p_seq)
print(output_file)

## Load input file
cluster_alloc <- readRDS(input_file)

## Setting up the boundary for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Setting up the inference
max_cluster_size_inference <- 10000

## Define the function to minimize
minus_log_lik <- function(param){
  R <- param[1]
  k <- param[2]
  
  ll <- get_loglik_obs(vec_cluster_size = cluster_alloc$df_size_distrib$cluster_size,
                       vec_n_clusters = cluster_alloc$df_size_distrib$count,
                       R = R, k = k,
                       p = p_trans_before_mut,
                       p_detect = p_seq, 
                       max_cluster_size = max_cluster_size_inference)
  
  
  return(-ll)
}

## Get maximum likelihood estimates
mle_estim <- optim(par = c(1.0, 0.1),
                   fn = minus_log_lik,
                   method = 'L-BFGS-B',
                   lower = c(R_min, k_min), upper = c(R_max, k_max))

if(compute_CI){
  
  ## Get likelihood-profile confidence intervals
  CI_R <- get_CI(1, minus_log_lik, mle_estim,
                 param_min = R_min, param_max = R_max,
                 increment_param = 0.01)
  CI_k <- get_CI(2, minus_log_lik, mle_estim,
                 param_min = k_min, param_max = k_max,
                 increment_param = 0.001)
  
  ## Merge results
  df_res <- bind_cols(tibble(mle_estim = mle_estim$par,
                             param = c('R', 'k'),
                             p_seq = rep(p_seq, 2),
                             p_trans_before_mut = rep(p_trans_before_mut, 2)),
                      bind_rows(CI_R, CI_k))
  
} else{
  df_res <- tibble(mle_estim = mle_estim$par,
                   param = c('R', 'k'),
                   p_seq = rep(p_seq, 2),
                   p_trans_before_mut = rep(p_trans_before_mut, 2))
}

## Save results
write.csv(df_res, output_file)