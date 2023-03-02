library(doParallel)
library(foreach)
library(tidyverse)

source('utils_inference.R')
source('utils_cluster_alloc.R')

name_variant <- 'Alpha'

## Input files
file_input_WA_analysis <- '../data/ncov_WA/input_WA_analysis.csv'
file_input_cluster_dist <- paste0('../data/ncov_WA/cluster_distrib_', name_variant, '.rds')
file_proba_trans_before_mut_post_omicron <- '../results/proba_trans_before_mut/p_ncov_post_omicron.txt'
file_proba_trans_before_mut_pre_omicron <- '../results/proba_trans_before_mut/p_ncov_pre_omicron.txt'

df_scenarios <- read_csv(file_input_WA_analysis)
if(! name_variant %in% df_scenarios$variant){
  print('Need to pick a variant whose name is in the df_scenarios dataframe!')
}

## Probability that transmission occurs before mutation
p_trans_before_mut_pre_omicron <- as.numeric(read.table(file_proba_trans_before_mut_pre_omicron))
p_trans_before_mut_post_omicron <-as.numeric(read.table(file_proba_trans_before_mut_post_omicron))


## Datasets for the analysis
list_df_cluster_distrib <- readRDS(file_input_cluster_dist)

vec_Nextstrain_clades_Omicron <- c("21K (Omicron)", "21L (Omicron)", 
                                   "21M (Omicron)", 
                                   "22A (Omicron)", "22B (Omicron)",
                                   "22C (Omicron)", "22D (Omicron)",
                                   "22E (Omicron)", "22F (Omicron)")

df_clusters_var_pre_omicron <- list_df_cluster_distrib$df_cluster_distrib %>% 
  filter(is_variant) %>% 
  filter(! Nextstrain_clade %in% vec_Nextstrain_clades_Omicron) %>% 
  group_by(n_in_same_clique) %>% 
  summarise(count = sum(n_clusters))

df_clusters_var_post_omicron <- list_df_cluster_distrib$df_cluster_distrib %>% 
  filter(is_variant) %>% 
  filter(Nextstrain_clade %in% vec_Nextstrain_clades_Omicron) %>% 
  group_by(n_in_same_clique) %>% 
  summarise(count = sum(n_clusters))

df_clusters_non_var_pre_omicron <- list_df_cluster_distrib$df_cluster_distrib %>% 
  filter(! is_variant) %>% 
  filter(! Nextstrain_clade %in% vec_Nextstrain_clades_Omicron) %>% 
  group_by(n_in_same_clique) %>% 
  summarise(count = sum(n_clusters))

df_clusters_non_var_post_omicron <- list_df_cluster_distrib$df_cluster_distrib %>% 
  filter(! is_variant) %>% 
  filter(Nextstrain_clade %in% vec_Nextstrain_clades_Omicron) %>% 
  group_by(n_in_same_clique) %>% 
  summarise(count = sum(n_clusters))

## Scenario for the proportion of infections detected
vec_p_detect_cases <-  c(0.1, 0.2, 0.5, 0.8) # Proportion of infections detected (as cases)
curr_prop_cases_sequenced <- # Proportion of cases sequenced during timewindow
  list_df_cluster_distrib$n_seq_during_timewindow / list_df_cluster_distrib$n_cases_during_timewindow
vec_p_detect <- vec_p_detect_cases * curr_prop_cases_sequenced # Scenarios for the proportion of infections sequenced

## Setting up the bounds for the optimization
R_min <- 0.01
R_max <- 10.0
k_min <- 0.001
k_max <- 10.0

## Setting up the inference
max_cluster_size_inference <- 10000
n_cores <- 4

vec_p_quantile <- c(0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11,
                    1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, seq(2e-5, 1.0, 1e-5))
vec_quantile_function_chisq <- qchisq(vec_p_quantile, df = 1)

## Running the inference
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
df_inference <- Reduce('bind_rows', foreach(i_p_detect = 1:length(vec_p_detect_cases), .packages = c('tidyverse')) %dopar% {
  
  curr_p_detect <- vec_p_detect[i_p_detect]
  
  ## Define functions to be minimized
  minus_log_lik_combined <- function(param){
    R <- param[1]
    k <- param[2]
    
    if(nrow(df_clusters_non_var_pre_omicron) == 0){
      curr_log_lik_non_var_pre_omicron <- 0.0
    } else{
      curr_log_lik_non_var_pre_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_non_var_pre_omicron$n_in_same_clique,
                                                         vec_n_clusters = df_clusters_non_var_pre_omicron$count,
                                                         R = R, k = k,
                                                         p = p_trans_before_mut_pre_omicron, p_detect = curr_p_detect,
                                                         max_cluster_size = max_cluster_size_inference)
    }
    
    if(nrow(df_clusters_non_var_post_omicron) == 0){
      curr_log_lik_non_var_post_omicron <- 0.0
    } else{
      curr_log_lik_non_var_post_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_non_var_post_omicron$n_in_same_clique,
                                                          vec_n_clusters = df_clusters_non_var_post_omicron$count,
                                                          R = R, k = k,
                                                          p = p_trans_before_mut_post_omicron, p_detect = curr_p_detect,
                                                          max_cluster_size = max_cluster_size_inference)
    }
    
    if(nrow(df_clusters_var_pre_omicron) == 0){
      curr_log_lik_var_pre_omicron <- 0.0
    } else{
      curr_log_lik_var_pre_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_var_pre_omicron$n_in_same_clique,
                                                     vec_n_clusters = df_clusters_var_pre_omicron$count,
                                                     R = R, k = k,
                                                     p = p_trans_before_mut_pre_omicron, p_detect = curr_p_detect,
                                                     max_cluster_size = max_cluster_size_inference)
      
    }
    
    if(nrow(df_clusters_var_post_omicron) == 0){
      curr_log_lik_var_post_omicron <- 0.0
      
    } else{
      curr_log_lik_var_post_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_var_post_omicron$n_in_same_clique,
                                                      vec_n_clusters = df_clusters_var_post_omicron$count,
                                                      R = R, k = k,
                                                      p = p_trans_before_mut_post_omicron, p_detect = curr_p_detect,
                                                      max_cluster_size = max_cluster_size_inference)
      
    }
    
    return(- (curr_log_lik_non_var_pre_omicron + curr_log_lik_non_var_post_omicron +
                curr_log_lik_var_pre_omicron + curr_log_lik_var_post_omicron))
  }
  
  minus_log_lik_split <- function(param){
    R_var <- param[1]
    R_non_var <- param[2]
    k <- param[3]
    
    if(nrow(df_clusters_non_var_pre_omicron) == 0){
      curr_log_lik_non_var_pre_omicron <- 0.0
    } else{
      curr_log_lik_non_var_pre_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_non_var_pre_omicron$n_in_same_clique,
                                                         vec_n_clusters = df_clusters_non_var_pre_omicron$count,
                                                         R = R_non_var, k = k,
                                                         p = p_trans_before_mut_pre_omicron, p_detect = curr_p_detect,
                                                         max_cluster_size = max_cluster_size_inference)
    }
    
    if(nrow(df_clusters_non_var_post_omicron) == 0){
      curr_log_lik_non_var_post_omicron <- 0.0
    } else{
      curr_log_lik_non_var_post_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_non_var_post_omicron$n_in_same_clique,
                                                          vec_n_clusters = df_clusters_non_var_post_omicron$count,
                                                          R = R_non_var, k = k,
                                                          p = p_trans_before_mut_post_omicron, p_detect = curr_p_detect,
                                                          max_cluster_size = max_cluster_size_inference)
    }
    
    if(nrow(df_clusters_var_pre_omicron) == 0){
      curr_log_lik_var_pre_omicron <- 0.0
    } else{
      curr_log_lik_var_pre_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_var_pre_omicron$n_in_same_clique,
                                                     vec_n_clusters = df_clusters_var_pre_omicron$count,
                                                     R = R_var, k = k,
                                                     p = p_trans_before_mut_pre_omicron, p_detect = curr_p_detect,
                                                     max_cluster_size = max_cluster_size_inference)
      
    }
    
    if(nrow(df_clusters_var_post_omicron) == 0){
      curr_log_lik_var_post_omicron <- 0.0
      
    } else{
      curr_log_lik_var_post_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_var_post_omicron$n_in_same_clique,
                                                      vec_n_clusters = df_clusters_var_post_omicron$count,
                                                      R = R_var, k = k,
                                                      p = p_trans_before_mut_post_omicron, p_detect = curr_p_detect,
                                                      max_cluster_size = max_cluster_size_inference)
      
    }
    
    return(- (curr_log_lik_non_var_pre_omicron + curr_log_lik_non_var_post_omicron +
                curr_log_lik_var_pre_omicron + curr_log_lik_var_post_omicron))
  }
  
  minus_log_lik_fully_split <- function(param){
    R_var <- param[1]
    R_non_var <- param[2]
    k_var <- param[3]
    k_non_var <- param[4]
    
    if(nrow(df_clusters_non_var_pre_omicron) == 0){
      curr_log_lik_non_var_pre_omicron <- 0.0
    } else{
      curr_log_lik_non_var_pre_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_non_var_pre_omicron$n_in_same_clique,
                                                         vec_n_clusters = df_clusters_non_var_pre_omicron$count,
                                                         R = R_non_var, k = k_non_var,
                                                         p = p_trans_before_mut_pre_omicron, p_detect = curr_p_detect,
                                                         max_cluster_size = max_cluster_size_inference)
    }
    
    if(nrow(df_clusters_non_var_post_omicron) == 0){
      curr_log_lik_non_var_post_omicron <- 0.0
    } else{
      curr_log_lik_non_var_post_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_non_var_post_omicron$n_in_same_clique,
                                                          vec_n_clusters = df_clusters_non_var_post_omicron$count,
                                                          R = R_non_var, k = k_non_var,
                                                          p = p_trans_before_mut_post_omicron, p_detect = curr_p_detect,
                                                          max_cluster_size = max_cluster_size_inference)
    }
    
    if(nrow(df_clusters_var_pre_omicron) == 0){
      curr_log_lik_var_pre_omicron <- 0.0
    } else{
      curr_log_lik_var_pre_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_var_pre_omicron$n_in_same_clique,
                                                     vec_n_clusters = df_clusters_var_pre_omicron$count,
                                                     R = R_var, k = k_var,
                                                     p = p_trans_before_mut_pre_omicron, p_detect = curr_p_detect,
                                                     max_cluster_size = max_cluster_size_inference)
      
    }
    
    if(nrow(df_clusters_var_post_omicron) == 0){
      curr_log_lik_var_post_omicron <- 0.0
      
    } else{
      curr_log_lik_var_post_omicron <- get_loglik_obs(vec_cluster_size = df_clusters_var_post_omicron$n_in_same_clique,
                                                      vec_n_clusters = df_clusters_var_post_omicron$count,
                                                      R = R_var, k = k_var,
                                                      p = p_trans_before_mut_post_omicron, p_detect = curr_p_detect,
                                                      max_cluster_size = max_cluster_size_inference)
      
    }
    
    return(- (curr_log_lik_non_var_pre_omicron + curr_log_lik_non_var_post_omicron +
                curr_log_lik_var_pre_omicron + curr_log_lik_var_post_omicron))
  }
  
  mle_combined <- optim(par = c(1.0, 0.1),
                        fn = minus_log_lik_combined,
                        method = 'L-BFGS-B',
                        lower = c(R_min, k_min), upper = c(R_max, k_max))
  
  mle_split <- optim(par = c(1.0, 1.0, 0.1),
                     fn = minus_log_lik_split,
                     method = 'L-BFGS-B',
                     lower = c(R_min, R_min, k_min), upper = c(R_max, R_max, k_max))
  
  mle_fully_split <- optim(par = c(1.0, 1.0, 0.1, 0.1),
                           fn = minus_log_lik_fully_split,
                           method = 'L-BFGS-B',
                           lower = c(R_min, R_min, k_min, k_min),
                           upper = c(R_max, R_max, k_max, k_max))
  
  
  
  ## Compute likelihood ratio test p-value for split vs combined
  LRT_split_vs_combined <- as.numeric(2 * (mle_combined$value - mle_split$value))
  p_val_split_vs_combined <- 1.0 - vec_p_quantile[which(vec_quantile_function_chisq > LRT_split_vs_combined)[1]]
  
  ## Compute likelihood ratio test p-value for fully split vs combined
  LRT_fully_split_vs_combined <- as.numeric(2 * (mle_combined$value - mle_fully_split$value))
  p_val_fully_split_vs_combined <- 1.0 - vec_p_quantile[which(vec_quantile_function_chisq_2 > LRT_fully_split_vs_combined)[1]]
  
  ## Compute likelihood ratio test p-value for fully split vs split
  LRT_fully_split_vs_split <- as.numeric(2 * (mle_split$value - mle_fully_split$value))
  p_val_fully_split_vs_split <- 1.0 - vec_p_quantile[which(vec_quantile_function_chisq > LRT_fully_split_vs_split)[1]]
  
  
  
  c('max_log_lik_combined' = - mle_combined$value,
    'max_log_lik_split' = - mle_split$value,
    'max_log_lik_fully_split' = - mle_fully_split$value,
    
    'R_MLE_combined' = mle_combined$par[1],
    'k_MLE_combined' = mle_combined$par[2],
    
    'R_var_MLE_split' = mle_split$par[1],
    'R_non_var_MLE_split' = mle_split$par[2],
    'k_MLE_split' = mle_split$par[3],
    
    'R_var_MLE_fully_split' = mle_fully_split$par[1],
    'R_non_var_MLE_fully_split' = mle_fully_split$par[2],
    'k_var_MLE_fully_split' = mle_fully_split$par[3],
    'k_non_var_MLE_fully_split' = mle_fully_split$par[4],
    
    'LRT_split_vs_combined' = LRT_split_vs_combined,
    'p_val_split_vs_combined' = p_val_split_vs_combined,
    
    'LRT_fully_split_vs_combined' = LRT_fully_split_vs_combined,
    'p_val_fully_split_vs_combined' = p_val_fully_split_vs_combined,
    
    'LRT_fully_split_vs_split' = LRT_fully_split_vs_split,
    'p_val_fully_split_vs_split' = p_val_fully_split_vs_split,
    
    'CV_combined' = mle_combined$convergence,
    'CV_split' = mle_split$convergence,
    'CV_fully_split' = mle_fully_split$convergence,
    
    'p_detect_cases' = vec_p_detect_cases[i_p_detect]
  )
})
stopCluster(cl)

print(df_inference)
