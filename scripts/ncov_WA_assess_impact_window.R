library(tidyverse)
library(stringr)
library(foreach)
library(doParallel)

source('utils_inference.R')
source('utils_cluster_alloc.R')

name_variant <- 'Alpha'
i_day <- 10 # Length of time window of analysis

## Input files
file_input_WA_analysis <- '../data/ncov_WA/input_WA_analysis_review.csv'
file_input_cluster_dist <- paste0('../data/ncov_WA/list_df_cluster_dist_by_period_', name_variant, '.rds')
file_proba_trans_before_mut <- '../results/proba_trans_before_mut/df_p_trans_before_mut_with_uncertainty.rds'

df_scenarios <- read_csv(file_input_WA_analysis)
if(! name_variant %in% df_scenarios$variant){
  print('Need to pick a variant whose name is in the df_scenarios dataframe!')
}

## Time window of analysis
date_beginning_window <- as.Date(unlist(df_scenarios[df_scenarios$variant == name_variant, 'date_reached_10_seq']))
date_first_end_window <- as.Date(unlist(df_scenarios[df_scenarios$variant == name_variant, 'date_reached_10_seq']))
date_end_window <- date_first_end_window + 60
vec_date_end_window <- seq.Date(date_first_end_window, date_end_window, by = 'day')

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

## Datasets for the analysis
list_df_cluster_dist_by_period <- readRDS(file_input_cluster_dist)

vec_Nextstrain_clades_Omicron <- c("21K (Omicron)", "21L (Omicron)", 
                                   "21M (Omicron)", 
                                   "22A (Omicron)", "22B (Omicron)",
                                   "22C (Omicron)", "22D (Omicron)",
                                   "22E (Omicron)", "22F (Omicron)")

## Define function used to run analysis
run_analysis_given_cluster_ditrib <- function(list_df_cluster_distrib, 
                                              p_trans_before_mut_pre_omicron, p_trans_before_mut_post_omicron, 
                                              n_cores = 4){
  ## Scenarios for the proportion of cases detected
  vec_p_detect_cases <- c(0.1, 0.2, 0.5, 0.8)
  curr_prop_cases_sequenced <- list_df_cluster_distrib$n_seq_during_timewindow / list_df_cluster_distrib$n_cases_during_timewindow
  vec_p_detect <- vec_p_detect_cases * curr_prop_cases_sequenced
  
  ## Dataset for the analysis
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
  
  ## Setting up the bounds for the optimization
  R_min <- 0.01
  R_max <- 5.0
  k_min <- 0.001
  k_max <- 10.0
  
  max_cluster_size_inference <- 10000
  
  cl <- makeForkCluster(n_cores)
  registerDoParallel(cl)
  
  
  df_inference <- Reduce('bind_rows', foreach(i_p_detect = 1:length(vec_p_detect_cases), .packages = c('dplyr')) %dopar% {
    
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
    
    mle_combined <- optim(par = c(1.0, 0.1),
                          fn = minus_log_lik_combined,
                          method = 'L-BFGS-B',
                          lower = c(R_min, k_min), upper = c(R_max, k_max))
    
    mle_split <- optim(par = c(1.0, 1.0, 0.1),
                       fn = minus_log_lik_split,
                       method = 'L-BFGS-B',
                       lower = c(R_min, R_min, k_min), upper = c(R_max, R_max, k_max))
    
    ## Compute likelihood ratio test p-value for split vs combined
    LRT_split_vs_combined <- as.numeric(2 * (mle_combined$value - mle_split$value))
    p_val_split_vs_combined <- 1.0 - vec_p_quantile[which(vec_quantile_function_chisq > LRT_split_vs_combined)[1]]
    
    
    c('max_log_lik_combined' = - mle_combined$value,
      'max_log_lik_split' = - mle_split$value,
      
      'R_MLE_combined' = mle_combined$par[1],
      'k_MLE_combined' = mle_combined$par[2],
      
      'R_var_MLE_split' = mle_split$par[1],
      'R_non_var_MLE_split' = mle_split$par[2],
      'k_MLE_split' = mle_split$par[3],
      
      'LRT_split_vs_combined' = LRT_split_vs_combined,
      'p_val_split_vs_combined' = p_val_split_vs_combined,
      
      'CV_combined' = mle_combined$convergence,
      'CV_split' = mle_split$convergence,
      
      'p_detect_cases' = vec_p_detect_cases[i_p_detect]
    )
  })
  
  stopCluster(cl)
  
  return(df_inference)
}

## Setting up the inference
max_cluster_size_inference <- 10000
n_cores <- 4

vec_p_quantile <- c(0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11,
                    1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, seq(2e-5, 1.0, 1e-5))
vec_quantile_function_chisq <- qchisq(vec_p_quantile, df = 1)

curr_day <- vec_date_end_window[i_day]
list_df_cluster_distrib <- list_df_cluster_dist_by_period[[i_day]]

df_inference <- run_analysis_given_cluster_ditrib(list_df_cluster_distrib, 
                                                  p_trans_before_mut_pre_omicron, p_trans_before_mut_post_omicron, 
                                                  n_cores)

print(df_inference)
