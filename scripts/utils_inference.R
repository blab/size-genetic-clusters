get_proba_size_cluster <- function(cluster_size, R, k, p){
  ## Probability for a cluster to be of size cluster_size 
  ## for a given value of the reproduction number R,
  ## the dispersion parameter k, and the probability p
  ## that transmission occurs before mutation
  
  proba <-
    gamma(k * cluster_size + cluster_size - 1) /
    ( gamma(k * cluster_size) * gamma(cluster_size + 1) ) *
    (p * R / k) ^ (cluster_size - 1) * (1 + p * R / k) ^ (1 - k * cluster_size - cluster_size) 
  
  return(proba)
}

get_log_proba_size_cluster <- function(cluster_size, R, k, p){
  
  ## Log-probability for a cluster to be of size cluster_size 
  ## for a given value of the reproduction number R,
  ## the dispersion parameter k, and the probability p
  ## that transmission occurs before mutation
  
  log_proba <- 
    lgamma(k * cluster_size + cluster_size - 1) - lgamma(k * cluster_size) - lgamma(cluster_size + 1) +
    (cluster_size - 1)*log(p * R / k) - (k * cluster_size + cluster_size - 1) * log(1 + p * R / k)
  
  return(log_proba)
}

## Vectorized version of get_proba_size_cluster (takes as an input a vector of clusters sizes)
get_vec_proba_cluster_size <- function(vec_cluster_size, R, k, p){
  
  vec_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
    get_proba_size_cluster(cluster_size, R, k, p)
  })
  return(vec_proba)
}

## Vectorized version of get_log_proba_size_cluster (takes as an input a vector of clusters sizes)
get_vec_log_proba_cluster_size <- function(vec_cluster_size, R, k, p){
  vec_log_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
    get_log_proba_size_cluster(cluster_size, R, k, p)
  })
  return(vec_log_proba)
}

## Same as above but accounting for partial sequencing of infections
get_vec_proba_obs_cluster_size <- function(vec_cluster_size, R, k, p, p_detect, max_cluster_size){
  ## Probability to observe a cluster of size cluster_size 
  ## for a given value of the reproduction number R,
  ## the dispersion parameter k, and the probability p
  ## that transmission occurs before mutation
  ## and the proportion of infection sequenced (p_detect)
  ## max_cluster_size corresponds to the cluster size threshold c_max
  ## used in the computation of the probabilities
  
  vec_all_clust_size <- 1:max_cluster_size # l
  vec_proba_clust_size <- exp(get_vec_log_proba_cluster_size(vec_all_clust_size, R, k, p)) # r_l
  
  vec_proba_obs_cluster_size <- sapply(vec_cluster_size, FUN = function(j){
    sum(sapply(j:max_cluster_size, FUN = function(l){
      if(p_detect < 1.0) {
        log_proba <- log(vec_proba_clust_size[l]) + lchoose(l, j) + j*log(p_detect) + (l-j)*log(1.0 - p_detect)
        proba <- exp(log_proba)
      }
      else{
        proba <- ifelse(j == l, vec_proba_clust_size[j], 0.)
      }
      
      #vec_proba_clust_size[l] * choose(l, j) * (p_detect ^ j) * ((1. - p_detect) ^ (l - j))
    }))
  })
  proba_obs_zero_element <- sum(sapply(1:max_cluster_size, FUN = function(l){
    vec_proba_clust_size[l] *  (1. - p_detect) ^ l
  }))
  
  return(vec_proba_obs_cluster_size/(1. - proba_obs_zero_element))
}

## Log-likelihood function
get_loglik_obs <- function(vec_cluster_size, vec_n_clusters, R, k, p, p_detect, max_cluster_size){
  ## Log-likelihood of observing the data D = (vec_cluster_size, vec_n_clusters)
  ## where vec_cluster size is a vector of cluster sizes
  ## and vec_n_clusters the number of clusters of each size in vec_cluster_size
  ## given the reproduction number R, the dispersion parameter k,
  ## the probability that transmission occurs before mutation p
  ## the fraction of infections sequenced p_detect
  ## and the cluster size threshold max_cluster_size (c_max)
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size(vec_cluster_size, R, k, p, p_detect, max_cluster_size))
  log_proba_obs <- ifelse(is.infinite(log_proba_obs), -800, log_proba_obs)
  
  loglik <- sum(
    vec_n_clusters * log_proba_obs
  )
  return(loglik)
}

## Same thing as above but without making the sum
get_loglik_obs_no_sum <- function(vec_cluster_size, R, k, p, p_detect, max_cluster_size){
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size(vec_cluster_size, R, k, p, p_detect, max_cluster_size))
  log_proba_obs <- ifelse(is.infinite(log_proba_obs), -800, log_proba_obs)
  
  return(log_proba_obs)
}


## Function running a grid search and computing for each combination of R and k 
## assuming a value for the fraction of infections sequences (p_detect),
## a value for the probability p_trans_before_mut that transmission occurs before mutation
## a cluster size threshold max_cluster_size (c_max)
## and given a cluster size distribution cluster_alloc

run_grid_search <- function(vec_R, vec_k, p_detect, max_cluster_size_inference,
                            cluster_alloc, p_trans_before_mut){
  df_grid <- expand.grid(R = vec_R, k = vec_k) %>% 
    mutate(p_seq_infection = p_detect)
  
  df_grid$log_lik <- sapply(1:nrow(df_grid), FUN = function(i){
    get_loglik_obs(cluster_alloc$df_dist_clique_size$cluster_size,
                   cluster_alloc$df_dist_clique_size$count,
                   df_grid[i, 'R'], df_grid[i, 'k'],
                   p = p_trans_before_mut, p_detect = p_detect,
                   max_cluster_size = max_cluster_size_inference)
  })
  
  return(df_grid)
  
}

## Function running a grid search and computing for each combination of R and k 
## assuming a value for the fraction of infections sequences (p_detect),
## a value for the probability p_trans_before_mut that transmission occurs before mutation
## a cluster size threshold max_cluster_size (c_max)
## and given a cluster size distribution (vec_clusters_sizes, vec_counts_cluster_sizes)

run_grid_search_from_vector <- function(vec_R, vec_k, p_detect,
                                        max_cluster_size_inference,
                                        vec_clusters_sizes, vec_counts_cluster_sizes,
                                        p_trans_before_mut){
  df_grid <- expand.grid(R = vec_R, k = vec_k) %>% 
    mutate(p_seq_infection = p_detect)
  
  df_grid$log_lik <- sapply(1:nrow(df_grid), FUN = function(i){
    get_loglik_obs(vec_clusters_sizes,
                   vec_counts_cluster_sizes,
                   df_grid[i, 'R'], df_grid[i, 'k'],
                   p = p_trans_before_mut, p_detect = p_detect,
                   max_cluster_size = max_cluster_size_inference)
  })
  
  return(df_grid)
  
}

## Function to obtain MLE estimate and likelihood profile confidence intervals
## from a grid search dataframe which has columns with names R, k and col_loglik

get_profile_likelihood_CI <- function(df_grid, col_loglik = 'log_lik_full'){
  # https://web.stat.tamu.edu/~suhasini/teaching613/chapter3.pdf
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)
  
  df_profile_R <- df_grid %>%
    rename(loglik = col_loglik) %>%
    group_by(R) %>% 
    summarise(loglik_profile = max(loglik)) %>% 
    mutate(max_loglik = max(loglik_profile)) %>% 
    mutate(LR = 2 * (max_loglik - loglik_profile),
           is_within_95 = (LR < chisq_stat_95),
           is_within_90 = (LR < chisq_stat_90),
           is_within_50 = (LR < chisq_stat_50))
  
  df_profile_k <- df_grid %>%
    rename(loglik = col_loglik) %>%
    group_by(k) %>% 
    summarise(loglik_profile = max(loglik)) %>% 
    mutate(max_loglik = max(loglik_profile)) %>% 
    mutate(LR = 2 * (max_loglik - loglik_profile),
           is_within_95 = (LR < chisq_stat_95),
           is_within_90 = (LR < chisq_stat_90),
           is_within_50 = (LR < chisq_stat_50))
  
  CI_R <- df_profile_R %>% 
    summarise(mle_estim = R[loglik_profile == max_loglik],
              lower_95 = min(R[is_within_95]),
              upper_95 = max(R[is_within_95]),
              lower_90 = min(R[is_within_90]),
              upper_90 = max(R[is_within_90]),
              lower_50 = min(R[is_within_50]),
              upper_50 = max(R[is_within_50])) %>% 
    mutate(param = 'R')
  
  CI_k <- df_profile_k %>% 
    summarise(mle_estim = k[loglik_profile == max_loglik],
              lower_95 = min(k[is_within_95]),
              upper_95 = max(k[is_within_95]),
              lower_90 = min(k[is_within_90]),
              upper_90 = max(k[is_within_90]),
              lower_50 = min(k[is_within_50]),
              upper_50 = max(k[is_within_50])) %>% 
    mutate(param = 'k')
  
  df_CI <- bind_rows(CI_R, CI_k)
  return(df_CI)
}

## Function to obtain MLE estimate and likelihood profile confidence intervals
## from a grid search dataframe which has columns with names names_param (for the different parameters)
## and col_loglik for the log-likelihood

get_profile_likelihood_CI_multiple_params <- function(df_grid, 
                                                      names_param = c('R', 'k'),
                                                      col_loglik = 'loglik'){
  
  # https://web.stat.tamu.edu/~suhasini/teaching613/chapter3.pdf
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)
  
  df_all_CI <- Reduce('bind_rows', lapply(1:length(names_param), FUN = function(i_param){
    df_profile <- df_grid %>% 
      rename(loglik = col_loglik) %>% 
      rename(curr_param = names_param[i_param]) %>% 
      group_by(curr_param, p_detect_cases) %>% 
      summarise(loglik_profile = max(loglik)) %>% 
      group_by(p_detect_cases) %>% 
      mutate(max_loglik = max(loglik_profile)) %>% 
      mutate(LR = 2 * (max_loglik - loglik_profile),
             is_within_95 = (LR < chisq_stat_95),
             is_within_90 = (LR < chisq_stat_90),
             is_within_50 = (LR < chisq_stat_50))
    
    df_CI <- df_profile %>% 
      group_by(p_detect_cases) %>% 
      summarise(mle_estim = curr_param[loglik_profile == max_loglik],
                lower_95 = min(curr_param[is_within_95]),
                upper_95 = max(curr_param[is_within_95]),
                lower_90 = min(curr_param[is_within_90]),
                upper_90 = max(curr_param[is_within_90]),
                lower_50 = min(curr_param[is_within_50]),
                upper_50 = max(curr_param[is_within_50])) %>% 
      mutate(param = names_param[i_param])
    
    df_CI
  }))
  
  return(df_all_CI)
}

