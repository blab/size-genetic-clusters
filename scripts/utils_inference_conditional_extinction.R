get_proba_size_cluster <- function(cluster_size, R, k, p){
  
  proba <-
    gamma(k * cluster_size + cluster_size - 1) /
    ( gamma(k * cluster_size) * gamma(cluster_size + 1) ) *
    (p * R / k) ^ (cluster_size - 1) * (1 + p * R / k) ^ (1 - k * cluster_size - cluster_size) 
  
  return(proba)
}

get_log_proba_size_cluster <- function(cluster_size, R, k, p){
  
  log_proba <- 
    lgamma(k * cluster_size + cluster_size - 1) - lgamma(k * cluster_size) - lgamma(cluster_size + 1) +
    (cluster_size - 1)*log(p * R / k) - (k * cluster_size + cluster_size - 1) * log(1 + p * R / k)
  
  return(log_proba)
}

get_vec_proba_cluster_size <- function(vec_cluster_size, R, k, p){
  vec_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
    get_proba_size_cluster(cluster_size, R, k, p)
  })
  return(vec_proba)
}

get_vec_log_proba_cluster_size <- function(vec_cluster_size, R, k, p){
  vec_log_proba <- sapply(vec_cluster_size, FUN = function(cluster_size){
    get_log_proba_size_cluster(cluster_size, R, k, p)
  })
  return(vec_log_proba)
}

get_loglik <- function(vec_cluster_size, vec_n_clusters, R, k, p){
  loglik <- sum(
    vec_n_clusters * get_vec_log_proba_cluster_size(vec_cluster_size, R, k, p)
  )
  return(loglik)
}


######## Same as above but accounting for imperfect observation
get_vec_proba_obs_cluster_size <- function(vec_cluster_size, R, k, p, p_detect, max_cluster_size){
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
get_vec_proba_obs_cluster_size_supercritical <- function(vec_cluster_size,
                                                         R, k, p, p_detect,
                                                         max_cluster_size,
                                                         alpha){
  subcritical_R <- R * exp(- alpha * (1 + 1/k))
  
  vec_all_clust_size <- 1:max_cluster_size # l
  vec_proba_clust_size <- exp(get_vec_log_proba_cluster_size(vec_all_clust_size, subcritical_R * p, k, p)) # r_l
  
  vec_proba_obs_cluster_size <- sapply(vec_cluster_size, FUN = function(j){
    sum(sapply(j:max_cluster_size, FUN = function(l){
      if(p_detect < 1.0) {
        log_proba <- - alpha + log(vec_proba_clust_size[l]) + lchoose(l, j) + j*log(p_detect) + (l-j)*log(1.0 - p_detect)
        proba <- exp(log_proba)
      }
      else{
        proba <- ifelse(j == l, vec_proba_clust_size[j], 0.)
      }
      
      #vec_proba_clust_size[l] * choose(l, j) * (p_detect ^ j) * ((1. - p_detect) ^ (l - j))
    }))
  })
  proba_obs_zero_element <- sum(sapply(1:max_cluster_size, FUN = function(l){
    exp(-alpha) * vec_proba_clust_size[l] *  (1. - p_detect) ^ l
  }))
  
  return(vec_proba_obs_cluster_size/(1. - proba_obs_zero_element))
}

estimate_minus_log_proba_extinction <- function(R, p, k){
  step <- 1e-5
  vec_alpha <- seq(step, 50, step)
  
  vec_func <- (
    p * R * (1. - exp(-vec_alpha)) -
    k * (exp(vec_alpha / k) - 1.)
  )
  
  estim_alpha <- mean(vec_alpha[max(which(vec_func > 0))],
                      vec_alpha[min(which(vec_func < 0))], na.rm = T
                      )
  return(estim_alpha)
}
estimate_minus_log_proba_extinction_vec <- function(vec_R, vec_p, vec_k){
  vec_estim <- sapply(1:length(vec_R), FUN = function(i){
    estimate_minus_log_proba_extinction(vec_R[i], vec_p[i], vec_k[i])
  })
  return(vec_estim)
}

get_loglik_obs <- function(vec_cluster_size, vec_n_clusters, R, k, p, p_detect, max_cluster_size){
  if(R * p > 1.0){
    minus_log_proba_extinction <- estimate_minus_log_proba_extinction(R, p, k)
    log_proba_obs <- log(get_vec_proba_obs_cluster_size_supercritical(vec_cluster_size,
                                                                      R, k, p, p_detect, max_cluster_size,
                                                                      minus_log_proba_extinction))
  } else{
    log_proba_obs <- log(get_vec_proba_obs_cluster_size(vec_cluster_size, R, k, p, p_detect, max_cluster_size))
    
  }
  log_proba_obs <- ifelse(is.infinite(log_proba_obs), -800, log_proba_obs)
  
  loglik <- sum(
    vec_n_clusters * log_proba_obs
  )
  
  return(loglik)
}

get_loglik_obs_no_sum <- function(vec_cluster_size, R, k, p, p_detect, max_cluster_size){
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size(vec_cluster_size, R, k, p, p_detect, max_cluster_size))
  log_proba_obs <- ifelse(is.infinite(log_proba_obs), -800, log_proba_obs)
  
  return(log_proba_obs)
}

get_loglik_obs_size_dep <- function(vec_cluster_size, vec_n_clusters, R, k, p, p_sent, max_cluster_size){
  
  log_proba_obs <- log(get_vec_proba_obs_cluster_size_dependent(vec_cluster_size, R, k, p, p_sent, max_cluster_size))
  
  loglik <- sum(
    vec_n_clusters * log_proba_obs
  )
  return(loglik)
}


######### Obtaining likelihood profile CI from the results of a grid search
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



get_CI_from_vec_param <- function(i_param,
                                  minus_log_lik, mle_estim,
                                  vec_param){
  
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)
  
  
  mle_param <- mle_estim$par
  
  ## Log-likelihood profile
  vec_log_lik <- sapply(vec_param, FUN = function(curr_param){
    mle_param[i_param] <- curr_param
    - minus_log_lik(mle_param)
  })
  
  ## Likelihood ratio
  vec_LR <- 2 * (- mle_estim$value - vec_log_lik)
  
  param_within_95 <- vec_param[vec_LR < chisq_stat_95]
  param_within_90 <- vec_param[vec_LR < chisq_stat_90]
  param_within_50 <- vec_param[vec_LR < chisq_stat_50]
  
  return(
    c('lower_95' = min(param_within_95),
      'upper_95' = max(param_within_95),
      'lower_90' = min(param_within_90),
      'upper_90' = max(param_within_90),
      'lower_50' = min(param_within_50),
      'upper_50' = max(param_within_50))
  )
  
  
}

