# Computes the probability p that transmission occurs before mutation
# for a given value of the ratio d_mut / d_gen where:
# d_mut is the mean delay before the appearance of a mutation
# d_gen is the mean generation time
get_proba_trans_before_mut <- function(ratio_delay_between_mutation_delay_generation){
  return(
    ratio_delay_between_mutation_delay_generation / (ratio_delay_between_mutation_delay_generation + 1)
    )
}

## Draw n_draws samples from a negative binomial distribution of dispersion parameter k and mean mean
draw_from_negbin <- function(n_draws, mean, k){
  return(rnbinom(n = n_draws, size = k, mu = mean))
}
## Draw one sample from a binomial distribution of parameter n, p
draw_from_binomial <- function(n, p){
  return(rbinom(n = 1, size = n, prob = p))
}


## Function to compute the theoretical distribution for the size of clusters
theoretical_distrib_clusters <- function(R, k, p, max_cluster_size){
  ## Input:
  ## R: the reproduction number
  ## k: the dispersion parameter
  ## p: the probability that transmission occurs before mutation
  ## max_cluster_size: the maximum size of cluster for which the distribution is computed
  
  vec_j <- 1:max_cluster_size
  vec_rj <- sapply(vec_j, FUN = function(j){
    exp(lgamma(k*j + j - 1) - lgamma(k*j) - lgamma(j + 1)) *
      (p*R/k)^(j -1) / exp((k*j + j  - 1)* log(1+p*R/k) )
  })
  
  
  return(tibble(cluster_size = vec_j,
                prob = vec_rj,
                cum_prob = cumsum(vec_rj)))
}

## Function to simulate clusters of identical sequences for a specific set of parameters
## The simulation of a cluster stops when no more offsprings are generated or when its size reaches max_cluster_size

simulate_branching_process_until_max_cluster_size_save_largest_gen <- function(n_sim, R, k, p, max_cluster_size){
  # Input: 
  # n_sim: Number of outbreaks simulated
  # R: Reproduction number
  # k: Overdispersion parameter of the offspring distribution
  # p: Probability that transmission occurs before mutation
  # max_cluster_size: maximum cluster size
  
  df_sim <- Reduce('bind_rows', lapply(1:n_sim, FUN = function(i_sim){
    # Initialisation
    tot_cluster_size <- 1
    largest_gen_within_cluster <- 0
    contribution_largest_transmission_event <- 0
    largest_gen <- 0
    curr_pop_size <- 1
    i_gen <- 0
    
    # Run 
    while(curr_pop_size > 0 && tot_cluster_size < max_cluster_size){
      n_offsprings <- sum(draw_from_negbin(n_draws = curr_pop_size, mean = R, k = k))
      n_offspring_within_cluster <- draw_from_binomial(n = n_offsprings, p = p)
      largest_gen_within_cluster <- max(largest_gen_within_cluster, n_offspring_within_cluster)
      if(n_offsprings > largest_gen){
        largest_gen <- n_offsprings
        contribution_largest_transmission_event <- n_offspring_within_cluster
      }
      
      curr_pop_size <- n_offspring_within_cluster
      tot_cluster_size <- tot_cluster_size + curr_pop_size
      i_gen <- i_gen + 1
    }
    # Format output
    c('i_sim' = i_sim,
      'gen_extinction' = i_gen,
      'largest_gen' = largest_gen_within_cluster,
      'contrib_largest_transmission_event' = contribution_largest_transmission_event,
      'tot_cluster_size' = tot_cluster_size, 
      'curr_pop_size' = curr_pop_size)
  }))
  
  # Output 
  return(df_sim)
}


