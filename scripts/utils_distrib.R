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



