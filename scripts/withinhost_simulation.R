## Install the SEEDY package from CRAN archive is not already installed
is_installed <- ('seedy' %in% .packages(all.available = TRUE))
if(! is_installed){
  seedy_url <- "https://cran.r-project.org/src/contrib/Archive/seedy/seedy_1.3.tar.gz"
  install.packages(seedy_url, repos = NULL, type = 'source')
}
## Load libraries
library(seedy)
library(tidyverse)
library(ggplot2)
library(doParallel)
library(foreach)

source('utils_seedy_analysis.R')

## Create function to define model
create_simulation_model <- function(mean_gen_time,
                                    mut_rate, genome_length, 
                                    bottleneck_size, Ne_equilibrium,
                                    output_freq){
  
  simulate_pairs <- function(n_pairs){
    # Generate timing of infection of the pairs
    vec_time_to_infection <- ceiling(rexp(n_pairs, rate = 1. / mean_gen_time))
    
    vec_infection_times <- rep(0, n_pairs + 1)
    vec_recovery_times <- rep(0, n_pairs + 1)
    
    for(i_pair in 1:n_pairs){
      vec_infection_times[i_pair + 1] <- vec_infection_times[i_pair] + vec_time_to_infection[i_pair]
      vec_recovery_times[i_pair] <- vec_infection_times[i_pair + 1] + 1
    }
    vec_recovery_times[n_pairs + 1] <- max(vec_recovery_times[1:n_pairs] + 1)
    
    
    # Simulate the pairs
    sim <- simfixoutbreak(ID = 1:(n_pairs + 1),
                          inf.times = vec_infection_times,
                          rec.time = vec_recovery_times,
                          inf.source = 0:n_pairs,
                          mut.rate = mut_rate, glen = genome_length,
                          inoc.size = bottleneck_size,
                          imp.var = 25,
                          samp.schedule = "random",
                          samples.per.time = 1,
                          samp.freq = 1,
                          full = TRUE,
                          feedback = output_freq,
                          ref.strain = NULL)
    
    return(sim)
  }
  
  
  return(
    list(
      simulate_pairs = simulate_pairs
    )
  )
}


## Define the characteristics of the scenario to be run
inf_duration <- 6. # Infectious period duration in days 
bottleneck_size <- 1 # Bottleneck size

n_generations_per_day <- 10 # Number of generations per time unit (day)
mean_gen_time <- inf_duration * n_generations_per_day
recovery_rate <- 1. / (inf_duration * n_generations_per_day) # Rate of recovery (per generation)
av_delay_between_mutations <- 12. # Average delay between mutations (in days)
mut_rate <- 1. / (av_delay_between_mutations * n_generations_per_day) # Mutation rate (per generation)
genome_length <- 10000 # Genome length
Ne_equilibrium <- 1000 # Equilibrium effective population size in the host
output_freq <- 500 # Frequency of feedback (verbose)

## Initialize model
model_sim <- create_simulation_model(mean_gen_time = mean_gen_time, 
                                     mut_rate = mut_rate, genome_length = genome_length, 
                                     bottleneck_size = bottleneck_size, Ne_equilibrium = Ne_equilibrium,
                                     output_freq = output_freq)

## Run simulations in parallel
n_pairs_per_sim <- 20
n_rep <- 60 ## (Can reduce to decrease computing time)

n_cores <- 4 # Number of threads used for multi-threading
cl <- makeForkCluster(n_cores)
registerDoParallel(cl)
df_SNPs_transm_pair <- foreach(i_rep = 1:n_rep, .combine = bind_rows, .packages = c('dplyr', 'tidyr', 'seedy')) %dopar% {
  set.seed(i_rep)
  sim <- model_sim$simulate_pairs(n_pairs_per_sim)
  get_freq_donor_recepient(sim) %>% 
    mutate(i_rep = i_rep)
}
stopCluster(cl)
t1 <- Sys.time()
print(t1 - t0)
