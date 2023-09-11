library(ape)
library(reshape2)
library(igraph)
library(tidyverse)

source('utils_cluster_alloc.R')

args <- commandArgs(trailingOnly = T)
input_file <- as.character(args[1])
input_file_metadata <- as.character(args[2])
output_file <- as.character(args[3])
save_csv <- as.logical(as.numeric(args[4]))

if(save_csv){
  output_file_csv <- paste0(substr(output_file, start = 1, stop = nchar(output_file) - 4), '.csv')
  output_file_cluster_info_csv <- paste0(substr(output_file, start = 1, stop = nchar(output_file) - 4), '_time.csv')
}

## Load alignment
sequence_data <- read.FASTA(input_file)

## Load sequence dates
sequence_metadata <- read.csv(input_file_metadata) %>% 
  mutate(date = as.Date(date))

## Compute matrix of pairwise distance (measure in number of mutations) from an alignment
dist_mat <- get_dist_mat(sequence_data)

## Compute clusters of identical sequences
cluster_alloc <- get_identical_clusters_from_dist_mat_with_dates(dist_mat, sequence_metadata)
saveRDS(cluster_alloc, output_file)

if(save_csv){
  write.csv(cluster_alloc$df_size_distrib, output_file_csv)
  write.csv(cluster_alloc$df_cluster_info, output_file_cluster_info_csv)
}

## This function returns a list of 2 objects:

## 1/ g_identical: a graph object (from the igraph library)
## where vertices correspond to sequences who have at least another identical sequence in the dataset
## and where vertices are connected if they were indicated as identical in the dist_mat matrix

#plot(cluster_alloc$g_identical)

## 2/ df_size_distrib: a dataframe with 2 columns (cluster_size and count)
## describing the size distribution of clusters of identical sequences

#plt_cluster_size_distrib(cluster_alloc) # Cluster size distribution

#plt_nb_sequences_per_cluster_size(cluster_alloc) # How sequences are spread across cluster sizes
