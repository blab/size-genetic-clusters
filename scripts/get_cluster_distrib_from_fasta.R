library(ape)
library(reshape2)
library(igraph)
library(tidyverse)

source('utils_cluster_alloc.R')

## Load alignment
sequence_data <- read.FASTA('../data/synthetic-fasta.fasta')

## Compute matrix of pairwise distance (measure in number of mutations) from an alignment
dist_mat <- get_dist_mat(sequence_data)

## Compute clusters of identical sequences
cluster_alloc <- get_identical_clusters_from_dist_mat(dist_mat)

## This function returns a list of 2 objects:

## 1/ g_identical: a graph object (from the igraph library)
## where vertices correspond to sequences who have at least another identical sequence in the dataset
## and where vertices are connected if they were indicated as identical in the dist_mat matrix

plot(cluster_alloc$g_identical)

## 2/ df_size_distrib: a dataframe with 2 columns (cluster_size and count)
## describing the size distribution of clusters of identical sequences

plt_cluster_size_distrib(cluster_alloc) # Cluster size distribution

plt_nb_sequences_per_cluster_size(cluster_alloc) # How sequences are spread across cluster sizes
