
## Get distance matrix from a list of sequences
get_dist_mat <- function(sequence_data){
  
  if(length(unique(sapply(sequence_data, length))) > 1){
    stop('Sequences should be aligned and have the same length!')
  }
  
  dist_mat <- dist.dna(sequence_data, model = 'N',
                       as.matrix = T, pairwise.deletion = T)
  
  return(dist_mat)
}

## Allocating clusters of identical sequences based on a distance matrix
get_identical_clusters_from_dist_mat <- function(dist_mat){
  
  ## This function returns a list of 2 objects:
  
  ## 1/ g_identical: a graph object (from the igraph library)
  ## where vertices correspond to sequences who have at least another identical sequence in the dataset
  ## and where vertices are connected if they were indicated as identical in the dist_mat matrix
  
  ## 2/ df_size_distrib: a dataframe with 2 columns (cluster_size and count)
  ## describing the size distribution of clusters of identical sequences
  
  ## Number of sequences without any other identical sequences
  n_singletons <- dist_mat %>%  melt() %>%  as_tibble() %>% 
    rename(seq_1 = Var1, seq_2 = Var2, n_mut = value) %>% 
    filter(n_mut == 0) %>% 
    group_by(seq_1) %>% summarise(n_identical = n()) %>% 
    filter(n_identical == 1) %>% 
    nrow()
  
  if(n_singletons == nrow(dist_mat)){
    g_identical <- NULL
    df_size_distrib <- tibble(cluster_size = 1, count = n_singletons)
    
    return(list(
      g_identical = g_identical,
      df_size_distrib = df_size_distrib
    ))
  }
  
  ## Get a dataframe with distinct identical pairs of sequences (edge list)
  df_identical <- dist_mat %>%  melt() %>%  as_tibble() %>% 
    rename(seq_1 = Var1, seq_2 = Var2, n_mut = value) %>% 
    filter(n_mut == 0, as.character(seq_1) > as.character(seq_2))
  
  ## Generate the graph asociated to the edge list 
  g_identical <- graph_from_edgelist(el = as.matrix(df_identical[, c('seq_1', 'seq_2')]), directed = F)
  
  ## Gett the maximal cliques in g_identical
  max_cliques_identical <- max_cliques(g_identical)
  
  ## Reorder the cliques by increasing size
  size_cliques <- unlist(lapply(max_cliques_identical, length))
  max_cliques_identical <- max_cliques_identical[order(size_cliques)]
  
  ## To ensure that all sequences within a cluster are identical, 
  ## clusters are defined by starting with the smallest maximal clique. 
  ## When considering a new maximal clique, the sequences of this clique
  ## who have not yet been allocated to a cluster of identical sequences
  ## are allocated to the same cluster.
  
  cluster_membership <- rep(NA, length(V(g_identical))) # Vector of cluster membership
  
  for(clique_id in 1:length(max_cliques_identical)){
    for(vertex_id in max_cliques_identical[[clique_id]]){
      if(is.na(cluster_membership[vertex_id])){
        cluster_membership[vertex_id] <- clique_id
      }
    }
  }
  
  V(g_identical)$clique_id <- cluster_membership
  
  df_clusters <- tibble(seq = V(g_identical)$name,
                        clique_id = V(g_identical)$clique_id)
  
  
  ## Get distribution of clique size
  df_size_distrib <- df_clusters %>% 
    group_by(clique_id) %>% 
    summarise(cluster_size = n()) %>% 
    group_by(cluster_size) %>% 
    summarise(count = n())
  
  df_size_distrib <-
    bind_rows(
      tibble(cluster_size = 1, count = n_singletons),
      df_size_distrib
    ) %>% 
    group_by(cluster_size) %>% 
    summarise(count = sum(count))
  
  return(list(
    g_identical = g_identical,
    df_size_distrib = df_size_distrib
  ))
}


## Functions to display the cluster size distribution
plt_cluster_size_distrib <- function(cluster_alloc, color_plot = 'darkslateblue'){
  plt_cluster_dist <- cluster_alloc$df_size_distrib %>% 
    ggplot(aes(x = as.factor(cluster_size), y = count)) +
    geom_bar(stat = 'identity', fill = color_plot) +
    theme_classic() +
    geom_text(aes(label = count), vjust = -0.5) +
    scale_x_discrete(name = 'Size of clusters of\nidentical sequences',) +
    scale_y_continuous(name = 'Number of clusters',
                       expand = expansion(mult = c(0.01, 0.1)))
  
  return(plt_cluster_dist)
}

plt_nb_sequences_per_cluster_size <- function(cluster_alloc, color_plot = 'darkslateblue'){
  
  plt_sequences_dist <- cluster_alloc$df_size_distrib %>% 
    mutate(n_sequences = count*cluster_size) %>% 
    ggplot(aes(x = as.factor(cluster_size), y = n_sequences)) +
    geom_bar(stat = 'identity', fill = color_plot) +
    theme_classic() +
    geom_text(aes(label = n_sequences), vjust = -0.5) +
    scale_x_discrete(name = 'Size of clusters of identical sequences') +
    scale_y_continuous(name = 'Number of sequences',
                       expand = expansion(mult = c(0.01, 0.1)))
  
  return(plt_sequences_dist)
}






