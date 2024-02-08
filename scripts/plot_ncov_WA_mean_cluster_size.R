library(tidyverse)

# Load data describing the size and timing of clusters of identical sequences in WA by Nextstrain clade
df_timing <- read_csv('../data/ncov_WA/df_timing_cluster_identical_sequences.csv')

# Load dataframe describing how we want to group Nextstrain clades by larger clades for the plot
df_grouping_clade <- read_csv('../data/ncov_WA/df_grouping_clade_for_plot.csv')
  
# Compute mean cluster size by clade and month of first cluster detection
df_mean_timing <- df_timing %>%
  left_join(df_grouping_clade, by = 'Nextstrain_clade') %>%
  mutate(month_first_date = format(first_date_cluster, '%y/%m')) %>% 
  group_by(month_first_date, new_grouping) %>%
  summarise(n_clusters = n(),
            mean_clust_size = mean(cluster_size)) %>% 
  ungroup()

# Only display mean cluster size that were computed from a least 5 clusters
df_mean_timing_for_plot <- df_mean_timing %>%
  filter(n_clusters >= 5) # Only display months where at least 5 clusters of a given Pango lineage were observed

# The next block of code is used to create a grouping feature (modif_group)
# that will be used to draw lines on the plot. This is not a necessary step
# but makes the plot nicer and avoids lines being drawn between two
# non-consecutive months

df_mean_timing_for_plot <- df_mean_timing_for_plot %>%
  mutate(corr_date = as.Date(paste0(month_first_date, '/01'), format = '%y/%m/%d')) %>%
  group_by(new_grouping) %>%
  mutate(modif_group = case_when(corr_date == min(corr_date) ~ paste0(new_grouping, '_1')),
         cur_id_group = case_when(corr_date == min(corr_date) ~ 1)) %>%
  arrange(new_grouping, corr_date) %>%
  ungroup()

for(i_row in 1:nrow(df_mean_timing_for_plot)){
  if(is.na(df_mean_timing_for_plot$modif_group[i_row])){
    curr_previous_date_in_data <- df_mean_timing_for_plot %>%
      filter(corr_date < df_mean_timing_for_plot$corr_date[i_row],
             new_grouping == df_mean_timing_for_plot$new_grouping[i_row]) %>%
      filter(corr_date == max(corr_date)) %>%
      mutate(corr_date = as.character(corr_date)) %>%
      select(corr_date) %>% unlist()
    
    curr_max_group_id <- df_mean_timing_for_plot %>%
      filter(corr_date < df_mean_timing_for_plot$corr_date[i_row],
             new_grouping == df_mean_timing_for_plot$new_grouping[i_row]) %>%
      filter(corr_date == max(corr_date)) %>%
      select(cur_id_group) %>% unlist()
    
    bool_previous_month_in_data <- (as.numeric(df_mean_timing_for_plot$corr_date[i_row] - as.Date(curr_previous_date_in_data))) <= 31
    
    if(bool_previous_month_in_data){
      df_mean_timing_for_plot$modif_group[i_row] <- paste0(df_mean_timing_for_plot$new_grouping[i_row], '_', curr_max_group_id)
      df_mean_timing_for_plot$cur_id_group[i_row] <- curr_max_group_id
    }  else{
      df_mean_timing_for_plot$modif_group[i_row] <- paste0(df_mean_timing_for_plot$new_grouping[i_row], '_', curr_max_group_id + 1)
      df_mean_timing_for_plot$cur_id_group[i_row] <- curr_max_group_id + 1
    }
  }
}


## Colour palette used for the clades
col_palettes <- c(brewer.pal(12, 'Paired')[1:4],
                  'darkgrey',
                  brewer.pal(9, 'YlOrRd')[-(1:2)],
                  'black')
## Corresponding clades
vec_clade_colors <- c('Alpha', 'Beta', 'Gamma', 'Delta',
                      'Other',
                      'BA.1',
                      'BA.2 (21L)',
                      'BA.2.12.1',
                      'BA.4', 'BA.5',
                      'BA.2.75', 'BQ.1', 'XBB')

# Plot of mean cluster size over time
plt_clust_over_time <- df_mean_timing_for_plot %>%
  ggplot(aes(x = month_first_date,
             y = mean_clust_size,
             group = modif_group,
             colour = new_grouping)) +
  geom_line(alpha = 0.5) +
  geom_point(size = 2, shape = 16) +
  scale_colour_manual(name = 'Variant',
                      values = col_palettes,
                      breaks = vec_clade_colors) +
  scale_x_discrete(name = '',
                   breaks = format(seq.Date(from = as.Date('2020-03-01'), to = as.Date('2023-01-01'), by = 'month'), '%y/%m'),
                   labels = format(seq.Date(from = as.Date('2020-03-01'), to = as.Date('2023-01-01'), by = 'month'), '%b\n%y')
  ) +
  scale_y_continuous(name = 'Mean polytomy size',
                     limits = c(1., NA),
                     expand = expansion(mult = c(0.01, 0.05))) +
  theme_classic() +
  guides(colour = guide_legend(ncol = 3)) +
  theme(legend.position = 'bottom')


plot(plt_clust_over_time)
