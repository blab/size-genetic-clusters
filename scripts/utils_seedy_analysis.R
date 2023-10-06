get_freq_donor_recepient <- function(W){
  
  mat_deepseq <- deepseqmat(W)
  
  sample_at_recovery <-  W$sampledata %>% 
    as_tibble() %>% 
    arrange(pID) %>% 
    group_by(pID) %>% 
    filter(sampletimes == max(sampletimes)) %>% 
    ungroup()
  
  mat_deepseq_at_recovery <- Reduce('cbind', lapply(1:nrow(sample_at_recovery), FUN = function(id_row){
    mat_deepseq[, sample_at_recovery$sampleID[id_row]]
  }))
  colnames(mat_deepseq_at_recovery) <- paste0('patient ', sample_at_recovery$pID)
  
  df_freq_mutant <- mat_deepseq_at_recovery %>% 
    as_tibble() %>% 
    mutate(locus = substr(rownames(mat_deepseq_at_recovery), start = 7, stop = 1e9),
           locus = ifelse(locus == 'NA', "0", locus),
           locus = as.numeric(locus)) %>% 
    pivot_longer(cols = - 'locus',
                 names_prefix = 'patient ', names_to = 'patient_id',
                 values_to = 'freq_mutant')
  
  df_ids_transm_pairs <- W$epidata %>% 
    as_tibble() %>% 
    filter(inf.source != 0) %>% 
    rename(donor_id = inf.source, recepient_id = ID) %>% 
    select(donor_id, recepient_id) %>% 
    mutate(id_pair = 1:n())
  
  df_SNPs_transm_pair <- Reduce('bind_rows', lapply(1:nrow(df_ids_transm_pairs), FUN = function(id_pair){
    id_recepient <- as.character(df_ids_transm_pairs[id_pair, 'recepient_id'])
    id_donor <- as.character(df_ids_transm_pairs[id_pair, 'donor_id'])
    
    df_freq_donor <- df_freq_mutant %>% 
      filter(patient_id == id_donor) %>% 
      select(-patient_id) %>% 
      rename(freq_donor = freq_mutant)
    
    df_freq_recepient <- df_freq_mutant %>% 
      filter(patient_id == id_recepient) %>% 
      select(-patient_id) %>% 
      rename(freq_recepient = freq_mutant)
    
    df_freq_donor_recepient <- df_freq_donor %>% 
      left_join(df_freq_recepient, by = c('locus')) %>% 
      filter(freq_donor + freq_recepient > 0.) %>% 
      mutate(id_pair = id_pair,
             id_recepient = id_recepient,
             id_donor = id_donor)
  })) %>% 
    filter(freq_donor + freq_recepient > 0.)
  
  return(df_SNPs_transm_pair)
}
