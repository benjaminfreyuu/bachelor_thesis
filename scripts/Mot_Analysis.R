






distance_msd <- function(df_sim, df_exp){
  df_sim <- msd_analysis(df_sim) %>% mutate(id = 'sim')
  df <- bind_rows(df_exp %>% mutate(id = 'exp'), df_sim) %>%
    mutate(sd = mean - lower) %>% select(mean,sd,id,i) %>% gather(measure,value,-id,-i) %>%
    unite(measure,c('measure','id')) %>% spread(measure,value) %>%
    summarise(distance=((mean_exp - mean_sim)/sd_exp)**2 %>% mean(na.rm=TRUE)) %>% mutate(measure='msd')
  return(df)
}

distance_mot <- function(df_sim, df_exp){
  df_sim <- motility_analysis(df_sim) %>% mutate(id = 'sim')
  df <- bind_rows(df_exp %>% mutate(id = 'exp'), df_sim) %>%
    select(straight:arrest, id, track) %>% distinct() %>% select(-track) %>%
    group_by(id) %>%
    summarise_all(list(mean = ~mean(x=. , na.rm = TRUE),
                       sig = ~sd(x=., na.rm = TRUE))) %>% ungroup() %>%
    gather(measure,value,-id)  %>% separate(measure,c('measure','stat'),sep = '_') %>%
    unite(stat,c('stat','id')) %>%
    spread(stat,value) %>% group_by(measure) %>%
    mutate(distance = (mean_exp - mean_sim)**2/sig_exp**2 + 
             ((sig_exp/mean_exp) - (sig_sim/mean_exp))**2) %>% ungroup()
  return(df)
}

distance_mot_msd <- function(df_sim, df_exp_mot, df_exp_msd){
  bind_rows(
    distance_mot(df_sim,df_exp_mot) %>% select(measure,distance),
    distance_msd(df_sim, df_exp_msd)
  )
}
sampling_best <- function(df_distances) {
  df_distances %>% group_by(measure) %>%
    mutate(distance = distance/max(distance)) %>% ungroup() %>% 
    group_by(column_label) %>%
    mutate(distance_agg = mean(distance)) %>% ungroup() %>% distinct(column_label,distance_agg) %>%
    filter(distance_agg == min(distance_agg)) %>% pull(column_label) %>% as.numeric()
}
# Function for applying motility analysis to whole dataframe 
motility_analysis <- function(df,scale_time = 1, scale_dim = 1){
  df <- df %>% mutate(X_scale = X * scale_dim,
                      Y_scale = Y * scale_dim,
                      time_scale = time * scale_time)
  files <- df %>% distinct(file) %>% pull()
  out_df <- data.frame()
  for (j in c(1:length(files))) {
    motile_df <- df %>%
      filter(file == files[j]) %>%
      group_by(track) %>%
      filter(n() > 2) %>% ungroup()
    tracks <- distinct(motile_df, track) %>% pull()
    tracks_list <- vector("list", length(tracks))
    names(tracks_list) <- tracks
    for (i in c(1:length(tracks))) {
      m <- motile_df %>% filter(track == tracks[i]) %>% select(time_scale,X_scale,Y_scale) %>% arrange(time_scale) %>% as.matrix()
      tracks_list[[i]] <-  m
    }
    out_df <- rbind.data.frame(out_df,
                               calc_stats(stats_list, tracks_list, motile_df) %>%
                                 mutate(file = files[j], rmsd = sqrt(square_displace),angle=coalesce(angle,0)))
  }
  return(out_df)
}


msd_analysis <- function(df, scale_time = 1, scale_dim = 1){
  df <- df %>% mutate(X_scale = X * scale_dim,
                      Y_scale = Y * scale_dim,
                      time_scale = time * scale_time)
  files <- df  %>% distinct(col,type)
  collagen <- files %>% select(col) %>% pull()
  type_exp <- files %>% select(type) %>% pull()
  out_df <- tibble()
  for (j in c(1:length(collagen))) {
    motile_df <- df %>%
      filter(col == collagen[j], type == type_exp[j]) %>%
      group_by(track,rep) %>%
      filter(n() > 2) %>% ungroup()
    replicates <- motile_df %>% distinct(rep) %>% pull()
    tracks_list_all <- list()
    for (k in c(1:length(replicates))) {
      motile_df_2 <- motile_df %>% filter(rep == replicates[k])
      tracks <- motile_df_2 %>% distinct(track) %>% pull()
      tracks_list <- vector("list", length(tracks))
      names(tracks_list) <- tracks
      for (i in c(1:length(tracks))) {
        m <- motile_df_2 %>%
          filter(track == tracks[i]) %>%
          select(time_scale,X_scale,Y_scale) %>%
          arrange(time_scale) %>%
          as.matrix()
        tracks_list[[i]] = m
      }
      tracks_list_all <- c(tracks_list_all, tracks_list)
      names(tracks_list_all) <- as.character(seq(1, length(tracks_list_all), 1))
      
    }
    bind_df <- aggregate(as.tracks(tracks_list_all), squareDisplacement, FUN = 'mean.se') %>%
      as_tibble() %>%
      mutate(col = collagen[j], type = type_exp[j])
    out_df <- rbind.data.frame(out_df, bind_df)
  }
  return(out_df)
}