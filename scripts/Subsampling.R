## ---------------------------
##
## Script name: Subsampling synthetic data
##
## Purpose of script: Try out subsampling on synthetic data
##
## Author: Benjamin Frey
##
## Date Created: 2021-09-04
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
library(tidyverse)
library(dplyr)
library(MotilityLab)

df_syn <- read.table("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Bachelorarbeit_Vorpraktikum/target_tracker.csv", header = T)

# create minute column 
df_syn <- df_syn %>% mutate(time_step = time / 60) %>% filter(time_step%%1==0)

"***********************************************************************"
"Data preprocessing"
"***********************************************************************"

#create list from dataframe
as.tracks.data.frame <- function (x, id.column = 1, time.column = 2, pos.columns = c(3, 4, 5), scale.t = 1, scale.pos = 1, ...) {
  if (ncol(x) < length(pos.columns) + 1) {
    stop("Data frame does not contain enough columns! (Perhaps you need to specify 'sep')")
  }
  if (length(pos.columns) < 1) {
    stop("At least one position column needs to be specified!")
  }
  if (length(pos.columns) == 2 && !is.finite(pos.columns[2])) {
    cx <- match(pos.columns[1], colnames(x))
    if (is.na(cx) && is.numeric(pos.columns[1])) {
      cx <- pos.columns[1]
    }
    pos.columns <- seq(cx, length(x))
  }
  cx <- as.character(c(id.column, time.column, pos.columns))
  cxc <- match(cx, colnames(x))
  cxi <- match(cx, seq_len(ncol(x)))
  cxi[is.na(cxi)] <- cxc[is.na(cxi)]
  if (any(is.na(cxi))) {
    stop("Column(s) not found: ", paste(cx[is.na(cxi)], collapse = ","))
  }
  r <- x[, as.integer(cxi)]
  if (ncol(r) <= 5) {
    colnames(r) <- c("id", "t", c("x", "y", "z")[seq_along(pos.columns)])
  }
  else {
    colnames(r) <- c("id", "t", paste0("x", seq_along(pos.columns)))
  }
  if (scale.t != 1) {
    r[, "t"] <- scale.t * r[, "t"]
  }
  if (any(scale.pos != 1)) {
    r[, -c(1, 2)] <- scale.t * r[, -c(1, 2)]
  }
  sort.tracks(as.tracks.list(split.data.frame(as.matrix(r[, -1]), r[, 1])))
}

#create dataframe as input for subsampling
df_sim <- tracksData.filtered %>% lapply(function(x){as.data.frame(x)}) %>% bind_rows(.id = 'track') %>% as_tibble() %>%
  dplyr::rename(time = t, X = x, Y = y)


#create dataframe from list from tracks 

motility_wrapper <- function(celltype,  # Celltype(str)
                             directory, # directory of the logger file 
                             timeVector,    # passed by setParameters
                             timeThreshold, # adaption phase
                             timeScale,     # passed by setParameters
                             sim = TRUE){
  log_list = list(target = "target_tracker.csv", 
                  infected = "infected_tracker.csv")
  file_ <- log_list[celltype][[1]]
  in.file <- paste(directory, file_, sep = "")
  
  tracksData <- myReadData(file_name = in.file,   # Read data
                           timeScale = timeScale) 
  tracksData.filtered <- myFilterTracks(trackData.orig = tracksData, # Filter
                                        timeVector = timeVector,     # adaption
                                        timeThreshold = timeThreshold)  # phase
  # sampling for cell drop out
  df_sim <-  tracksData.filtered %>% lapply(function(x){as.data.frame(x)}) %>% bind_rows(.id = 'track') %>% as_tibble() %>%
    dplyr::rename(time = t, X = x, Y = y) %>%
    mutate(file = 'sim', col = 'loose', type = 'dnef', rep = '1',time=time-2, track = as.numeric(track)) %>%
    mutate_at(c('file','col','type','rep'), as.factor) %>%
    filter(time <= nrow(df_msd))
  df_samples <- suppressWarnings(lapply(vector(mode = 'list',length = 150),function(x){sample_tracks(df_mot, df_sim)}))
  df_distances <- suppressWarnings(lapply(df_samples, function(x){distance_mot_msd(x, df_mot, df_msd)}) %>% bind_rows(.id = "column_label"))
  best_sample <- sampling_best(df_distances) 
  tracksData.sampled <- df_samples[[best_sample[1]]] %>% select(-file, -col, -type, -rep) %>%
    dplyr::rename(t=time, x=X, y=Y) %>% arrange(t) %>% 
    as.tracks.data.frame(id.column = 'track', time.column ='t', pos.columns = c(3,4))
  
  tracksData.analyzed <- myMotilityAnalysis(trackData = tracksData.sampled)
  # Split output
  msdAnalysis         <- data.frame(tracksData.analyzed[[1]])
  otherMotilityParams <- data.frame(tracksData.analyzed[[2]])
  
  # Fuse to named sumStat
  sumStat        <- list(msdAnalysis, otherMotilityParams)
  names(sumStat) <- return_names(celltype, sim)
  return(sumStat)
}

"************************"
"Subsampling functions"
"************************"

"*****Weitere Variable adden mit Zeitpunkt wo bestimmter %-Satz abfällt.!****"


#give dataset and percentage at the end and function automatically decreases exponentially

subsample2 <- function(df_sync, perc) {
  time_steps <- df_syn %>% select(time) %>% distinct() %>% pull()
  frac <- ((perc / 100)) ** (1/length(time_steps))
  for (i in seq_along(time_steps)) {
    frac_t <- 1 - frac ** i
    drop_tracks <- df_syn %>% filter(time == time_steps[i]) %>% sample_frac(frac_t) %>% select(cell.id) %>% pull()
    if (time_steps[i] == 0) {
      df_syn <- df_syn %>% filter(!cell.id %in% drop_tracks)
      next
    }
    else{
      #sample_tracks <- df_syn %>% filter(cell.id %in% drop_tracks & time < time_steps[i])
      df_drop <- df_syn %>% filter(cell.id %in% drop_tracks, time == time_steps[i])
      df_syn <- anti_join(df_syn, df_drop)
      #df_syn <- bind_rows(df_syn %>% filter(!cell.id %in% drop_tracks), sample_tracks)
    }
  }
  return (df_syn)
}

#latest version
subsample4 <- function(df_syn, perc_vec, time_vec) {
  time_steps <- df_syn %>% select(time_step) %>% unique() %>% pull()
  for (i in seq_along(time_steps)) {
    drop_n <- (1 - perc_vec[i]) 
    drop_tracks <- df_syn %>% filter(time_step == time_vec[i]) %>% sample_frac(drop_n) %>% select(cell.id) %>% pull()
    print(drop_tracks)
    if (time_vec[i] == 0) {
      df_syn <- df_syn %>% filter(!cell.id %in% drop_tracks)
      next
    }
      #   #sample_tracks <- df_syn %>% filter(cell.id %in% drop_tracks & time < time_steps[i])
      #df_drop <- df_syn %>% filter(cell.id %in% drop_tracks, time_step == time_vec[i])
    sampled_track <- df_syn %>% filter(cell.id %in% drop_tracks & time_step < time_vec[i])
    df_syn <- bind_rows(df_syn%>% filter(!cell.id %in% drop_tracks), sampled_track)
      #print(df_syn)
    }
  return (df_syn)
}


subsample5 <- function(df_syn, perc_vec, time_vec) {
  for (i in seq_along(time_vec)) {
    drop_n <- (1 - perc_vec[i]) * 64
    drop_n <- drop_n %>% round(digits = 0)
    print(drop_n)
    drop_tracks <- df_syn %>% filter(time_step == time_vec[i]) %>% sample_n(drop_n) %>% select(cell.id) %>% pull()
    print(drop_tracks)
    if (time_vec[i] == 0) {
       df_syn <- df_syn %>% filter(!cell.id %in% drop_tracks)
       next
    }
    #   #sample_tracks <- df_syn %>% filter(cell.id %in% drop_tracks & time < time_steps[i])
      #df_drop <- df_syn %>% filter(cell.id %in% drop_tracks, time_step == time_vec[i])
    sampled_track <- df_syn %>% filter(cell.id %in% drop_tracks & time_step < time_vec[i])
    df_syn <- bind_rows((df_syn %>% filter(!cell.id %in% drop_tracks)), sampled_track)
    print(head(df_syn))
  }  
  return (df_syn)
}


sampled_track2 <- df_syn %>% filter(cell.id %in% drop_tracks2 & time_step < time_vec[10])
drop_tracks2 <- df_syn %>% filter(time_step == time_vec[10]) %>% sample_frac(1 - mock_frac[10]) %>% select(cell.id) %>% pull()
df_lol <- bind_rows(df_syn %>% filter(!cell.id %in% drop_tracks2), sampled_track2)

#generates random track dropouts at each time step in the magnitude of deviation of data 2 from simulated data (df_1)
#potentially has to be done for n times
sample_tracks2 <- function(df_exp, df_sim){
  df_frac_exp <- calc_frac(df_exp)
  time_steps <- df_frac_exp %>% select(time) %>% pull()
  for (i in seq_along(time_steps)) {
    df_frac_sim <- calc_frac(df_sim)
    n_exp <- df_frac_exp %>% filter(time == time_steps[i]) %>% select(n) %>% pull() + 3
    n_sim <- df_frac_sim %>% filter(time == time_steps[i]) %>% select(n) %>% pull()
    n_dif <- n_sim - n_exp 
      drop_track <- df_sim %>% filter(time == time_steps[i]) %>%
        distinct(track) %>% sample_n(n_dif) %>% select(track) %>% pull()
      if (time_steps[i] == 0) {
        df_sim <- df_sim %>% filter(!track %in% drop_track)
        next
      }
      sampled_track <- df_sim %>% filter(track %in% drop_track & time < time_steps[i])
      df_sim <- bind_rows(df_sim %>% filter(!track %in% drop_track), sampled_track)
    }
  return(df_sim)
}

"*********************************************************************************"
sample_tracks3 <- function(df_sim, frac_vec, time_vec){
  df_frac_sim <- calc_frac(df_sim)
  #time_steps <- df_frac_sim %>% select(time) %>% pull()
  n_vec <- frac_vec * 64 
  n_vec <- n_vec %>% round(digits = 0)
  for (i in seq_along(time_vec)) {
    df_frac_sim <- calc_frac(df_sim)
    n_frac <- n_vec[i]
    n_sim <- df_frac_sim %>% filter(time == time_vec[i]) %>% select(n) %>% pull()
    n_dif <- n_sim - n_frac 
      drop_track <- df_sim %>% filter(time == time_vec[i]) %>% distinct(track) %>% sample_n(n_dif) %>% select(track) %>% pull()
      if (time_vec[i] == 0) {
        df_sim <- df_sim %>% filter(!track %in% drop_track)
        next
    }
    sampled_track <- df_sim %>% filter(track %in% drop_track & time < time_vec[i])
    df_sim <- bind_rows(df_sim %>% filter(!track %in% drop_track), sampled_track)
  }
  return(df_sim)
}


"*********************************************************************************"

n_vec <- mock_frac *64
n_vec <- n_vec  %>% round(digits = 0)
# after calculating distances 

sampling_best <- function(df_distances) {
  df_distances %>% group_by(measure) %>%
    mutate(distance = distance/max(distance)) %>% ungroup() %>% 
    group_by(column_label) %>%
    mutate(distance_agg = mean(distance)) %>% ungroup() %>% distinct(column_label,distance_agg) %>%
    filter(distance_agg == min(distance_agg)) %>% pull(column_label) %>% as.numeric()
}

distance_mot_msd <- function(df_sim, df_exp_mot, df_exp_msd){
  bind_rows(
    distance_mot(df_sim,df_exp_mot) %>% select(measure,distance),
    distance_msd(df_sim, df_exp_msd)
  )
}

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


#loop over sample_tracks n time, inputs simulated data as well as number of iterations, outputs

df_samples <- suppressWarnings(lapply(vector(mode = 'list',length = 150),function(x){sample_tracks3(df_syn_500_2, mock_frac, time_vec)}))
df_distances <- suppressWarnings(lapply(df_samples, function(x){distance_mot_msd(x, df_mot, df_msd)}) %>% bind_rows(.id = "column_label"))
best_sample <- sampling_best(df_distances) 

