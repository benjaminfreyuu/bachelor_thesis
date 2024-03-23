# New abc_r_blueprint.R, where we store all functions as the atm we cant work
# with relative paths for automatic job submission.

# Load Libraries
library(MotilityLab)
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
#library(parallel)

home <- "/home/bq_charmel/dnef/data/"
sumMSDtargetExp <- read.table(paste0(home, "msd_par.csv"), 
                              sep = ",", header = T)
sumMottargetExp <- read.table(paste0(home, "mot_par.csv"),
                              sep = ",", header = T) %>% filter(rep != 'mock')
# sumMSDinfectedExp <- read.table(paste0(home, "loose_MSD_infected.csv"),
#                                 sep=",", header=T)
# sumMotinfectedExp <- read.table(paste0(home, "loose_Mot_infected.csv"),
#                                 sep=",", header=T)


# create dataframes for sampling cell drop out
df_mot <- sumMottargetExp %>%
  mutate(track = sumMottargetExp %>% group_by(track,file) %>% group_indices()) %>%
  filter(col == "loose", type == 'dnef')
df_msd <- sumMSDtargetExp %>% filter(col == "loose", type == 'dnef')




filter_mot <- function(df,collagen,type_exp){
  df %>%
    filter(col == collagen, type == type_exp) %>%
    group_by(col,type) %>%
    select(track,straight, speed, angle, arrest,col,type) %>%
    distinct()
}


filter_msd <- function(df,collagen,type_exp){
 df %>% filter(type == type_exp, col == collagen) 
}


Exp <- list(sumMSDtargetExp = filter_msd(sumMSDtargetExp,'loose','dnef'),
 sumMottargetExp = filter_mot(sumMottargetExp,'loose','dnef'))#, 
           # sumMSDinfectedExp=sumMSDinfectedExp, sumMotinfectedExp=sumMotinfectedExp)

# Set Parameters
celltypes <- c("target")
time_Shift <- 120
time_Step <- 60
time_End   <- nrow(Exp[['sumMSDtargetExp']])*time_Step + time_Shift
use_dist <- "normal" # "ks" or "normal"
dist_to_index <- list("normal" = 1,
                      "ks"     = 2)


setParameters <- function(simulation=TRUE,   
                        samplingStep = 30, # "measure" every 30 seconds 
                        timeEnd = 3720,    # Simulation end
                        timeShift = 120){  # Adaption phase 
  # Returns scale and a vector including all timepoints
  # Timepoints are rounded to 5 digits bc of num instability
  if (simulation)
  {
    time_scale  <- 1/60   # To convert to minutes
    time_vector <- seq(time_scale*timeShift,
                       timeEnd*time_scale,
                       samplingStep*time_scale) %>% round(5)
  } else{
    time_scale  <- 0.5        
    time_vector <- NULL
  }
  return(list("t_scale" = time_scale,
              "t_vector" = time_vector))
}

# sampling functions 
calc_frac <- function(df){
  df %>%
    group_by(time) %>%
    distinct(track) %>% 
    summarise(n = n()) %>%
    mutate(frac = n/max(n))
}
# arrest function 
arrest <- function(x){
  tmp  <- sapply(subtracks(x,1),speed)
  if (length(tmp) > 0) {  
    return(length(tmp[tmp < 2])/length(tmp))
  } else {
    return(NA)
  }
} 
# Define motility statistics to be calculated 
stats_list <- list('straight' = straightness,
                   'speed' = speed,
                   'angle' = meanTurningAngle,
                   'arrest' = arrest,
                   'square_displace' = squareDisplacement,
                   'duration' = duration,
                   'length' = trackLength)
# Function for applying motility lab functions to dataframe 
calc_stats <- function(stats_list,tracks_list,df){
  for (i in c(1:length(stats_list))) {
    stats <- lapply(tracks_list, stats_list[[i]]) %>%
      enframe(name = "track") %>%
      mutate(value = unlist(value),
             track = as.numeric(track))
    names(stats) <- c('track',names(stats_list)[i])
    df <- right_join(df,stats,by = "track")
  }
  return(df)
}


sample_tracks <- function(df_exp, df_sim){
  df_frac_exp <- calc_frac(df_exp)
  time_steps <- df_frac_exp %>% select(time) %>% pull()
  for (i in seq_along(time_steps)) {
    df_frac_sim <- calc_frac(df_sim)
    n_exp <- df_frac_exp %>% filter(time == time_steps[i]) %>% select(n) %>% pull()
    n_sim <- df_frac_sim %>% filter(time == time_steps[i]) %>% select(n) %>% pull()
    if (n_exp >= n_sim) {
      next
    } else {
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
  }
  return(df_sim)
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


myReadData <- function(file_name, timeScale = 1, dim = 2){
  # Uses the read.tracks.csv function from motility lab
  # timeScale should be inherited from setParameters
  assertthat::assert_that(dim == 2 | dim == 3, 
                          msg = "Dimension has to be 2 or 3")
  if (dim == 2) {
    df <- read.tracks.csv(file_name,  
                          sep = "\t", 
                          id.column = "cell.id", 
                          time.column = "time", 
                          pos.columns = c(3,4), # c(3,NA) would take all 
                          # colums from 3rd to last
                          header = TRUE, 
                          scale.t = timeScale)}
  if (dim == 3) {
    df <- read.tracks.csv(file_name,  
                          sep = "\t", 
                          id.column = "cell.id", 
                          time.column = "time", 
                          pos.columns = c(3,4,5), # c(3,NA) would take all 
                          # colums from 3rd to last
                          header = TRUE, 
                          scale.t = timeScale)
  }
  return(df)
}


myFilterTracks <- function(trackData.orig, 
                         timeVector = 0, 
                         timeThreshold = 0){
  # Filter timeThreshold data points (is to be given in the correct timescale)
  # or filter for all the given timeVector points
  trackData <- trackData.orig
  
  # 1. take out the first timeThreshold of simulation 
  #(time thought for the sim to reach steady-state)
  if (timeThreshold > 0) {
    trackData <- as.tracks(lapply(trackData, 
                                function(x) x[x[,1] >= timeThreshold, ]))  
  }
  
  # 2. extract the relevant data points according to the timeVector
  if (length(timeVector) > 1) {
    trackData <- as.tracks(lapply(trackData, 
                                  function(x) x[round(x[,1], 5) %in% timeVector, ]))
    # Assert if we lost some timepoints due to numerics
    assertthat::assert_that((trackData[[1]] %>% dim)[1] == (length(timeVector)), 
                            msg = "Filtered data lost some data along the way...")
  }
  return(trackData)
}


myMotilityAnalysis <- function(trackData, velocityThreshold = 2){
  # returns a list with two data frames: the first one containing the MSD analysis 
  # and the second one containing speed, straightness, turning angles, and arrest coefficient
  msd            <- aggregate(trackData, squareDisplacement, FUN = "mean.se")
  inst.velocity  <- as.data.frame(sapply(trackData,speed))
  straight.ness  <- as.data.frame(sapply(trackData, straightness))
  turning.angles <- as.data.frame(sapply(trackData,meanTurningAngle)) %>%
                        mutate_all(coalesce, 0) # nan to 0

  arrest <- sapply(trackData, function(x){tmp  <- sapply(subtracks(x,1),speed)
    if (length(tmp) > 0) {  
      arrest <- length(tmp[tmp < velocityThreshold])/length(tmp)
    }  else{
      arrest <- NA
    }
  })
  
  out_else <- data.frame(inst.velocity, 
                         straight.ness, 
                         turning.angles, 
                         arrest)
  
  colnames(out_else) <- c("speed", 
                          "straight", 
                          "angle", 
                          "arrest")
  
  return(list("msd"      = msd,
              "motility" = out_else))
}

check_end_slash <- function(path){
  # Checks that the path ends with a slash. If not a "/" is appended
  #Check that directory ends with "/"
  if (strsplit(path,"")[[1]][length(strsplit(path,"")[[1]])] == "/") {
    path <- path
  } else{
    path <- paste0(path, "/")
  }
  return(path)
}


return_logger <- function(celltype){
  # Given the celltype as string it returns a given logger name which is 
  # predefined in this function.
  log_list = list(target = "target_tracker.csv", 
                  infected = "infected_tracker.csv")
  assertthat::assert_that(any(names(log_list) == celltype), 
                          msg = "Celltype has no associated logger file.
                          Please check return_logger function.")
  return(log_list[celltype][[1]])
}


return_names <- function(celltype, sim){
  # Return names for the summaryStatistic List given a celltype(str)
  # If Sim == True, Sim is appended else it's Exp
  if (sim == TRUE) {
    mot <- paste0("sumMot", celltype, "Sim")
    msd <- paste0("sumMSD", celltype, "Sim")
  } else {
    mot <- paste0("sumMot", celltype, "Exp")
    msd <- paste0("sumMSD", celltype, "Exp")
  }
  return(c(msd, mot))
}

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

motility_wrapper <- function(celltype,  # Celltype(str)
                             directory, # directory of the logger file 
                             timeVector,    # passed by setParameters
                             timeThreshold, # adaption phase
                             timeScale,     # passed by setParameters
                             sim = TRUE){
  file_ <- return_logger(celltype = celltype)
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


mySummaryStatistics <- function(directory, 
                                celltype=celltypes, 
                                timeEnd=time_End, 
                                timeShift=time_Shift, 
                                samplingStep=time_Step, 
                                simulation=TRUE){
  #' Calculates summary statistics from cell position logger (so far 2D)
  #' 
  #' It is important that the summary statistics have names to store it correctly 
  #' in pyABC's database
  #' @param celltype vector or char that has all the celltypes to be analysed
  #' @param directory Directory where to find logger files
  #' @param timeEnd Until which timepoint to evaluate
  #' @param timeShift Length of adaption phase not evaluated
  #' @param samplingStep Interval length of sampling.
  #' @simulation Data from simulation
  #' @return Named list of summary statistics
  
  Sim <- list()
  
  # timeScale (in hours): if the data is coming from the experiments then each 
  # time point is 0.5; 
  # if the data is from the simulations then 1/60 (assuming each MCS is one sec) 
  # timeVector: vector with the time points measrued in the simulations, 
  # i.e. one measurement every 30 s
  time <- setParameters(samplingStep = samplingStep, 
                        timeEnd      = timeEnd, 
                        timeShift    = timeShift,
                        simulation   = simulation)
  timeScale  <- time[[1]] 
  timeVector <- time[[2]]
  
  adaption_phase <- timeShift*timeScale
  directory <- check_end_slash(directory[['loc']])
  
  # Create sumStat for all celltypes
  for (t in celltype) {
    output <- motility_wrapper(celltype = t, 
                               directory = directory, 
                               timeVector = timeVector, 
                               timeThreshold = adaption_phase,
                               timeScale = timeScale)
    
    Sim <- append(Sim, output)
  }
  return(Sim)
}


distance_normal <- function(Sim, Exp){
  # Calculates a distance measure between experimental data and simulation,
  # given the mean + the sigma of experimental/simulated motility distribution.
  # To use this distance measure, we assume the distributions to be normal.
  mean_exp <- Exp %>% mean()
  mean_sim <- Sim %>% mean()
  sig_exp <- Exp %>% sd()
  sig_sim <- Sim %>% sd()
  
  dist <- (mean_exp - mean_sim)**2/sig_exp**2 + 
    ((sig_exp/mean_exp) - (sig_sim/mean_exp))**2
  return(dist)
}


distance_kolmogorov <- function(Sim, Exp){
  dist <- ks.test(Sim, Exp)[["statistic"]] %>% unname()
  return(dist)
}

distance_list <- list("normal" = distance_normal,
                      "ks"     = distance_kolmogorov)

dist_fun <- distance_list[[dist_to_index[[use_dist]]]]

distance_speed_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["speed"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["speed"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                   Exp = sumMotExp)
  
	return(dist)
}

distance_speed_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["speed"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["speed"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  return(dist)
}


distance_str_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["straight"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["straight"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)

  return(dist)
}

distance_str_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["straight"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["straight"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  
  return(dist)
}


distance_trn_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]]
  sumMotSim <- Sim[[names_sim[[2]]]]
  
  # Cells that are not moving at all == nan!
  nan_count <- sum(is.nan(sumMotSim[["angle"]]))
  sumMotSim_trn <- sumMotSim %>%
    select(angle) %>%
    filter(!is.nan(angle))
  
  if (dim(sumMotSim_trn)[1] > 0) {
    dist <- dist_fn(Sim = sumMotSim_trn[["angle"]],
                    Exp = sumMotExp[["angle"]])
  } else {
    dist <- 0  # Initialize dist
  }
  count_NaN <- sumMotSim %>%
    select(angle) %>%
    filter(is.nan(angle)) %>%
    count() %>% pull()
  fraction_NaN <- count_NaN/dim(sumMotSim)[1]
  dist <- dist + fraction_NaN*2*pi  # Fraction of NaN add maximal distance
  
  return(dist)
}

distance_trn_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]]
  sumMotSim <- Sim[[names_sim[[2]]]]
  
  # Cells that are not moving at all == nan!
  nan_count <- sum(is.nan(sumMotSim[["angle"]]))
  sumMotSim_trn <- sumMotSim %>%
    select(angle) %>%
    filter(!is.nan(angle))
  
  if (dim(sumMotSim_trn)[1] > 0) {
    dist <- dist_fn(Sim = sumMotSim_trn[["angle"]],
                    Exp = sumMotExp[["angle"]])
  } else {
    dist <- 0  # Initialize dist
  }
  count_NaN <- sumMotSim %>%
    select(angle) %>%
    filter(is.nan(angle)) %>%
    count() %>% pull()
  fraction_NaN <- count_NaN/dim(sumMotSim)[1]
  dist <- dist + fraction_NaN*2*pi  # Fraction of NaN add maximal distance
  
  return(dist)
}



distance_arr_tar <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["arrest"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["arrest"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  
  return(dist)
}


distance_arr_inf <- function(Sim, Exp, dist_fn = dist_fun, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMotExp <- Exp[[names_exp[[2]]]][["arrest"]]
  sumMotSim <- Sim[[names_sim[[2]]]][["arrest"]]
  
  dist <- dist_fn(Sim = sumMotSim,
                  Exp = sumMotExp)
  return(dist)
}


distance_msd_tar <- function(Sim, Exp, celltype = c("target")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMSDExp <- Exp[[names_exp[[1]]]]
  sumMSDSim <- Sim[[names_sim[[1]]]]
  
  mn_e <- sumMSDExp[["mean"]]
  mn_s <- sumMSDSim[["mean"]]
  sd_e <- sumMSDExp[["mean"]] - sumMSDExp[["lower"]]
  sd_s <- sumMSDSim[["mean"]] - sumMSDSim[["lower"]]
  
  dist <- ((mn_e - mn_s)/sd_e)**2 %>% mean()
  
  return(dist)
}

distance_msd_inf <- function(Sim, Exp, celltype = c("infected")){
  
  names_exp <- return_names(celltype = celltype, sim = FALSE) # 1st MSD 2nd Mot
  names_sim <- return_names(celltype = celltype, sim = TRUE)  # 1st MSD 2nd Mot
  sumMSDExp <- Exp[[names_exp[[1]]]]
  sumMSDSim <- Sim[[names_sim[[1]]]]
  
  mn_e <- sumMSDExp[["mean"]]
  mn_s <- sumMSDSim[["mean"]]
  sd_e <- sumMSDExp[["mean"]] - sumMSDExp[["lower"]]
  sd_s <- sumMSDSim[["mean"]] - sumMSDSim[["lower"]]
  
  dist <- ((mn_e - mn_s)/sd_e)**2 %>% mean()
  
  return(dist)
}






