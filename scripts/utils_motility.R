#### This script contains functions that perform motility analysis of 
#### tracked cells (results coming out form simulations or real data)
#### 1. mean square displacement (MSD: stacked based) 
#### 2. instant velocity 
#### 3. straightness
#### 4. turning angles
#### 5. arrest coefficient
#### It uses build-in funcitons from the package Motility Lab

# These functions were given by Frederik from Jana/Paola/Christopher
# There were some numerical issues if the time vectors didnt align perfectly
# Therefore here I add rounding to the set parameters function as well as 
# myFilterTracks function to round timestamps to 5 digits.

#load all necessary packages/libraries
library(MotilityLab)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape))

##### necessary functions -----

setParameters<-function(simulation=TRUE,   
                        samplingStep = 60, # "measure" every 30 seconds 
                        timeEnd = 7200,    # Simulation end
                        timeShift = 240){  # Adaption phase 
  # Returns scale and a vector including all timepoints
  # Timepoints are rounded to 5 digits bc of num instability
  if(simulation)
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


myReadData<-function(file_name, timeScale = 1, dim = 2){
  # Uses the read.tracks.csv function from motility lab
  # timeScale should be inherited from setParameters
  assertthat::assert_that(dim == 2 | dim == 3, 
                          msg = "Dimension has to be 2 or 3")
  if(dim == 2){
    DF <- read.tracks.csv(file_name,  
                          sep = "\t", 
                          id.column = "cell.id", 
                          time.column = "time", 
                          pos.columns = c(3,4), # c(3,NA) would take all 
                                                # colums from 3rd to last
                          header = TRUE, 
                          scale.t = timeScale)}
  if(dim == 3){
    DF <- read.tracks.csv(file_name,  
                          sep = "\t", 
                          id.column = "cell.id", 
                          time.column = "time", 
                          pos.columns = c(2,3,4), # c(3,NA) would take all 
                          # colums from 3rd to last
                          header = TRUE, 
                          scale.t = timeScale)
  }
  
  return(DF)
}


myFilterTracks<-function(trackData.orig, 
                         timeVector = 0, 
                         timeThreshold = 0){
  # Filter timeThreshold data points (is to be given in the correct timescale)
  # or filter for all the given timeVector points
  trackData <- trackData.orig
  
  # 1. take out the first timeThreshold of simulation 
  #(time thought for the sim to reach steady-state)
  if(timeThreshold>0){
    trackData<-as.tracks(lapply(trackData, 
                                function(x) x[x[,1] >= timeThreshold, ]))  
  }
  
  # 2. extract the relevant data points according to the timeVector
  if (length(timeVector) > 1){
    trackData <- as.tracks(lapply(trackData, 
                                function(x) x[round(x[,1], 5) %in% timeVector, ]))
    # Assert if we lost some timepoints due to numerics
    assertthat::assert_that((trackData[[3]] %>% dim)[1] == (length(timeVector)), 
                             msg = "Filtered data lost some data along the way")
  }
  
  return(trackData)
}


myMotilityAnalysis<-function(trackData, velocityThreshold = 2){
  # returns a list with two data frames: the first one containing the MSD analysis 
  # and the second one containing speed, straightness, turning angles, and arrest coefficient
  
  # start Analysis
#  print("calculating the MSD")
  msd <- aggregate(trackData, squareDisplacement, FUN = "mean.se")
  
#  print("calculating the instant velocity for all cells")
  inst.velocity <- as.data.frame(sapply(trackData,speed))
  
#  print("calculating straightness for all cells")
  straight.ness <- as.data.frame(sapply(trackData, straightness))
  
#  print("calculating turning angles")
  turning.angles <- as.data.frame(sapply(trackData,meanTurningAngle))
  
#  print("calculating arrest coefficient")
  arrest <- sapply(trackData, function(x){
    tmp  <- sapply(subtracks(x,1),speed)
    
    if(length(tmp) > 0){  
      arrest <- length(tmp[tmp < velocityThreshold])/length(tmp)
    }  else{
      arrest <- NA
    }
  })
  
  out_else <- data.frame(inst.velocity,straight.ness, turning.angles, arrest)
  colnames(out_else) <- c("inst.speed", 
                          "straightness", 
                          "turning.angles", 
                          "arrest.coeff")
  
  return(list("msd"      = msd,
              "motility" = out_else))
}
# 
# addParamTypeToDF <- function(DF, param, type){
#   # Add column Parameter value to data frame
#   DF <- data.frame(DF, rep(param, nrow(DF)),  
#                   rep(type, nrow(DF)))
#   
#   colnames(DF)[c(5,6)] <- c("param_value", "type")
#   # Add Param_value to data.frame
#   param_value  <- as.character(DF[,5])
#   param_value2 <- c()
#   for (i in param_value){
#     tmp          <- strsplit(i, "_")
#     param_value2 <- c(param_value2, tmp[[1]][length(tmp[[1]])])
#   }
#   DF[,5] <- as.character(param_value2)
#   return(DF)
# }


check_end_slash <- function(path){
  # Checks that the path ends with a slash. If not a "/" is appended
  #Check that directory ends with "/"
  if(strsplit(path,"")[[1]][length(strsplit(path,"")[[1]])] == "/"){
    path = path
  } else{
    path = paste0(path, "/")
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
  if(sim == TRUE){
    mot <- paste0("sumMot", celltype, "Sim")
    msd <- paste0("sumMSD", celltype, "Sim")
  } else {
    mot <- paste0("sumMot", celltype, "Exp")
    msd <- paste0("sumMSD", celltype, "Exp")
  }
  return(c(msd, mot))
}

motility_wrapper <- function(celltype, 
                             directory, 
                             timeVector, 
                             timeThreshold, 
                             timeScale,
                             sim = TRUE){
  # Wrapper function for motility analysis
  # Celltype(str)
  # directory = where is the logger file 
  # timeVector = vector of all time points supplied by setParameters
  # timeThreshold = adaption phase
  file_ <- return_logger(celltype = celltype)
  
  tmp     <- strsplit(directory,"/")
  #param   <- tmp[[1]][length(tmp[[1]])]       # Last element
  in.file <- paste(directory, file_, sep = "")
  
  tracksData <- myReadData(file_name = in.file,   # Read data
                           timeScale = timeScale) 
  tracksData.filtered <- myFilterTracks(trackData.orig = tracksData, # Filter
                                        timeVector = timeVector,     # adaption
                                        timeThreshold = timeThreshold)  # phase
  tracksData.analyzed <- myMotilityAnalysis(trackData = tracksData.filtered)
  
  msdAnalysis <- data.frame(tracksData.analyzed[[1]])
  otherMotilityParams <- data.frame(tracksData.analyzed[[2]])
  sumStat <- list(msdAnalysis, otherMotilityParams)
  names(sumStat) <- return_names(celltype, sim)
  return(sumStat)
}
