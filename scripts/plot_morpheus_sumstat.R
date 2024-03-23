### Callable Rscrip to transform a logger file which contains 
### time, cell id, x/y position will lateron be extended to 3D.
library(argparse)
suppressPackageStartupMessages(library(tidyverse))

source("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Bachelorarbeit_Vorpraktikum/utils_motility.R")

parser <- ArgumentParser(description = 
                           "Parameters for processing motility data")
parser$add_argument("-st", "--sample_time", default = 30, type = "double",
                    help = "Every which points were data extracted [s].")
parser$add_argument("-te", "--time_end", default = 3600, type = "double",
                    help = "How long was measured [s].")
parser$add_argument("-ts", "--time_shift", default = 240, type = "double",
                    help = "How long is the adaption phase [s].")
parser$add_argument("-vt", "--velocity_th", default = 2, type = "double",
                    help = "Velocity threshold for arrest calculation [µm/min]")
args <- parser$parse_args()


# call setParameters
myParms <- setParameters(samplingStep = args$sample_time, 
                         timeEnd = args$time_end, 
                         timeShift = args$time_shift)

# Retrieve all tracker files
setwd("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Bachelorarbeit_Vorpraktikum")

track_files <- "target_tracker.csv"

mot <- list()
msd <- list()

# Calc motility for each celltype
# for(fn in track_files){
#   celltype <- (str_split(fn, "_") %>% unlist)[1]
  
  # call myReadData
data <- myReadData(file_name = track_files, timeScale = myParms$t_scale) 
  
  # Remove adaption phase by filtering
data <- myFilterTracks(trackData.orig = data, 
                         timeVector = myParms$t_vector)
  
  # Calculate motility data
data <- myMotilityAnalysis(trackData = data, 
                             velocityThreshold = args$velocity_th)
  
  mot[[celltype]] <- data$motility
  msd[[celltype]] <- data$msd

df_mot <- bind_rows(mot, .id = "type") 
df_msd <- bind_rows(msd, .id = "type") 

colnames(df_msd)[2] <- "timepoint"
df_msd$time <- df_msd$timepoint*args$sample_time/60 # In minutes

df_mot %>% gather(inst.speed, straightness, turning.angles, arrest.coeff,
                  key = "Measure", value = "Measurement") -> df_mot

motility <- ggplot(df_mot) + facet_wrap(~Measure, scales = "free_y") + 
              geom_boxplot(mapping = aes(x = type, y = Measurement))

msd <- ggplot(df_msd) + facet_wrap(~type, scales = "free_y") + 
        geom_line(mapping = aes(x = timepoint, y = mean), size = 1.5) +
        geom_ribbon(mapping = aes(x = timepoint, ymin = lower, ymax = upper))

write.csv2(df_mot, file = "tidy_motility.csv")
write.csv2(df_msd, file = "tidy_msd.csv")

ggsave("Motility.png", plot = motility)
ggsave("MSD.png", plot=msd)
