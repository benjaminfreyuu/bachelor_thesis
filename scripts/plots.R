"***********************************************************************"
"generate combined plot for msd with different parameters and dimensions"
"***********************************************************************"

df_250_1 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_6/tidy_msd.csv")
df_250_2 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_11/tidy_msd.csv")
df_250_3 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_10/tidy_msd.csv")

df_350_1 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_13/tidy_msd.csv")
df_350_2 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_14/tidy_msd.csv")
df_350_3 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_15/tidy_msd.csv")

df_500_1 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_9/tidy_msd.csv")
df_500_2 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_8/tidy_msd.csv")
df_500_3 <- read.csv2("C:/Users/jodap/morpheus/Test_Model_1_16/tidy_msd.csv")

df_msd_all <- bind_rows("35_15_250_1" = df_250_1, "35_15_250_2" = df_250_2, "35_15_250_3" = df_250_3, "35_15_350_1" = df_350_1, "35_15_350_2" = df_350_2, "35_15_350_3" = df_350_3,"35_15_500_1" = df_500_1,"35_15_500_2" = df_500_2,"35_15_500_3" = df_500_3, .id = "parms")
df_msd_all <- df_msd_all %>% separate(parms, c("PS", "DT", "dim", "replicate"), sep = "_")

msd <- df_msd_all %>% 
  ggplot(aes(x = timepoint, y = mean/1000, col = replicate)) + 
  geom_ribbon(aes(x = timepoint, ymin = lower/1000, ymax = upper/1000, fill = replicate), alpha = 0.3) +
  facet_wrap(~dim, scales = "free_y") + 
  geom_line(size = 1) +
  theme_bw(base_size = 22) + 
  ylab("MSD [µm² x 10³]") + xlab("Time [min]") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        legend.title = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.text = element_text(size=15))
print(msd)

"**************************************"
"compare subsampled and normal motility"
"**************************************"

subsample <- read.csv("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Data_files/29.4/Subsampled/target_mot_sub.csv")
subsample <- subsample %>% select(-X)
subsample$type <- plyr::mapvalues(subsample$type, from = "1", to = "target")
mot <- read.csv2("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Data_files/29.4/tidy_motility.csv", sep = ";")
mot <- mot %>% select(-X)
mot$ID <- rep(1:64, 4) 
mot <- mot %>% spread(key = "Measure", value = "Measurement") 
mot <- mot %>% rename(arrest = arrest.coeff, straight = straightness, angles = turning.angles, speed = inst.speed)



combine_mot2 <- function(exp, sim){
  df <- bind_rows('Sub' = exp, 'Sim' = sim, .id = 'data')
  return(df)
}

combine_mot2(subsample, mot)

compare_mot2 <- function(exp, sim){
  df_mot <- combine_mot2(exp = exp,
                         sim = sim)
  df_mot <- df_mot %>% gather(key = 'sumstat', value = 'value', arrest:angles)
  plot <- ggplot(df_mot) +
    facet_wrap(~sumstat, scale = 'free_y', 
               strip.position = "left", 
               labeller = as_labeller(c(speed = "Mean velocity [µm/min]", 
                                        angles = "Mean turning angle [rad]",
                                        arrest="Arrest coefficient", 
                                        straight="Straightness"))) + 
    stat_boxplot(mapping = aes(x = data, y = value, col = data), show.legend = T) +
    geom_point(mapping = aes(x = data, y = value, col = data), position=position_jitterdodge(.1)) +
    ylab("") + xlab("") +
    theme_bw(base_size = 22) + 
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank(),
          legend.position = "none",
          panel.spacing = unit(2, "lines"))
  return(plot)
}
compare_mot2(subsample, mot)

df_mot <- combine_mot2(exp = subsample,
                       sim = mot)
df_mot <- df_mot %>%  gather(key = 'sumstat', value = 'value', arrest:angles)


"**************************************"
"compare subsampled and normal msd"
"**************************************"
subsample_msd <- read.csv("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Data_files/29.4/Subsampled/target_msd_sub.csv")
subsample_msd <- subsample_msd %>% rename(time = i)
msd <- read.csv2("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Data_files/29.4/tidy_msd.csv")

combine_msd2 <- function(exp, sim){
  df <- bind_rows('Sub' = exp, 'Sim' = sim, .id = 'data')
  return(df)
}

compare_msd2 <- function(exp, sim){
  df_msd <- combine_msd2(exp = exp,
                        sim = sim)
  plot <- df_msd %>% 
    ggplot(aes(x=time, y=mean/1000, col = data)) +
    geom_ribbon(aes(ymin=lower/1000,ymax=upper/1000, fill=data),alpha=0.3) +
    facet_wrap(~type, scales = "free_y", ncol = 1) + 
    geom_line() + theme_bw(base_size = 22) + 
    ylab("MSD [µm² x 10³]") + xlab("Time [min]") +
    theme(strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank(),
          panel.spacing = unit(3, "lines"),
          strip.text = element_text(size=22))
  return(plot)
}

compare_msd2(subsample_msd, msd)
df_msd <-combine_msd2(subsample_msd, msd)
