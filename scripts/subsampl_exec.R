df_track_test <- df_track %>% filter(col == "Dense Collagen", type == "Mock")
df_track_frac <- calc_frac(df_track_test)

sub_sample_tracks_500 <- sub_sample_tracks_500 %>% mutate(time = time * 60)


sub_sample_tracks_503 <- sample_tracks3(df_syn_500_2, mock_frac, time_vec) 
sub_sample_tracks_503_sort <- sub_sample_tracks_503 %>% arrange(cell.id)
sub_sample_tracks_503_frac <- calc_frac(sub_sample_tracks_503)

sub_sample_tracks_500_frac <- calc_frac(sub_sample_tracks_500)

sub_sample_tracks_502 <- sample_tracks2(df_track_test, df_syn_500_2)

sub_sample_tracks_500_frac <- calc_frac(sub_sample_tracks_502)

sub_sample_tracks_502 <- sub_sample_tracks_502 %>% mutate(time = time * 60) 
#df_test_2 <- myReadData("/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Bachelorarbeit_Vorpraktikum/df_test.csv")


sub_sample_tracks_503 <- sample_tracks2(df_track_test, df_syn_3_2) %>% arrange(cell.id)
sub_sample_tracks_503 <- sub_sample_tracks_503 %>% mutate(time = time * 60)
sub_sample_tracks_503_frac <- calc_frac(sub_sample_tracks_503)



df_feather_500 <- load_tidy_feather("C:/Users/jodap/OneDrive/Dokumente/Vopraktikum_BA_21/Benjamin_Frey/Bachelorarbeit_Vorpraktikum/test_model_3_500_m.feather", generation = "last")

write.csv(df_track_frac, "C:\\Users\\jodap\\OneDrive\\Dokumente\\Vopraktikum_BA_21\\Benjamin_Frey\\df_frack_mock.csv")
