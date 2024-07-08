direct = "~/Documents/CASSIA_Calibration/Raw_Data/ancillary_ecosystem/"

####
# Karake
####

karike_list = lapply(paste0(direct, list.files(path = direct)[3:29]), readxl::read_excel, sheet = "karike")
karike_list = lapply(karike_list, function(x) {if (sum(c("pvm.", "pvm....2") %in% names(x)) != 0) {
  names(x)["pvm." == names(x) | "pvm....2" == names(x)] = "pvm"
  }
  return(x)})

karike_df_raw_1 = Reduce(function(x, y) merge(x, y, all=TRUE), karike_list[1:6])
karike_df_raw_1 = karike_df_raw_1[order(karike_df_raw_1$pvm),]
karike_df_raw_2 = Reduce(function(x, y) merge(x, y, all=TRUE), karike_list[7:12])
karike_df_raw_2 = karike_df_raw_2[order(karike_df_raw_2$pvm),]
karike_df_raw_3 = Reduce(function(x, y) merge(x, y, all=TRUE), karike_list[13:18])
karike_df_raw_3 = karike_df_raw_3[order(karike_df_raw_3$pvm),]
karike_df_raw_4 = Reduce(function(x, y) merge(x, y, all=TRUE), karike_list[19:24])
karike_df_raw_4 = karike_df_raw_4[order(karike_df_raw_4$pvm),]
karike_df_raw_5 = Reduce(function(x, y) merge(x, y, all=TRUE), karike_list[25:27])
karike_df_raw_5 = karike_df_raw_5[order(karike_df_raw_5$pvm),]
karike_df_list = list(karike_df_raw_1, karike_df_raw_2, karike_df_raw_3, karike_df_raw_4, karike_df_raw_5)

karike_df_all = dplyr::bind_rows(dplyr::bind_rows(dplyr::bind_rows(dplyr::bind_rows(karike_df_raw_1, karike_df_raw_2), karike_df_raw_3), karike_df_raw_4), karike_df_raw_5)

plot(karike_df_all$pvm, karike_df_all$oksa)
plot(karike_df_all$pvm, karike_df_all$karike)
plot(karike_df_all$pvm, karike_df_all$neulanen)
plot(karike_df_all$pvm, karike_df_all$kuori)
plot(karike_df_all$pvm, karike_df_all$tikku)
plot(karike_df_all$pvm, karike_df_all$käpy)
plot(karike_df_all$pvm, karike_df_all$lehti)

save(karike_df_all, file = "./data/karike_df_all.RData")
