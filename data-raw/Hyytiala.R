load(paste0(raw.directory, "Hyde_weather_CASSIA.RData"))

FMI_Data <- lapply(paste0(raw.directory, list.files(path = raw.directory, pattern = "CASSIA-csv-"))[c(11:12, c(9, 5), c(6, 1), c(13, 10), c(7, 3), c(2, 8), c(14, 4))], read.csv, na.string  = "-")
FMI_df <- dplyr::bind_rows(FMI_Data)
names(FMI_df)[6:7] <- c("Rain", "RH")
FMI_df$Date <- as.POSIXct(paste(paste(FMI_df$Year, FMI_df$m, FMI_df$d, sep = "-"), FMI_df$Time), format = "%Y-%m-%d %H:%M", tz = "BST")
FMI_df <- FMI_df[order(FMI_df$Date),]
FMI_df$YMD <- as.Date(FMI_df$Date)

FMI_daily <- aggregate(cbind(Rain, RH) ~ YMD, data = FMI_df, mean, na.rm = T, na.action = NULL)
FMI_daily <- FMI_daily[-nrow(FMI_daily),]

### NOTE CO2168 is downloaded in the Hyytiala_Data_Creation function
CO2_average <- length(aggregate(HYY_META.CO2168 ~ lubridate::yday(Date), data = CO2168, mean, na.rm = T, na.action = NULL)$HYY_META.CO2168)

data_format$Rain[is.na(data_format$Rain)] <- FMI_daily$Rain[is.na(data_format$Rain)] - mean(FMI_daily$Rain, na.rm = T) + mean(data_format$Rain, na.rm = T)
data_format$RH[is.na(data_format$RH)] <- FMI_daily$RH[is.na(data_format$RH)] - mean(FMI_daily$RH, na.rm = T) + mean(data_format$RH, na.rm = T)
data_format$VPD[is.na(data_format$VPD)] <- bigleaf::rH.to.VPD(0.01*FMI_daily$RH[is.na(data_format$VPD)], data_format$T[is.na(data_format$VPD)])
data_format$CO2 <- zoo::na.approx(rep(CO2_average, nrow(data_format)))

save(data_format, file = paste0(raw.directory, "Hyde_weather_CASSIA.RData"))
