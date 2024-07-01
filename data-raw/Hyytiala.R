raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/"
# load(paste0(raw.directory, "Hyde_weather_CASSIA.RData")) or just build the package

FMI_Data <- lapply(paste0(raw.directory, list.files(path = raw.directory, pattern = "CASSIA-csv-")), read.csv, na.string  = "-")
# [c(17, 18, 11, 9, 13, 14, 10, 5, 6, 1, 12, 15, 3, 7, 8, 2, 4, 16)]
# 1 = 2014 RH, 2 = 2020 Precip, 3 = 2018 Precip, 4 = 2022 Precip, 5 = 2012 Precip, 6 = 2014 Precip, 7 = 2018 RH, 8 = 2020 RH,
# 9 = 2005 Precip, 10 = 2012 RH, 11 = 2005 RH, 12 = 2016 Precip, 13 = 1010 Precip, 14 = 2010 RH, 15 = 2016 RH, 16 = 2022 RH

FMI_df <- dplyr::bind_rows(FMI_Data)
FMI_df$Date <- as.POSIXct(paste(paste(FMI_df$Year, FMI_df$m, FMI_df$d, sep = "-"), FMI_df$Time), format = "%Y-%m-%d %H:%M", tz = "BST")
FMI_df <- FMI_df[order(FMI_df$Date),]
FMI_df$YMD <- as.Date(FMI_df$Date)
names(FMI_df)[6:7] <- c("RH", "Rain")

FMI_daily <- aggregate(cbind(Rain, RH) ~ YMD, data = FMI_df, mean, na.rm = T, na.action = NULL)
FMI_daily <- FMI_daily[-nrow(FMI_daily),] # Delete the first day of 2023

### NOTE CO2168 is downloaded in the Hyytiala_Data_Creation function
CO2_average <- aggregate(HYY_META.CO2168 ~ lubridate::yday(Date), data = CO2168, mean, na.rm = T, na.action = NULL)$HYY_META.CO2168

data_format$Rain[is.na(data_format$Rain)] <- FMI_daily$Rain[is.na(data_format$Rain)] - mean(FMI_daily$Rain, na.rm = T) + mean(data_format$Rain, na.rm = T)
data_format$RH[is.na(data_format$RH)] <- FMI_daily$RH[is.na(data_format$RH)] - mean(FMI_daily$RH, na.rm = T) + mean(data_format$RH, na.rm = T)
data_format$VPD[is.na(data_format$VPD)] <- bigleaf::rH.to.VPD(0.01*FMI_daily$RH[is.na(data_format$VPD)], data_format$T[is.na(data_format$VPD)])
data_format$CO2[is.na(data_format$CO2)] <- rep(CO2_average, length.out = nrow(data_format))[is.na(data_format$CO2)]

data_IIASA_smear$VPD[is.na(data_IIASA_smear$VPD)] <- bigleaf::rH.to.VPD(0.01*FMI_daily$RH[is.na(data_IIASA_smear$VPD)], data_IIASA_smear$T168[is.na(data_IIASA_smear$VPD)])

## Plot
par(mfrow = c(2, 2))
plot(data_format$Date, data_format$T, main = "Temperature", xlab = "Date", ylab = "Temperature")
plot(data_format$Date, data_format$TSA, main = "Temperature, Soil Depth A", xlab = "Date", ylab = "Temperature")
plot(data_format$Date, data_format$TSB, main = "Temperature, Soil Depth B", xlab = "Date", ylab = "Temperature")
plot(data_format$Date, data_format$MB, main = "Moisture", xlab = "Date", ylab = "Moisture")

plot(data_format$Date, data_format$Rain, main = "Rain", xlab = "Date", ylab = "Rain")
plot(data_format$Date, data_format$PAR, main = "PAR", xlab = "Date", ylab = "PAR")
plot(data_format$Date, data_format$VPD, main = "VPD", xlab = "Date", ylab = "VPD")
plot(data_format$Date, data_format$RH, main = "RH", xlab = "Date", ylab = "RH")

par(mfrow = c(1, 2))
plot(data_format$Date, data_format$CO2, main = "CO2", xlab = "Date", ylab = "CO2")
plot(data_format$Date, data_format$fAPAR, main = "fAPAR", xlab = "Date", ylab = "fAPAR")

# Save file
save(data_format, file = paste0(raw.directory, "Hyde_weather_CASSIA.RData"))

# Test
install.packages("/home/joanna/Asiakirjat/Rprebasso-master/", repos = NULL, type="source")

library(CASSIA)
library(Rprebasso)

hello = PRELES(PAR = data_format$PAR,
               TAir = data_format$T,
               VPD = data_format$VPD,
               Precip = data_format$Rain,
               CO2 = data_format$CO2,
               fAPAR = data_format$fAPAR,
               p = c(100, 0.250, 0.07, 2, rep(NA, 26)),
               returncols = c("GPP", "ET", "SW", "fD", "fW", "fS"))

weather_original_2015 = read.csv(file = "./data/weather_original_2015.csv", header = T, sep = ",")
weather_original_2016 = read.csv(file = "./data/weather_original_2016.csv", header = T, sep = ",")
weather_original_2017 = read.csv(file = "./data/weather_original_2017.csv", header = T, sep = ",")
weather_original = rbind(rbind(weather_original_2015, weather_original_2016), weather_original_2017)

extras = data.frame(Nitrogen = rep(0.012, length = nrow(weather_original)),
                    PAR = data_format[substring(data_format$Date, 1, 4) %in% 2015:2017,c("PAR")],
                    VPD = data_format[substring(data_format$Date, 1, 4) %in% 2015:2017,c("VPD")],
                    CO2 = data_format[substring(data_format$Date, 1, 4) %in% 2015:2017,c("CO2")],
                    fAPAR = rep(0.7, length = nrow(weather_original)))
weather_original <- cbind(weather_original, extras)

save(weather_original, file = "./data/weather_original.RData")

par(mfrow = c(3, 4))
plot(data_format$T[substring(data_format$Date, 1, 4) %in% c(2015:2017)], weather_original$T, ylab = "From Pauliiina", xlab = "My Data", main = "Temperature")
abline(0, 1, col = "red")
plot(data_format$T[substring(data_format$Date, 1, 4) %in% c(2015:2017)] - weather_original$T, ylab = "From Pauliiina", xlab = "My Data")

plot(data_format$Rain[substring(data_format$Date, 1, 4) %in% c(2015:2017)], weather_original$Rain, ylab = "From Pauliiina", xlab = "My Data", main = "Precipitation")
abline(0, 1, col = "red")
plot(data_format$Rain[substring(data_format$Date, 1, 4) %in% c(2015:2017)] - weather_original$Rain, ylab = "From Pauliiina", xlab = "My Data")

plot(data_format$PAR[substring(data_format$Date, 1, 4) %in% c(2015:2017)], weather_original$PAR, ylab = "From Pauliiina", xlab = "My Data", main = "PAR")
abline(0, 1, col = "red")
plot(data_format$PAR[substring(data_format$Date, 1, 4) %in% c(2015:2017)] - weather_original$PAR, ylab = "From Pauliiina", xlab = "My Data")
# TODO: is this to to do with the two different heights?

plot(data_format$MB[substring(data_format$Date, 1, 4) %in% c(2015:2017)], weather_original$MB, ylab = "From Pauliiina", xlab = "My Data", main = "MB")
abline(0, 1, col = "red")
plot(data_format$MB[substring(data_format$Date, 1, 4) %in% c(2015:2017)] - weather_original$MB, ylab = "From Pauliiina", xlab = "My Data")

plot(data_format$TSA[substring(data_format$Date, 1, 4) %in% c(2015:2017)], weather_original$TSA, ylab = "From Pauliiina", xlab = "My Data", main = "TSA")
abline(0, 1, col = "red")
plot(data_format$TSA[substring(data_format$Date, 1, 4) %in% c(2015:2017)] - weather_original$TSA, ylab = "From Pauliiina", xlab = "My Data")

plot(data_format$TSB[substring(data_format$Date, 1, 4) %in% c(2015:2017)], weather_original$TSB, ylab = "From Pauliiina", xlab = "My Data", main = "TSA")
abline(0, 1, col = "red")
plot(data_format$TSB[substring(data_format$Date, 1, 4) %in% c(2015:2017)] - weather_original$TSB, ylab = "From Pauliiina", xlab = "My Data")
