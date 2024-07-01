#' @export
Hyytiala_Data_Creation <- function(raw.directory,
                                   year_start,
                                   year_end,
                                   download = F,
                                   clean_data = F,
                                   save = T) {

  raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/"
  year_start = 1995
  year_end = 2023
  download = T
  clean_data = T
  save = T

  ####
  # Download the data from the SMEAR database
  ####

  if (download) {
    http.origin = "https://smear-backend.rahtiapp.fi/search/timeseries/csv?tablevariable=HYY_"
    for (variable in c(paste0("META.", c("RH672", "RH1250", "RHTd", "PAR", "CO2168", "T168", "T336", "Precip", "tsoil_5", "tsoil_10", "wsoil_B1", "wsoil_B2", "Glob", "Glob67", "Pamb336", "wpsoil_A", "wpsoil_B")), "EDDY233.GPP")[16:17]) {
      from = "&from="
      year1 = seq(year_start, year_end, by = 2)
      to = "-01-01T00%3A00%3A00.000&to="
      year2 = seq(year_start+1, year_end, by = 2)
      if (length(year2) != length(year1)) {year2 <- c(year2, year_end)}
      end = "-12-31T23%3A59%3A59.999&quality=ANY&aggregation=NONE&interval=1"

      files = as.list(paste0(http.origin, variable, from, year1, to, year2, end))
      imported_file = as.list(paste0(raw.directory, variable, from, year1, "to", year2))

      mapply(download.file, files, imported_file)
    }
  }

  ####
  # Clean the data with other Hyytiala data
  ####

  if (clean_data) {
    ####
    # Read the data
    #
    # NOTE! Don't run if no new data downloaded!
    ####
    RH672 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RH672")), read.csv, dec = "."))
    RH1250 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RH1250")), read.csv, dec = "."))
    RHTd <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RHTd")), read.csv, dec = "."))
    PAR <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "PAR")), read.csv, dec = "."))
    CO2168 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "CO2168")), read.csv, dec = "."))
    T336 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "T336")), read.csv, dec = "."))
    T168 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "T168")), read.csv, dec = "."))
    Precip <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Precip")), read.csv, dec = "."))
    tsoil_5 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "tsoil_5")), read.csv, dec = "."))
    tsoil_10 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "tsoil_10")), read.csv, dec = "."))
    wsoil_B1 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wsoil_B1")), read.csv, dec = "."))
    wsoil_B2 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wsoil_B2")), read.csv, dec = "."))
    wpsoil_A <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wpsoil_A")), read.csv, dec = "."))
    wpsoil_B <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wpsoil_B")), read.csv, dec = "."))
    GPP <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "GPP")), read.csv, dec = "."))
    GPP$Date <- as.Date(paste(GPP$Year, GPP$Month, GPP$Day, sep = "-"), "%Y-%m-%e")
    Glob <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Glob&")), read.csv, dec = "."))
    Glob$Date <- as.Date(paste(Glob$Year, Glob$Month, Glob$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")
    Glob67 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Glob67")), read.csv, dec = "."))
    Glob67$Date <- as.Date(paste(Glob67$Year, Glob67$Month, Glob67$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")
    CO2168$Date <- as.Date(paste(CO2168$Year, CO2168$Month, CO2168$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")

    RH = merge(RH1250, RHTd, by = c("Year", "Month", "Day", "Hour", "Minute", "Second"), all = T)
    RH$Date = as.POSIXct(paste(paste(RH$Year, RH$Month, RH$Day, sep = "-"), paste(RH$Hour, RH$Minute, sep = ":")),
                         format = "%Y-%m-%d %H:%M", tz = "BST")
    RH <- RH[order(RH$Date),]

    all_data <- as.data.frame(cbind(RH672$HYY_META.RH672, RH$HYY_META.RH1250, RH$HYY_META.RHTd,
                                    PAR$HYY_META.PAR, CO2168$HYY_META.CO2168,
                                    T336$HYY_META.T336, T168$HYY_META.T168,
                                    Precip$HYY_META.Precip,
                                    tsoil_5$HYY_META.tsoil_5, tsoil_10$HYY_META.tsoil_10,
                                    wsoil_B1$HYY_META.wsoil_B1, wsoil_B2$HYY_META.wsoil_B2))
    names(all_data) <- c("RH672", "RH1250", "RHTd", "PAR", "CO2168", "T336", "T168",
                         "Precip", "tsoil_5", "tsoil_10", "wsoil_B1", "wsoil_B2")
    all_data$Date <- as.Date(paste(RH672$Year, RH672$Month, RH672$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")
    all_data$Date_Time <- as.POSIXct(paste(RH672$Year, RH672$Month, RH672$Day, RH672$Hour, RH672$Minute, RH672$Second, sep = "-"), format = "%Y-%m-%d-%H-%M-%S", tz = "BST")

    wpsoil_A$Date_Time <- as.POSIXct(paste(wpsoil_A$Year, wpsoil_A$Month, wpsoil_A$Day, wpsoil_A$Hour, wpsoil_A$Minute, wpsoil_A$Second, sep = "-"), format = "%Y-%m-%d-%H-%M-%S", tz = "BST")
    wpsoil_B$Date_Time <- as.POSIXct(paste(wpsoil_B$Year, wpsoil_B$Month, wpsoil_B$Day, wpsoil_B$Hour, wpsoil_B$Minute, wpsoil_B$Second, sep = "-"), format = "%Y-%m-%d-%H-%M-%S", tz = "BST")

    all_data <- merge(all_data, wpsoil_A, all = T)
    all_data <- merge(all_data, wpsoil_B, all = T)

    ### Save data
    data_storage <- "/home/josimms/Documents/Austria/Plant-FATE/tests/"
    # TODO: maybe make this in the CASSIA package at some point!
    write.csv(all_data[1:(nrow(all_data)/2),], file = paste0(data_storage, "all_data.csv"))
    write.csv(all_data[(nrow(all_data)/2+1):nrow(all_data),], file = paste0(data_storage, "all_data_2.csv"))
    # LAST RUN 01.07.2024
      # Don't need to run if you haven't downloaded new data

    ### Read data
    data_storage <- "/home/josimms/Documents/Austria/Plant-FATE/tests/"
    all_data_1 <- read.csv(paste0(data_storage, "all_data.csv"))
    all_data_2 <- read.csv(paste0(data_storage, "all_data_2.csv"))
    all_data <- rbind(all_data_1, all_data_2)
    names(all_data) <- gsub("HYY_META.", "", names(all_data))

    # W m⁻²
    # Make into daily values!
    all.daily = aggregate(RH672 ~ Date, data = all_data, mean, na.rm = T)
    all.daily.1 = aggregate(cbind(RH1250, RHTd, CO2168, T336, T168) ~ Date, data = all_data, mean, na.rm = T, na.action = NULL)
    all.daily.2 = aggregate(cbind(tsoil_5, tsoil_10, wsoil_B1, wsoil_B2) ~ Date, data = all_data, mean, na.rm = T, na.action = NULL)
    # µmol m⁻² s⁻¹
    all.daily.3 = aggregate(HYY_EDDY233.GPP ~ Date, data = GPP, mean, na.rm = T, na.action = NULL)
    all.daily.4 = aggregate(cbind(wpsoil_A, wpsoil_B) ~ Date, data = all_data, mean, na.rm = T, na.action = NULL)
    all.daily.sum = aggregate(cbind(PAR, Precip) ~ Date, data = all_data, sum, na.rm = F, na.action = NULL)
    # Here the data looks strange for 2017-2018 two years, so made into NA values
    all.daily.sum$Precip[lubridate::year(all.daily.sum$Date) %in% 2017:2018] <- NA
    all.daily.sum$Glob = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob~Date, data = Glob, sum, na.rm = T, na.action = NULL)$HYY_META.Glob)
    all.daily.sum$Glob67 = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob67, sum, na.rm = T, na.action = NULL)$HYY_META.Glob67)
    all.daily.mean = aggregate(HYY_META.Glob~Date, data = Glob, mean, na.rm = T, na.action = NULL)
    all.daily.mean$Glob_mean = bigleaf::Rg.to.PPFD(all.daily.mean$HYY_META.Glob)
    all.daily.mean$Glob67_mean = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob67, mean, na.rm = T, na.action = NULL)$HYY_META.Glob67)
    all.daily.max = aggregate(HYY_META.Glob~Date, data = Glob, max, na.rm = T, na.action = NULL)
    all.daily.max$Glob_max = bigleaf::Rg.to.PPFD(all.daily.max$HYY_META.Glob)
    all.daily.max$Glob67_max = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob67, max, na.rm = T, na.action = NULL)$HYY_META.Glob67)

    install.packages("ProfoundData")

    weather_original_2015 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2015.csv", header = T, sep = ",")
    weather_original_2016 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2016.csv", header = T, sep = ",")
    weather_original_2017 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2017.csv", header = T, sep = ",")
    weather_original_2018 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2018.csv", header = T, sep = ",")
    weather_original = rbind(rbind(rbind(weather_original_2015, weather_original_2016), weather_original_2017), weather_original_2018)

    par(mfrow = c(3, 4))
    # TODO: this should be below as data_format is definied later!
    plot(all.daily.1$T336[substring(all.daily.1$Date, 1, 4) %in% c(2015:2018)], weather_original$T, ylab = "From Pauliiina", xlab = "My Data", main = "Temperature")
    points(all.daily.1$T168[substring(all.daily.1$Date, 1, 4) %in% c(2015:2018)], weather_original$T, col = "blue")
    abline(0, 1, col = "red")
    plot(all.daily.1$T336[substring(all.daily.1$Date, 1, 4) %in% c(2015:2018)] - weather_original$T, ylab = "From Pauliiina", xlab = "My Data")
    points(all.daily.1$T168[substring(all.daily.1$Date, 1, 4) %in% c(2015:2018)] - weather_original$T, ylab = "From Pauliiina", xlab = "My Data", col = "blue")

    # Doesn't plot as values are NA throughout
    plot(all.daily.sum$Precip[substring(all.daily.sum$Date, 1, 4) %in% c(2015:2018)], weather_original$Rain, ylab = "From Pauliiina", xlab = "My Data", main = "Precipitation")
    abline(0, 1, col = "red")
    plot(all.daily.sum$Precip[substring(all.daily.sum$Date, 1, 4) %in% c(2015:2018)] - weather_original$Rain, ylab = "From Pauliiina", xlab = "My Data")

    # Weather_original doesn't have PAR
    #plot(all.daily.sum$Glob67[substring(all.daily.sum$Date, 1, 4) %in% c(2015:2018)], weather_original$PAR, ylab = "From Pauliiina", xlab = "My Data", main = "PAR")
    #points(all.daily.sum$Glob[substring(all.daily.sum$Date, 1, 4) %in% c(2015:2018)], weather_original$PAR, col = "blue")
    #abline(0, 1, col = "red")
    #plot(all.daily.sum$Glob67[substring(all.daily.sum$Date, 1, 4) %in% c(2015:2018)] - weather_original$PAR, ylab = "From Pauliiina", xlab = "My Data")
    #points(all.daily.sum$Glob[substring(all.daily.sum$Date, 1, 4) %in% c(2015:2018)] - weather_original$PAR, col = "blue")

    plot(all.daily.2$wsoil_B2[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)], weather_original$MB, ylab = "From Pauliiina", xlab = "My Data", main = "MB")
    points(all.daily.2$wsoil_B1[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)], weather_original$MB, col = "blue")
    abline(0, 1, col = "red")
    plot(all.daily.2$wsoil_B2[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)] - weather_original$MB, ylab = "From Pauliiina", xlab = "My Data")
    points(all.daily.2$wsoil_B1[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)] - weather_original$MB, col = "blue")

    plot(all.daily.2$tsoil_5[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)], weather_original$TSA, ylab = "From Pauliiina", xlab = "My Data", main = "TSA")
    abline(0, 1, col = "red")
    plot(all.daily.2$tsoil_5[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)] - weather_original$TSA, ylab = "From Pauliiina", xlab = "My Data")

    plot(all.daily.2$tsoil_10[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)], weather_original$TSB, ylab = "From Pauliiina", xlab = "My Data", main = "TSB")
    abline(0, 1, col = "red")
    plot(all.daily.2$tsoil_10[substring(all.daily.2$Date, 1, 4) %in% c(2015:2018)] - weather_original$TSB, ylab = "From Pauliiina", xlab = "My Data")

    # PAR data sorting
    PAR$Missing <-is.nan(PAR$HYY_META.PAR)
    Glob$MissingGlob <-is.nan(Glob$HYY_META.Glob)
    Glob$MissingGlob67 <-is.nan(Glob67$HYY_META.Glob67)
    Glob$HYY_META.Glob67 <- Glob67$HYY_META.Glob67
    PAR$YMD <- as.Date(paste(PAR$Year, PAR$Month, PAR$Day, sep = "-"))
    Glob$YMD <- as.Date(paste(Glob$Year, Glob$Month, Glob$Day, sep = "-"))
    Daily_Missing <- aggregate(Missing~YMD, PAR, sum)
    Daily_Missing_Glob <- aggregate(MissingGlob67~YMD, Glob, sum)
    Daily_Missing_Glob$MissingGlob <- aggregate(MissingGlob~YMD, Glob, sum)$MissingGlob
    Daily_Missing$Percentage <- Daily_Missing$Missing/(24*60)
    # Some days are repeated, but if there is 200% of entries then none of the entries exist anyway.
    Daily_Missing$Percentage[Daily_Missing$Percentage > 1] <- 1
    Daily_Missing_Glob$Percentage <- Daily_Missing_Glob$MissingGlob/(24*60)
    Daily_Missing_Glob$Percentage67 <- Daily_Missing_Glob$MissingGlob67/(24*60)
    second_multiplier <- function(x) {
      multi <- (24 - 1) * (x) + 1
      return(multi)
    }

    # Unit correction, µmol m⁻² s⁻¹ to sum mmol m⁻² day-1
    all.daily.sum$PAR <- 0.000001 * 24 * second_multiplier(Daily_Missing$Percentage) * all.daily.sum$PAR
    all.daily.sum$Glob <- 0.000001 * 24 * second_multiplier(Daily_Missing_Glob$Percentage) * all.daily.sum$Glob
    all.daily.sum$Glob67 <- 0.000001 * 24 * second_multiplier(Daily_Missing_Glob$Percentage67) * all.daily.sum$Glob67
    all.daily.sum$PAR[Daily_Missing$Missing >= 0.75] <- NA
    all.daily.sum$Glob[Daily_Missing_Glob$Missing >= 0.75] <- NA
    all.daily.sum$Glob67[Daily_Missing_Glob$Missing67 >= 0.75] <- NA

    all.gapfill.1 <- merge(all.daily, all.daily.sum, all = T)
    all.gapfill.2 <- merge(all.daily.1, all.daily.2, all = T)
    all.gapfill.1.2 <- merge(all.gapfill.1, all.gapfill.2, all = T)
    all.gapfill.3.4 <- merge(all.daily.3, all.daily.4, all = T)
    all.gapfill.1.4 <- merge(all.gapfill.1.2, all.gapfill.3.4, all = T)
    all.gapfill.mean.max <- merge(all.daily.mean, all.daily.max, all = T)
    all.gapfill <- merge(all.gapfill.mean.max, all.gapfill.1.4, all = T)

    # Gap fill with the other level if possible
    all.gapfill$RH672[is.na(all.gapfill$RH672)] <- all.gapfill$RH1250[is.na(all.gapfill$RH672)] - mean(all.gapfill$RH1250, na.rm = T) + mean(all.gapfill$RH672, na.rm = T)
    all.gapfill$RH1250[is.na(all.gapfill$RH1250)] <- all.gapfill$RH672[is.na(all.gapfill$RH1250)] - mean(all.gapfill$RH672, na.rm = T) + mean(all.gapfill$RH1250, na.rm = T)
    all.gapfill$RH672[is.na(all.gapfill$RH672)] <- all.gapfill$RHTd[is.na(all.gapfill$RH672)] - mean(all.gapfill$RHTd, na.rm = T) + mean(all.gapfill$RH672, na.rm = T)
    all.gapfill$RH1250[is.na(all.gapfill$RH1250)] <- all.gapfill$RHTd[is.na(all.gapfill$RH1250)] - mean(all.gapfill$RHTd, na.rm = T) + mean(all.gapfill$RH1250, na.rm = T)
    all.gapfill$T336[is.na(all.gapfill$T336)] <- all.gapfill$T168[is.na(all.gapfill$T336)] - mean(all.gapfill$T336, na.rm = T) + mean(all.gapfill$T168, na.rm = T)
    all.gapfill$T168[is.na(all.gapfill$T168)] <- all.gapfill$T336[is.na(all.gapfill$T168)] - mean(all.gapfill$T168, na.rm = T) + mean(all.gapfill$T336, na.rm = T)
    all.gapfill$tsoil_5[is.na(all.gapfill$tsoil_5)] <- all.gapfill$tsoil_10[is.na(all.gapfill$tsoil_5)] - mean(all.gapfill$tsoil_10, na.rm = T) + mean(all.gapfill$tsoil_5, na.rm = T)
    all.gapfill$tsoil_10[is.na(all.gapfill$tsoil_10)] <- all.gapfill$tsoil_5[is.na(all.gapfill$tsoil_10)] - mean(all.gapfill$tsoil_5, na.rm = T) + mean(all.gapfill$tsoil_10, na.rm = T)
    all.gapfill$wsoil_B1[is.na(all.gapfill$wsoil_B1)] <- all.gapfill$wsoil_B2[is.na(all.gapfill$wsoil_B1)] - mean(all.gapfill$wsoil_B2, na.rm = T) + mean(all.gapfill$wsoil_B1, na.rm = T)
    all.gapfill$wsoil_B2[is.na(all.gapfill$wsoil_B2)] <- all.gapfill$wsoil_B1[is.na(all.gapfill$wsoil_B2)] - mean(all.gapfill$wsoil_B1, na.rm = T) + mean(all.gapfill$wsoil_B2, na.rm = T)
    all.gapfill$wpsoil_A[is.na(all.gapfill$wpsoil_A)] <- all.gapfill$wpsoil_B[is.na(all.gapfill$wpsoil_A)] - mean(all.gapfill$wpsoil_B, na.rm = T) + mean(all.gapfill$wpsoil_A, na.rm = T)
    all.gapfill$wpsoil_B[is.na(all.gapfill$wpsoil_B)] <- all.gapfill$wpsoil_A[is.na(all.gapfill$wpsoil_B)] - mean(all.gapfill$wpsoil_A, na.rm = T) + mean(all.gapfill$wpsoil_B, na.rm = T)
    # Units for the Glob value already changed earlier in the code so not changed again here!
    all.gapfill$PAR[is.na(all.gapfill$PAR)] <- all.gapfill$Glob67[is.na(all.gapfill$PAR)] - mean(all.gapfill$Glob67, na.rm = T) + mean(all.gapfill$PAR, na.rm = T)
    all.gapfill$PAR[is.na(all.gapfill$PAR)] <- all.gapfill$Glob[is.na(all.gapfill$PAR)] - mean(all.gapfill$Glob, na.rm = T) + mean(all.gapfill$PAR, na.rm = T)
    all.gapfill$Glob_mean[is.na(all.gapfill$Glob_mean)] <- all.gapfill$Glob67_mean[is.na(all.gapfill$Glob_mean)] - mean(all.gapfill$Glob67_mean, na.rm = T) + mean(all.gapfill$Glob_mean, na.rm = T)
    all.gapfill$Glob_max[is.na(all.gapfill$Glob_max)] <- all.gapfill$Glob67_max[is.na(all.gapfill$Glob_max)] - mean(all.gapfill$Glob67_max, na.rm = T) + mean(all.gapfill$Glob_max, na.rm = T)

    # If the gap is smaller than 6 days apply linear interpolation
    all.gapfill[,-1] <- zoo::na.approx(all.gapfill[,-1], na.rm = F, maxgap = 6)

    data_IIASA_smear <- all.gapfill[,c("Date", "T168", "Glob_mean", "Glob_max", "wpsoil_B")]
    data_IIASA_smear$VPD <- bigleaf::rH.to.VPD(0.01*all.gapfill$RH672, all.gapfill$T168)
    data_IIASA_smear$wpsoil_B <- -data_IIASA_smear$wpsoil_B
    if (save) {
      save(data_IIASA_smear, file = paste0(raw.directory, "data_IIASA_smear.RData"))
    }

    data_format = all.gapfill[,c("Date", "T168", "tsoil_5", "tsoil_10", "wsoil_B1", "Precip", "PAR", "CO2168", "RH672")]
    names(data_format) = c("Date", "T", "TSA", "TSB", "MB", "Rain", "PAR", "CO2", "RH")
    data_format$VPD <- bigleaf::rH.to.VPD(0.01*all.gapfill$RH672, all.gapfill$T168) # VPD
    data_format$fAPAR <- rep(0.7, nrow(data_format)) # fAPAR

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

    if (sum(unlist(lapply(data_format, is.na))) != 0) {warning(paste("The dataset still has NAs in. Add some suplimentary data processing!"))}

    if (save) {
      # Saves the data in the same place as the download files
      save(data_format, file = paste0(raw.directory, "Hyde_weather_CASSIA.RData"))
    }

    ## TODO: gapfill the hyytiälä data still

    return(data_format)
  }
}

