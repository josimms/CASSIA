#' @export
Hyytiala_Data_Creation <- function(raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/",
                                   year_start,
                                   year_end,
                                   download = F,
                                   clean_data = F,
                                   save = T) {


  ####
  # Clean the data with other Hyytiala data
  ####

  if (clean_data) {
    # Here the data looks strange for 2017-2018 two years, so made into NA values
    all.daily.sum$Glob = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob~Date, data = Glob, sum, na.rm = T, na.action = NULL)$HYY_META.Glob)
    all.daily.sum$Glob67 = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob67, sum, na.rm = T, na.action = NULL)$HYY_META.Glob67)
    all.daily.mean = aggregate(HYY_META.Glob~Date, data = Glob, mean, na.rm = T, na.action = NULL)
    all.daily.mean$Glob_mean = bigleaf::Rg.to.PPFD(all.daily.mean$HYY_META.Glob)
    all.daily.mean$Glob67_mean = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob67, mean, na.rm = T, na.action = NULL)$HYY_META.Glob67)
    all.daily.max = aggregate(HYY_META.Glob~Date, data = Glob, max, na.rm = T, na.action = NULL)
    all.daily.max$Glob_max = bigleaf::Rg.to.PPFD(all.daily.max$HYY_META.Glob)
    all.daily.max$Glob67_max = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob67, max, na.rm = T, na.action = NULL)$HYY_META.Glob67)

    # Weather from Hyytiälä used for the original CASSIA runs - new weather tested against these old files to make sure that they are consistent
    weather_original_2015 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2015.csv", header = T, sep = ",")
    weather_original_2016 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2016.csv", header = T, sep = ",")
    weather_original_2017 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2017.csv", header = T, sep = ",")
    weather_original_2018 = read.csv(file = "/home/josimms/Documents/CASSIA/data/weather_original_2018.csv", header = T, sep = ",")
    weather_original = rbind(rbind(rbind(weather_original_2015, weather_original_2016), weather_original_2017), weather_original_2018)

    ### Plot the data to make sure it makes sense
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
    # As there are missing values and we aim to have an aaverage we need to make the average for the seconds that are missing per day
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
    all.daily.sum$Glob_umol <- all.daily.sum$Glob
    all.daily.sum$Glob67_umol <- all.daily.sum$Glob67
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

