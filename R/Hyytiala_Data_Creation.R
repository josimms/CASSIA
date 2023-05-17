#' @export
Hyytiala_Data_Creation <- function(raw.directory,
                                   year_start,
                                   year_end,
                                   download = F,
                                   clean_data = F,
                                   save = T) {

  ####
  # Download the data from the SMEAR database
  ####

  if (download) {
    http.origin = "https://smear-backend.rahtiapp.fi/search/timeseries/csv?tablevariable=HYY_"
    for (variable in c(paste0("META.", c("RH672", "RH1250", "PAR", "CO2168", "T336", "T168", "Precip", "tsoil_5", "tsoil_10", "wsoil_B1", "Glob", "Pamb336")), "EDDY233.GPP")) {
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
    ####
    RH672 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RH672")), read.csv, dec = "."))
    RH1250 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RH1250")), read.csv, dec = "."))
    PAR <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "PAR")), read.csv, dec = "."))
    CO2168 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "CO2168")), read.csv, dec = "."))
    T336 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "T336")), read.csv, dec = "."))
    T168 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "T168")), read.csv, dec = "."))
    Precip <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Precip")), read.csv, dec = "."))
    tsoil_5 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "tsoil_5")), read.csv, dec = "."))
    tsoil_10 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "tsoil_10")), read.csv, dec = "."))
    wsoil_B1 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wsoil_B1")), read.csv, dec = "."))
    wsoil_B2 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wsoil_B2")), read.csv, dec = "."))
    GPP <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "GPP")), read.csv, dec = "."))
    GPP$Date <- as.Date(paste(GPP$Year, GPP$Month, GPP$Day, sep = "-"), "%Y-%m-%e")
    Glob <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Glob")), read.csv, dec = "."))
    Glob$Date <- as.Date(paste(Glob$Year, Glob$Month, Glob$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")
    CO2168$Date <- as.Date(paste(CO2168$Year, CO2168$Month, CO2168$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")

    all_data <- as.data.frame(cbind(RH672$HYY_META.RH672, RH1250$HYY_META.RH1250, PAR$HYY_META.PAR,
                                    CO2168$HYY_META.CO2168, T336$HYY_META.T336, T168$HYY_META.T168,
                                    Precip$HYY_META.Precip,
                                    tsoil_5$HYY_META.tsoil_5, tsoil_10$HYY_META.tsoil_10,
                                    wsoil_B1$HYY_META.wsoil_B1, wsoil_B2$HYY_META.wsoil_B2))

    names(all_data) <- c("RH672", "RH1250", "PAR", "CO2168", "T336", "T168",
                         "Precip", "tsoil_5", "tsoil_10", "wsoil_B1", "wsoil_B2")
    all_data$Date <- as.Date(paste(RH672$Year, RH672$Month, RH672$Day, sep = "-"), format = "%Y-%m-%d", tz = "BST")
    # Make into daily values!
    all.daily = aggregate(RH672 ~ Date, data = all_data, mean, na.rm = T, na.action = NULL)
    all.daily.1 = aggregate(cbind(RH1250, CO2168, T336, T168) ~ Date, data = all_data, mean, na.rm = T, na.action = NULL)
    all.daily.2 = aggregate(cbind(tsoil_5, tsoil_10, wsoil_B1, wsoil_B2) ~ Date, data = all_data, mean, na.rm = T, na.action = NULL)
    # µmol m⁻² s⁻¹
    all.daily$GPP = aggregate(HYY_EDDY233.GPP ~ Date, data = GPP, mean, na.rm = T, na.action = NULL)$HYY_EDDY233.GPP
    all.daily.sum = aggregate(cbind(PAR, Precip) ~ Date, data = all_data, sum, na.action = NULL)
    # Here the data looks strange for 2017-2018 two years, so made into NA values
    all.daily.sum$Precip[lubridate::year(all.daily.sum$Date) %in% 2017:2018] <- NA
    all.daily.sum$Glob = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob~Date, data = Glob, sum, na.rm = T, na.action = NULL)$HYY_META.Glob)
    all.daily.sum$Glob67 = bigleaf::Rg.to.PPFD(aggregate(HYY_META.Glob67~Date, data = Glob, sum, na.rm = T, na.action = NULL)$HYY_META.Glob67)

    # PAR data sorting
    PAR$Missing <-is.nan(PAR$HYY_META.PAR)
    Glob$MissingGlob <-is.nan(Glob$HYY_META.Glob)
    Glob$MissingGlob67 <-is.nan(Glob$HYY_META.Glob67)
    PAR$YMD <- as.Date(paste(PAR$Year, PAR$Month, PAR$Day, sep = "-"))
    Glob$YMD <- as.Date(paste(Glob$Year, Glob$Month, Glob$Day, sep = "-"))
    Daily_Missing <- aggregate(Missing~YMD, PAR, sum)
    Daily_Missing_Glob <- aggregate(MissingGlob~YMD, Glob, sum)
    Daily_Missing_Glob$MissingGlob67 <- aggregate(MissingGlob67~YMD, Glob, sum)$MissingGlob67
    Daily_Missing$Percentage <- Daily_Missing$Missing/(24*60*60)
    Daily_Missing_Glob$Percentage <- Daily_Missing_Glob$MissingGlob/(24*60*60)
    Daily_Missing_Glob$Percentage67 <- Daily_Missing_Glob$MissingGlob67/(24*60*60)
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
    all.gapfill <- merge(all.gapfill.1, all.gapfill.2, all = T)

    # Gap fill with the other level if possible
    all.gapfill$RH672[is.na(all.gapfill$RH672)] <- all.gapfill$RH1250[is.na(all.gapfill$RH672)] - mean(all.gapfill$RH1250, na.rm = T) + mean(all.gapfill$RH672, na.rm = T)
    all.gapfill$RH1250[is.na(all.gapfill$RH1250)] <- all.gapfill$RH672[is.na(all.gapfill$RH1250)] - mean(all.gapfill$RH672, na.rm = T) + mean(all.gapfill$RH1250, na.rm = T)
    all.gapfill$T336[is.na(all.gapfill$T336)] <- all.gapfill$T168[is.na(all.gapfill$T336)] - mean(all.gapfill$T336, na.rm = T) + mean(all.gapfill$T168, na.rm = T)
    all.gapfill$T168[is.na(all.gapfill$T168)] <- all.gapfill$T336[is.na(all.gapfill$T168)] - mean(all.gapfill$T168, na.rm = T) + mean(all.gapfill$T336, na.rm = T)
    all.gapfill$tsoil_5[is.na(all.gapfill$tsoil_5)] <- all.gapfill$tsoil_10[is.na(all.gapfill$tsoil_5)] - mean(all.gapfill$tsoil_10, na.rm = T) + mean(all.gapfill$tsoil_5, na.rm = T)
    all.gapfill$tsoil_10[is.na(all.gapfill$tsoil_10)] <- all.gapfill$tsoil_5[is.na(all.gapfill$tsoil_10)] - mean(all.gapfill$tsoil_5, na.rm = T) + mean(all.gapfill$tsoil_10, na.rm = T)
    all.gapfill$wsoil_B1[is.na(all.gapfill$wsoil_B1)] <- all.gapfill$wsoil_B2[is.na(all.gapfill$wsoil_B1)] - mean(all.gapfill$wsoil_B2, na.rm = T) + mean(all.gapfill$wsoil_B1, na.rm = T)
    all.gapfill$wsoil_B2[is.na(all.gapfill$wsoil_B2)] <- all.gapfill$wsoil_B1[is.na(all.gapfill$wsoil_B2)] - mean(all.gapfill$wsoil_B1, na.rm = T) + mean(all.gapfill$wsoil_B2, na.rm = T)
    # Units for the Glob value already changed earlier in the code so not changed again here!
    all.gapfill$PAR[is.na(all.gapfill$PAR)] <- all.gapfill$Glob[is.na(all.gapfill$PAR)] - mean(all.gapfill$Glob, na.rm = T) + mean(all.gapfill$PAR, na.rm = T)
    all.gapfill$PAR[is.na(all.gapfill$PAR)] <- all.gapfill$Glob67[is.na(all.gapfill$PAR)] - mean(all.gapfill$Glob67, na.rm = T) + mean(all.gapfill$PAR, na.rm = T)

    # If the gap is smaller than 6 days apply linear interpolation
    all.gapfill[,-1] <- zoo::na.approx(all.gapfill[,-1], na.rm = F, maxgap = 6)

    data_format = all.gapfill[,c("Date", "T168", "tsoil_5", "tsoil_10", "wsoil_B1", "Precip", "PAR", "CO2168", "RH672")]
    names(data_format) = c("Date", "T", "TSA", "TSB", "MB", "Rain", "PAR", "CO2", "RH")
    data_format$VPD <- bigleaf::rH.to.VPD(0.01*all.gapfill$RH672, all.gapfill$T168) # VPD
    data_format$fAPAR <- rep(0.7, nrow(data_format)) # fAPAR

    if (sum(unlist(lapply(data_format, is.na))) != 0) {warning(paste("The dataset still has NAs in. Add some suplimentary data processing!"))}

    if (save) {
      # Saves the data in the same place as the download files
      save(data_format, file = paste0(raw.directory, "Hyde_weather_CASSIA.RData"))
    }
    return(data_format)
  }
}
