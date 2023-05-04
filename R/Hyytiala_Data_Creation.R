#' @export

Hyytiala_Data_Creation <- function(raw.directory = "/home/joanna/Asiakirjat/Hyytiälä/",
                                   year_start = 2010,
                                   year_end = 2022,
                                   download = F,
                                   clean_data = F) {

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
    RH672 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RH672")), read.csv))
    RH1250 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "RH1250")), read.csv))
    PAR <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "PAR")), read.csv))
    CO2168 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "CO2168")), read.csv))
    T336 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "T336")), read.csv))
    T168 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "T168")), read.csv))
    Precip <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Precip")), read.csv))
    tsoil_5 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "tsoil_5")), read.csv))
    tsoil_10 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "tsoil_10")), read.csv))
    wsoil_B1 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "wsoil_B1")), read.csv))
    Glob <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Glob")[1:3]), read.csv))
    Pamb336 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Pamb336")), read.csv))
    Pamb0 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Pamb0")), read.csv))
    GPP <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "GPP")), read.csv))
    GPP$YMD <- as.Date(paste(GPP$Year, GPP$Month, GPP$Day, sep = "-"), "%Y-%m-%e")
    Globmast <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Globmast")), read.csv))
    Glob67 <- dplyr::bind_rows(lapply(paste0(raw.directory, list.files(raw.directory, "Glob67")), read.csv))

    all_data <- as.data.frame(cbind(RH672$Year, RH672$Month, RH672$Day, RH672$Hour, RH672$Minute, RH672$Second,
                                    RH672$HYY_META.RH672, RH1250$HYY_META.RH1250, PAR$HYY_META.PAR,
                                    CO2168$HYY_META.CO2168, T336$HYY_META.T336, T168$HYY_META.T168,
                                    Glob$HYY_META.Glob, Pamb336$HYY_META.Pamb336, Pamb0$HYY_META.Pamb0,
                                    Precip$HYY_META.Precip, tsoil_5$HYY_META.tsoil_5, tsoil_10$HYY_META.tsoil_10,
                                    wsoil_B1$HYY_META.wsoil_B1, GPP$HYY_EDDY233.GPP))
    names(all_data) <- c("Year", "Month", "Day", "Hour", "Minute", "Second",
                         "RH672", "RH1250", "PAR", "CO2168", "T336", "T168",
                         "Glob", "Pamb336", "Pamb0", "Precip",
                         "tsoil_5", "tsoil_10", "wsoil_B1", "GPP")
    # Make into daily values!
    all_data$YMD = as.Date(paste(all_data$Year, all_data$Month, all_data$Day, sep = "-"), "%Y-%m-%e")
    all.daily = aggregate(cbind(RH672, RH1250, CO2168, T336, T168, tsoil_5, tsoil_10, wsoil_B1) ~ YMD, data = all_data, mean, na.rm = T, na.action = NULL)
    # µmol m⁻² s⁻¹
    all.daily$GPP = aggregate(HYY_EDDY233.GPP ~ YMD, data = GPP, mean, na.rm = T, na.action = NULL)$HYY_EDDY233.GPP
    all.daily.sum = aggregate(cbind(PAR, Precip) ~ YMD, data = all_data, sum, na.rm = T, na.action = NULL)
    # Here the data looks strange for the first two years, so made into NA values
    all.daily.sum$Precip[1:(365*2)] <- NA

    # PAR data sorting
    PAR$Missing <-is.nan(PAR$HYY_META.PAR)
    PAR$YMD <- as.Date(paste(PAR$Year, PAR$Month, PAR$Day, sep = "-"))
    Daily_Missing <- aggregate(Missing~YMD, PAR, sum)
    Daily_Missing$Percentage <- Daily_Missing$Missing/(24*60*60)
    second_multiplier <- function(x) {
      multi <- (24 - 1) * (x) + 1
      return(multi)
    }
    # Unit correction, µmol m⁻² s⁻¹ to sum mmol m⁻² day-1
    all.daily.sum$PAR <- 0.000001 * 24 * second_multiplier(Daily_Missing$Percentage) * all.daily.sum$PAR
    all.daily.sum$PAR[Daily_Missing$Missing >= 0.75] <- NA

    all.gapfill <- merge(all.daily, all.daily.sum, all = T)

    # Gap fill with the other level if possible
    all.gapfill$RH672[is.na(all.gapfill$RH672)] <- all.gapfill$RH1250[is.na(all.gapfill$RH672)]
    all.gapfill$RH1250[is.na(all.gapfill$RH1250)] <- all.gapfill$RH672[is.na(all.gapfill$RH1250)]

    # If the gap is smaller than 4 days apply linear interpolation - Pauliina
    all.gapfill[,-c(1, 10)] <- zoo::na.approx(all.gapfill[,-c(1, 10)], na.rm = F, maxgap = 6)

    save(all.gapfill, file = paste0(raw.directory, "Hyde_weather.RData"))
  }


}
