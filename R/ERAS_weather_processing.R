# DEW POINT

calculate_VPD <- function(dew_point_temp, air_temp) {
  # Calculate saturation vapor pressure at air temperature (in hPa)
  e_s <- 6.11 * 10^(7.5 * air_temp / (237.3 + air_temp))

  # Calculate vapor pressure at dew point temperature (in hPa)
  e_d <- 6.11 * 10^(7.5 * dew_point_temp / (237.3 + dew_point_temp))

  # Calculate relative humidity (in percentage)
  RH <- (e_d / e_s)

  # Calculate VPD
  VPD <- 0.1 * bigleaf::rH.to.VPD(RH, air_temp) # Units kPa to hPa
  return(list("RH" = RH,
              "VPD" = VPD))
}

# READING NC FILE

# Plotting function (reuse for both monthly and daily)
plot_data <- function(data, title_prefix, monthly) {
  par(mfrow = c(3, 3))
  if (monthly) {
    data$YM <- paste0(data$YM, "-01")
    data$date <- as.Date(data$YM, format = "%Y-%m-%d")
  } else {
    data$date <- as.Date(data$YMD, format = "%Y-%m-%d")
  }

  # "vpd",
  plot_vars <- c("Temp", "PPFD", "PPFD_max", "swvl2", "VPD", "SWP", "co2")
  for (var in plot_vars) {
    plot(data$date, data[[var]], main = paste(title_prefix, var),
         xlab = "Dates", ylab = switch(var,
                                       Temp = "Degrees C",
                                       PPFD = "umol m-2 s-1",
                                       PPFD_max = "umol m-2 s-1",
                                       swvl2 = "-kPa",
                                       VPD = "hPa",
                                       SWP = "- Ma",
                                       co2 = "ppm"))
    title(sub = sprintf("Percentage missing: %.2f%%", 100 * sum(is.na(data[[var]])) / nrow(data)))
  }
  plot(data$PPFD, data$PPFD_max, main = "Global: Mean vs Max", xlab = "Mean", ylab = "Max")
  abline(lm(PPFD_max ~ PPFD, data = data), col = "red", lwd = 2)
  title(sub = sprintf("Gradient: %.2f", coef(lm(PPFD_max ~ PPFD, data = data))[2]))
}

###
# ERAS Reading NC
###

ERAS_reading_nc <- function(path_nc = "/home/josimms/Documents/Austria/eras_data",
                            path_test = "/home/josimms/Documents/Austria/Plant-FATE/tests/data",
                            raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/") {

  library(dplyr)

  # "2m_dewpoint_temperature", 2d
  # "2m_temperature", 2t
  # "total_precipitation", tp
  # "surface_solar_radiation_downwards", ssrd
  # "soil_temperature_level_1", stl1
  # "soil_temperature_level_2", stl2
  # "volumetric_soil_water_layer_1", swvl1
  # "volumetric_soil_water_layer_2" swvl2
  variables <- c("tp", "ssrd", "d2m", "t2m", "stl1", "stl2", "swvl1", "swvl2")

  # Pre-allocate lists
  nc_data <- vector("list", length(variables))
  names(nc_data) <- c(variables)

  # Import the files
  nc_files_accum <- list.files(path = path_nc, pattern = "_era5_accum.nc", full.names = TRUE)
  nc_files_inst <- list.files(path = path_nc, pattern = "_era5_instant.nc", full.names = TRUE)

  # Coordinates
  lon_target <- 24.29477
  lat_target <- 61.84741

  # Extract the data for the right location in the files
  dataset_cds_raw <- list()
  for (i in 1:length(nc_files_accum)) {
    nc_accum <- ncdf4::nc_open(nc_files_accum[i]) # ssrd, tp
    nc_inst <- ncdf4::nc_open(nc_files_inst[i]) # others

    lon <- ncdf4::ncvar_get(nc, "longitude")
    lat <- ncdf4::ncvar_get(nc, "latitude")
    lon_index <- which.min(abs(lon - lon_target))
    lat_index <- which.min(abs(lat - lat_target))

    dates_accum <- ncdf4::ncvar_get(nc_accum, nc_accum$dim$valid_time)
    dates_inst <- ncdf4::ncvar_get(nc_inst, nc_inst$dim$valid_time)

    # cat("Year", i+1959)
    # print(length(dates_accum) == length(dates_inst))
    for (var in variables[1:2]) {
      nc_data[[var]] <- ncdf4::ncvar_get(nc_accum, varid = var,
                                         start = c(lon_index, lat_index, 1),
                                         count = c(1, 1, -1))
      # print(length(nc_data[[var]]))
    }
    for (var in variables[3:length(variables)]) {
      nc_data[[var]] <- ncdf4::ncvar_get(nc_inst, varid = var,
                                         start = c(lon_index, lat_index, 1),
                                         count = c(1, 1, -1))
      # print(length(nc_data[[var]]))
    }

    ncdf4::nc_close(nc_accum)
    ncdf4::nc_close(nc_inst)

    # Create data.table
    dataset_cds_raw_year <- data.table::data.table(
      Temp = nc_data$t2m  - 273.15, # 'C
      Temp_Dew = nc_data$d2m  - 273.15, # 'C
      Precip = 1000 * nc_data$tp, # m to mm
      PPFD = bigleaf::Rg.to.PPFD(nc_data$ssrd/(60*60)), # W m-2 to umol m-2 s-1, PPFD (daily 24-hr mean)
      Temp_Soil_1 = nc_data$stl1 - 273.15, # 'C
      Temp_Soil_2 = nc_data$stl2 - 273.15, # 'C
      TotGlob = nc_data$ssrd/(60*60), # W m-2 # TODO: units
      PAR_preles = 0.000001 * 86400 * bigleaf::Rg.to.PPFD(nc_data$ssrd/(60*60)), # W m-2 to sum umol m-2 day-1
      swvl1 = nc_data$swvl1, # m**3 m**-3
      swvl2 = nc_data$swvl2, # m**3 m**-3
      date = as.POSIXct(dates_accum, origin = "1970-01-01 00:00:00", tz = "GMT")
    )
    dataset_cds_raw[[i]] <- dataset_cds_raw_year
  }

  # Make this list into a datatable
  dataset_cds_raw_all <- data.table::rbindlist(dataset_cds_raw)

  # Use the VPD function above and format the dates so monthly and daily averages can be made
  dataset_cds_raw_all[, VPD := calculate_VPD(Temp_Dew, Temp)$VPD] # hPa
  dataset_cds_raw_all[, RH := calculate_VPD(Temp_Dew, Temp)$RH] # decimal percentage
  dataset_cds_raw_all[, YM := format(date, "%Y-%m")]
  dataset_cds_raw_all[, YMD := format(date, "%Y-%m-%d")]
  dataset_cds_raw_all[, Year := as.numeric(format(date, "%Y"))]
  dataset_cds_raw_all[, Month := as.numeric(format(date, "%m"))]

  # Monthly aggregation
  monthy_dataset <- dataset_cds_raw_all[, lapply(.SD, mean), by = YM, .SDcols = -c("YMD", "date")]
  monthy_dataset[, PPFD_max := dataset_cds_raw_all[, .(PPFD_max = max(PPFD)), by = YM]$PPFD_max]

  # Daily aggregation
  daily_dataset <- dataset_cds_raw_all[, lapply(.SD, mean), by = YMD, .SDcols = -c("YM", "date")]
  daily_dataset[, PPFD_max := dataset_cds_raw_all[, .(PPFD_max = max(PPFD)), by = YMD]$PPFD_max]

  # Save datasets for only the ERA5 data
  data.table::fwrite(monthy_dataset, file = file.path(path_test, "montly_dataset.csv"))
  data.table::fwrite(daily_dataset, file = file.path(path_test, "daily_dataset.csv"))

  ####
  # Generate the weather file in the right format
  ####

  # Generate a soil water function to make up for the lack of data for the soil water
    # This baseline is from Hyyitälä data

    # Read the files, make into a dataset not a list then add the MD date aggregation to form a baseline
  soil_water_potential_list <- list()
  count = 1
  for (var in c("wpsoil_A", "wpsoil_B", "GPP")) {
    soil_water_potential_list[[count]] <- data.table::rbindlist(lapply(paste0(raw.directory, list.files(raw.directory, var)), data.table::fread))
    count = count + 1
  }
  soil_water_potential <- data.table::rbindlist(soil_water_potential_list, fill = TRUE)

  soil_water_potential[, MD := paste(soil_water_potential$Month,
                                     soil_water_potential$Day, sep = "-")]
  # Gapfil one level with the other level
  soil_water_potential$HYY_META.wpsoil_B[is.na(soil_water_potential$HYY_META.wpsoil_B)] = soil_water_potential$HYY_META.wpsoil_A[is.na(soil_water_potential$HYY_META.wpsoil_B)] -
    mean(soil_water_potential$HYY_META.wpsoil_A, na.rm = T) +
    mean(soil_water_potential$HYY_META.wpsoil_B, na.rm = T)

  # Monthly aggregation
  soil_water_potential[, YM := paste(soil_water_potential$Year, sprintf("%02d", soil_water_potential$Month), sep = "-")]
  soil_water_potential[, YMD := paste(soil_water_potential$Year, soil_water_potential$Month, soil_water_potential$Day, sep = "-")]
  gpp <- soil_water_potential[, lapply(.SD, mean, na.rm = T), by = YM, .SDcols = -c("MD", "YM", "YMD")]
  names(gpp) <- gsub("HYY_EDDY233.", "", names(gpp))

  # Monthly aggregation - for only one year!
  soil_water_potential_montly <- soil_water_potential[, lapply(.SD, mean, na.rm = T), by = Month, .SDcols = -c("MD", "YM", "YMD")]

  # Daily aggregation
  soil_water_potential_daily <- soil_water_potential[, lapply(.SD, mean, na.rm = T), by = MD,  .SDcols = -c("YM", "YMD")]

  # Load the data
  Hyytiala <- data.table::fread("./data/daily_dataframe.csv")

  # Extract year and month, and calculate the monthly max Glob for each year
  Hyytiala_monthly <- Hyytiala %>%
    group_by(Year, Month) %>%
    summarise(Mean_Temp = mean(T336_mean, na.rm = T), # 'C
              Mean_VPD = mean(VPD, na.rm = T), # kPa
              Mean_PPFD = mean(PAR_mean, na.rm = T), # µmol m⁻² s⁻¹
              Mean_PPFD_max = mean(PAR_max, na.rm = T)) %>% # max µmol m⁻² s⁻¹
    mutate_all(~replace(., is.infinite(.), NA)) %>%
    group_by(Month) %>%
    summarise(Mean_Temp = mean(Mean_Temp, na.rm = TRUE), # 'C
              Mean_VPD = mean(Mean_VPD, na.rm = TRUE), # kPa
              Mean_PPFD = mean(Mean_PPFD, na.rm = TRUE), # µmol m⁻² s⁻¹
              Mean_PPFD_max = mean(Mean_PPFD_max, na.rm = TRUE)) # max µmol m⁻² s⁻¹

  Hyytiala_daily <- Hyytiala %>%
    group_by(Year, Month, Day) %>%
    summarise(Mean_Temp = mean(T336_mean, na.rm = T), # 'C
              Mean_STemp = mean(tsoil_5_mean, na.rm = T), # 'C
              Mean_VPD = mean(VPD, na.rm = T), # hPa
              Mean_RH = mean(RH672_mean, na.rm = T), # %
              Mean_SWC = mean(wsoil_B1_mean, na.rm = T), # m³ m⁻³
              Mean_PAR = mean(PAR_mean, na.rm = T), # µmol m⁻² s⁻¹
              Mean_PPFD = mean(PAR_mean, na.rm = T), # max µmol m⁻² s⁻¹
              Mean_PPFD_max = mean(PAR_max, na.rm = T), # max µmol m⁻² s⁻¹
              Mean_Glob = mean(Glob_mean, na.rm = T), # W m⁻²
              Mean_Precip = mean(Precip_sum, na.rm = T)) %>% # mm
    mutate_all(~replace(., is.infinite(.), NA)) %>%
    group_by(Month, Day) %>%
    summarise(Mean_Temp = mean(Mean_Temp, na.rm = TRUE), # 'C
              Mean_STemp = mean(Mean_STemp, na.rm = T), # 'C
              Mean_VPD = mean(Mean_VPD, na.rm = TRUE), # hPa
              Mean_RH = mean(Mean_RH, na.rm = T), # %
              Mean_SWC = mean(Mean_SWC, na.rm = T), # m³ m⁻³
              Mean_PAR = mean(Mean_PAR, na.rm = T), # µmol m⁻² s⁻¹
              Mean_PPFD = mean(Mean_PPFD, na.rm = TRUE), # µmol m⁻² s⁻¹
              Mean_PPFD_max = mean(Mean_PPFD_max, na.rm = TRUE), # max µmol m⁻² s⁻¹
              Mean_Glob = mean(Mean_Glob, na.rm = T), # W m⁻²
              Mean_Precip = mean(Mean_Precip, na.rm = T)) # mm

  Hyytiala_Prebas <- data.table::fread("./data/preles_smear_CASSIA_ready.csv")
  Hyytiala_Prebas[, Year := substring(dates, 1, 4)]
  Hyytiala_Prebas[, Month := substring(dates, 6, 7)]
  Hyytiala_Prebas[, Day := substring(dates, 9, 10)]
  Hyytiala_daily_prebas <- Hyytiala_Prebas %>%
    group_by(Year, Month, Day) %>%
    summarise(Mean_Temp = mean(T, na.rm = T), # 'C
              Mean_VPD = mean(VPD, na.rm = T), # kPa
              Mean_PAR = mean(PAR, na.rm = T), # µmol m⁻² s⁻¹
              Mean_Precip = mean(Rain, na.rm = T)) %>% # mm
    mutate_all(~replace(., is.infinite(.), NA)) %>%
    group_by(Month, Day) %>%
    summarise(Mean_Temp = mean(Mean_Temp, na.rm = TRUE), # 'C
              Mean_VPD = mean(Mean_VPD, na.rm = TRUE), # kPa
              Mean_PAR = mean(Mean_PAR, na.rm = T), # µmol m⁻² s⁻¹
              Mean_Precip = mean(Mean_Precip, na.rm = T)) # mm

  # Monthly

  plantfate_monthy_dataset <- monthy_dataset
  plantfate_monthy_dataset$co2 <- 380 # ppm TODO: if time get the values from Hyytiala like in the SWP
  plantfate_monthy_dataset$SWP <- - 0.001 * rep(soil_water_potential_montly$HYY_META.wpsoil_B,
                                                length.out = nrow(plantfate_monthy_dataset)) # Soil water potential kPa to - MPa
  plantfate_monthy_dataset$Decimal_year <- seq(plantfate_monthy_dataset$Year[1],
                                               plantfate_monthy_dataset$Year[nrow(plantfate_monthy_dataset)],
                                               length.out = nrow(plantfate_monthy_dataset))
  plantfate_monthy_dataset <- merge(plantfate_monthy_dataset,
                                    gpp[gpp$Year < plantfate_monthy_dataset$Year[nrow(plantfate_monthy_dataset)],c("YM", "GPP")],
                                    by = "YM",
                                    all.x = T)

  monthly_means <- plantfate_monthy_dataset %>%
    group_by(Month) %>%
    summarise(Mean_Temp = mean(Temp, na.rm = TRUE),
              Mean_VPD = mean(VPD, na.rm = TRUE),
              Mean_PPFD = mean(PPFD, na.rm = TRUE),
              Mean_PPFD_max = mean(PPFD_max, na.rm = TRUE))

  cols = c("Mean_Temp", "Mean_VPD", "Mean_PPFD", "Mean_PPFD_max")
  error = Hyytiala_monthly[,cols] - monthly_means[,cols]

  ## Daily

  plantfate_daily_dataset <- daily_dataset
  plantfate_daily_dataset$co2 <- 380
  plantfate_daily_dataset$SWP <- -0.001 * rep(soil_water_potential_daily$HYY_META.wpsoil_B,
                                              length.out = nrow(plantfate_daily_dataset))
  plantfate_daily_dataset$Decimal_year <- seq(plantfate_daily_dataset$Year[1],
                                              plantfate_daily_dataset$Year[nrow(plantfate_daily_dataset)],
                                              length.out = nrow(plantfate_daily_dataset))
  soil_water_potential_daily$YMD <- paste(soil_water_potential_daily$Year, soil_water_potential_daily$Month, soil_water_potential_daily$Day, sep = "-")
  plantfate_daily_dataset <- merge(plantfate_daily_dataset,
                                   soil_water_potential_daily[soil_water_potential_daily$Year < plantfate_daily_dataset$Year[nrow(plantfate_daily_dataset)],c("YMD", "HYY_EDDY233.GPP")],
                                   by = "YMD",
                                   all.x = T)

  plantfate_daily_dataset$Day <- substring(plantfate_daily_dataset$YMD, 9, 10)
  plantfate_daily_dataset$Month <- substring(plantfate_daily_dataset$YMD, 6, 7)

  daily_means <- plantfate_daily_dataset %>%
    group_by(Month, Day) %>%
    summarise(Mean_Temp = mean(Temp, na.rm = TRUE),
              Mean_VPD = mean(VPD, na.rm = TRUE),
              Mean_PPFD = mean(PPFD, na.rm = TRUE),
              Mean_PPFD_max = mean(PPFD_max, na.rm = TRUE))

  cols = c("Mean_Temp", "Mean_VPD", "Mean_PPFD", "Mean_PPFD_max")
  error_daily = Hyytiala_daily[,cols] - daily_means[,cols]

  # Preles - As preles data needs to be in different units etc.

  preles_daily_dataset <- daily_dataset
  preles_daily_dataset$VPD <- 10 * daily_dataset$VPD # kPa
  preles_daily_dataset$co2 <- 380
  preles_daily_dataset$fAPAR <- 0.7
  preles_daily_dataset$Precip <- daily_dataset$Precip # mm

  soil_water_potential_daily$YMD <- paste(soil_water_potential_daily$Year, soil_water_potential_daily$Month, soil_water_potential_daily$Day, sep = "-")
  preles_daily_dataset <- merge(preles_daily_dataset,
                                   soil_water_potential_daily[soil_water_potential_daily$Year < preles_daily_dataset$Year[nrow(preles_daily_dataset)],c("YMD", "HYY_EDDY233.GPP")],
                                   by = "YMD",
                                   all.x = T)

  preles_daily_dataset$Day <- substring(preles_daily_dataset$YMD, 9, 10)

  daily_means_preles <- preles_daily_dataset %>%
    group_by(Month, Day) %>%
    summarise(Mean_Temp = mean(Temp, na.rm = TRUE),
              Mean_VPD = mean(VPD, na.rm = TRUE),
              Mean_PAR = mean(PAR_preles, na.rm = TRUE),
              Mean_Precip = mean(Precip, na.rm = TRUE))

  cols = c("Mean_Temp", "Mean_VPD", "Mean_PAR", "Mean_Precip")
  error_daily_preles = Hyytiala_daily_prebas[,cols] - daily_means_preles[,cols]

  # SPP

  # Number of Columns: 13, units (daily_dataset name) [hyytiälä column]

  # CO2:                6, ppm
  # TotGlobal:          7, W m-2 (TotGlob)                      [Glob_mean]
  # TotPAR:             8, umol m-2 s-1 (PPFD)                  [PAR_mean] # TODO: check units
  # TAir:               9, 'C (Temp)                            [Mean_Temp]
  # Precip:             10, mm (Precip)                         [Mean_Precip]
  # Press:              11, kPa (Not here)                      [Mean_Press] ?? TODO: should I make a baseine as for the water potential?
  # VPD:                12, kPa (VPD but in hPa)                [Mean_VPD]
  # RH:                 13, % (RH 100 * decimal percentage)     [Mean_RH]
  # TSoil               'C (Temp_Soil_1)                        [Mean_STemp]
  # H2O                 ppth / mol m-3                          [Mean_SWC]

  spp_daily_dataset <- daily_dataset
  spp_daily_dataset$VPD <- 10 * daily_dataset$VPD
  spp_daily_dataset$RH <- 100 * daily_dataset$RH
  spp_daily_dataset$co2 <- 380 # TODO: real value?
  spp_daily_dataset$Press <- 1000 # TODO: real value?
  spp_daily_dataset$swvl1 <- daily_dataset$swvl1 # Although it say that m3 m-3 to mol m-3, the values look like they are the same as Hyytiälä without the correction

  spp_daily_dataset$Day <- substring(spp_daily_dataset$YMD, 9, 10)

  daily_means_spp <- spp_daily_dataset %>%
    group_by(Month, Day) %>%
    summarise(Mean_Glob = mean(TotGlob, na.rm = TRUE),
              Mean_PPFD = mean(PPFD, na.rm = TRUE),
              Mean_Temp = mean(Temp, na.rm = TRUE),
              Mean_Precip = mean(Precip, na.rm = TRUE),
              Mean_VPD = mean(10 * VPD, na.rm = TRUE),
              Mean_RH = mean(RH, na.rm = TRUE),
              Mean_STemp = mean(Temp_Soil_1, na.rm = TRUE),
              Mean_SWC = mean(swvl1, na.rm = TRUE))

  # TODO: Mean SWC is wrong units
  cols <- c("Mean_Glob", "Mean_PPFD", "Mean_Temp", "Mean_Precip", "Mean_VPD", "Mean_RH", "Mean_STemp", "Mean_SWC")
  error_spp_daily = Hyytiala_daily[,cols] - daily_means_spp[,cols]

  ### Bias correct

  # PlantFATE Monthly
  plantfate_monthy_dataset$Temp = plantfate_monthy_dataset$Temp + rep(error$Mean_Temp, length.out = nrow(plantfate_monthy_dataset))
  plantfate_monthy_dataset$VPD = plantfate_monthy_dataset$VPD + rep(error$Mean_VPD, length.out = nrow(plantfate_monthy_dataset))
  plantfate_monthy_dataset$PPFD = plantfate_monthy_dataset$PPFD + rep(error$Mean_PPFD, length.out = nrow(plantfate_monthy_dataset))
  plantfate_monthy_dataset$PPFD_max = plantfate_monthy_dataset$PPFD_max + rep(error$Mean_PPFD_max, length.out = nrow(plantfate_monthy_dataset))

  data.table::fwrite(plantfate_monthy_dataset[,c("Year", "Month", "Decimal_year", "Temp", "VPD", "PPFD", "PPFD_max", "SWP", "GPP")],
         file = file.path(path_test, "ERAS_Monthly.csv"))
  data.table::fwrite(plantfate_monthy_dataset[,c("Year", "Month", "Decimal_year", "Temp", "VPD", "PPFD", "PPFD_max", "SWP", "GPP")],
         file = file.path("./data/ERAS_Monthly.csv"))


  # PlantFATE Daily
  plantfate_daily_dataset$Temp = plantfate_daily_dataset$Temp + rep(error_daily$Mean_Temp, length.out = nrow(plantfate_daily_dataset))
  plantfate_daily_dataset$VPD = plantfate_daily_dataset$VPD + rep(error_daily$Mean_VPD, length.out = nrow(plantfate_daily_dataset))
  plantfate_daily_dataset$PPFD = plantfate_daily_dataset$PPFD + rep(error_daily$Mean_PPFD, length.out = nrow(plantfate_daily_dataset))
  plantfate_daily_dataset$PPFD_max = plantfate_daily_dataset$PPFD_max + rep(error_daily$Mean_PPFD_max, length.out = nrow(plantfate_daily_dataset))

  data.table::fwrite(plantfate_daily_dataset[,c("Year", "Month", "Decimal_year", "Temp", "VPD", "PPFD", "PPFD_max", "SWP")],
                     file = file.path(path_test, "ERAS_dataset_plantfate.csv"))
  data.table::fwrite(plantfate_daily_dataset[,c("Year", "Month", "Decimal_year", "Temp", "VPD", "PPFD", "PPFD_max", "SWP")],
                     file = file.path("./data/ERAS_dataset_plantfate.csv"))

  # Preles
  preles_daily_dataset$Temp <- preles_daily_dataset$Temp + rep(error_daily_preles$Mean_Temp, length.out = nrow(preles_daily_dataset))
  preles_daily_dataset$VPD <- preles_daily_dataset$VPD + rep(error_daily_preles$Mean_VPD, length.out = nrow(preles_daily_dataset))
  preles_daily_dataset$PAR <- preles_daily_dataset$PAR + rep(error_daily_preles$Mean_PAR, length.out = nrow(preles_daily_dataset))
  preles_daily_dataset$Temp <- preles_daily_dataset$Temp + rep(error_daily_preles$Mean_Temp, length.out = nrow(preles_daily_dataset))

  data.table::fwrite(preles_daily_dataset[,c("Temp", "VPD", "PAR", "Precip", "co2")],
                     file = file.path(path_test, "ERAS_dataset_preles.csv"))
  data.table::fwrite(preles_daily_dataset[,c("Temp", "VPD", "PAR", "Precip", "co2")],
                     file = file.path("./data/ERAS_dataset_preles.csv"))

  # SPP

  # CO2:                6, ppm
  # TotGlobal:          7, W m-2 (TotGlob)                      [Mean_Glob]
  # TotPAR:             8, umol m-2 s-1 (PPFD)                  [Mean_PPFD] # TODO: check units
  # TAir:               9, 'C (Temp)                            [Mean_Temp]
  # Precip:             10, mm (Precip)                         [Mean_Precip]
  # Press:              11, kPa (Not here)                      [Mean_Press] ?? TODO: should I make a baseine as for the water potential?
  # VPD:                12, kPa (VPD but in hPa)                [Mean_VPD]
  # RH:                 13, % (RH 100 * decimal percentage)     [Mean_RH]
  # TSoil               'C (Temp_Soil_1)                        [Mean_STemp]
  # H2O                 ppth / mol m-3 (swvl1)                  [Mean_SWC]

  spp_daily_dataset$TotGlob <- spp_daily_dataset$TotGlob + rep(error_spp_daily$Mean_Glob, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$PPFD <- spp_daily_dataset$PPFD + rep(error_spp_daily$Mean_PPFD, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$Temp <- spp_daily_dataset$Temp + rep(error_spp_daily$Mean_Temp, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$Precip <- spp_daily_dataset$Precip + rep(error_spp_daily$Mean_Precip, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$VPD <- spp_daily_dataset$VPD + rep(error_spp_daily$Mean_VPD, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$RH <- spp_daily_dataset$RH + rep(error_spp_daily$Mean_RH, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$Temp_Soil_1 <- spp_daily_dataset$Temp_Soil_1 + rep(error_spp_daily$Mean_STemp, length.out = nrow(spp_daily_dataset))
  spp_daily_dataset$swvl1 <- spp_daily_dataset$swvl1 + rep(error_spp_daily$Mean_SWC, length.out = nrow(spp_daily_dataset))

  data.table::fwrite(spp_daily_dataset[,c("Year", "Month", "Day", "co2", "TotGlob", "PPFD", "Temp", "Precip", "Press", "VPD", "RH", "Temp_Soil_1", "swvl1")],
                     file = file.path(path_test, "ERAS_dataset_spp.csv"))
  data.table::fwrite(spp_daily_dataset[,c("Year", "Month", "Day", "co2", "TotGlob", "PPFD", "Temp", "Precip", "Press", "VPD", "RH", "Temp_Soil_1", "swvl1")],
                     file = file.path("./data/ERAS_dataset_spp.csv"))

  # TODO: create them in the format that is needed to run the model?

  ####
  # Plot monthly and daily data
  ####
  plot_data(plantfate_monthy_dataset, "Monthly", monthly = TRUE)
  plot_data(plantfate_daily_dataset, "Daily", monthly = FALSE)

  return(paste("Your code is done - check the generated files in", path_test))
}
