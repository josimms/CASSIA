###
# Fast weather processing
###

rbinddatatable <- function(table_list) {
  # Define the common columns
  common_columns <- c("Year", "Month", "Day", "Hour", "Minute", "Second", "Date", "Monthly")

  # Process each table in the input list
  processed_tables <- lapply(table_list, function(table) {
    # Ensure the table is a data.table
    if (!is.data.table(table)) {
      table <- as.data.table(table)
    }

    # Reorder columns to put common columns first
    cols <- c(intersect(common_columns, names(table)),
              setdiff(names(table), common_columns))
    setcolorder(table, cols)

    return(table)
  })

  # Combine all processed tables
  result <- rbindlist(processed_tables, fill = TRUE)

  return(result)
}

###
# Downloading Data
###

downloading_data <- function(raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/",
                             year_start = 1995,
                             year_end = 2023) {

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

  # Doesn't return as it is downloading the files
}

###
# Downloading RAD data
###

downloading_rad_data <- function(
    raw.directory = "/home/josimms/Documents/Austria/Light Data/",
    year_start = 1995,
    year_end = 2026
) {
  http.origin <- "https://smear-backend-avaa-smear-prod.2.rahtiapp.fi/search/timeseries/csv?tablevariable=HYY_RAD."

  variables <- c("PAR", "Glob67", "maaPAR")

  # Build year pairs for chunked downloads (2-year chunks as in original)
  year1 <- seq(year_start, year_end, by = 2)
  year2 <- seq(year_start + 1, year_end, by = 2)
  if (length(year2) != length(year1)) year2 <- c(year2, year_end)

  for (variable in variables) {
    urls <- paste0(
      http.origin, variable,
      "&from=", year1, "-01-01T00%3A00%3A00.000",
      "&to=",   year2, "-12-31T23%3A59%3A59.999",
      "&quality=ANY&aggregation=NONE&interval=1"
    )
    dest_files <- paste0(raw.directory, "HYY_RAD.", variable, "_", year1, "to", year2, ".csv")

    mapply(download.file, urls, dest_files)
  }
}

###
# Read and plot RAD data: above-canopy PAR vs understory maaPAR (0.6 m)
###

plot_rad_comparison <- function(
    raw.directory = "/home/josimms/Documents/Austria/Light Data/",
    year_start    = NULL,
    year_end      = NULL
) {
  read_rad_variable <- function(dir, variable) {
    files <- list.files(dir, pattern = paste0("HYY_RAD\\.", variable, "_"), full.names = TRUE)
    if (length(files) == 0) stop("No files found for variable: ", variable)

    dt <- data.table::rbindlist(lapply(files, function(f) {
      d <- data.table::fread(f)
      data.table::setnames(d, gsub("HYY_RAD\\.", "", names(d)))
      d
    }), fill = TRUE)

    # Build a POSIXct timestamp for merging and plotting
    dt[, datetime := as.POSIXct(
      paste(Year, Month, Day, Hour, Minute, sep = "-"),
      format = "%Y-%m-%d-%H-%M", tz = "UTC"
    )]
    dt
  }

  par_dt    <- read_rad_variable(raw.directory, "PAR")
  maapar_dt <- read_rad_variable(raw.directory, "maaPAR")

  # Merge on the five time columns so we keep only matched rows
  merged <- merge(
    par_dt[,    .(Year, Month, Day, Hour, Minute, datetime, PAR)],
    maapar_dt[, .(Year, Month, Day, Hour, Minute, maaPAR)],
    by = c("Year", "Month", "Day", "Hour", "Minute")
  )

  # Optional year filter
  if (!is.null(year_start)) merged <- merged[Year >= year_start]
  if (!is.null(year_end))   merged <- merged[Year <= year_end]

  # Drop rows where both sensors are NA or negative (sensor off / night)
  merged <- merged[!is.na(PAR) & !is.na(maaPAR) & PAR >= 0 & maaPAR >= 0]

  # Daily means for cleaner time-series panels
  daily <- merged[, .(PAR    = mean(PAR,    na.rm = TRUE),
                      maaPAR = mean(maaPAR, na.rm = TRUE)),
                  by = .(Year, Month, Day)]
  daily[, date := as.Date(paste(Year, Month, Day, sep = "-"))]
  daily[, transmission := maaPAR / PAR]  # fraction of above-canopy light reaching 0.6 m
  data.table::setorder(daily, date)

  # --- Plots ---
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # 1. Time series – daily mean PAR above canopy
  plot(daily$date, daily$PAR,
       type = "l", col = "darkorange",
       xlab = "Date", ylab = expression(PAR ~ (mu*mol ~ m^{-2} ~ s^{-1})),
       main = "Above-canopy PAR (daily mean)")

  # 2. Time series – daily mean maaPAR at 0.6 m
  plot(daily$date, daily$maaPAR,
       type = "l", col = "steelblue",
       xlab = "Date", ylab = expression(maaPAR ~ (mu*mol ~ m^{-2} ~ s^{-1})),
       main = "Understory maaPAR at 0.6 m (daily mean)")

  # 3. Scatter: maaPAR vs PAR (all matched half-hourly / hourly points)
  plot(merged$PAR, merged$maaPAR,
       pch = 16, cex = 0.3, col = adjustcolor("forestgreen", alpha.f = 0.15),
       xlab = expression(PAR ~ (mu*mol ~ m^{-2} ~ s^{-1})),
       ylab = expression(maaPAR ~ (mu*mol ~ m^{-2} ~ s^{-1})),
       main = "maaPAR vs PAR (all observations)")
  abline(lm(maaPAR ~ PAR, data = merged), col = "red", lwd = 2)

  # 4. Daily canopy light-transmission fraction (maaPAR / PAR)
  # Exclude near-zero PAR to avoid division noise at dawn/dusk
  daily_filt <- daily[PAR > 10]
  plot(daily_filt$date, daily_filt$transmission,
       type = "l", col = "purple",
       ylim = c(0, 1),
       xlab = "Date", ylab = "maaPAR / PAR  (fraction)",
       main = "Canopy light transmission to 0.6 m")
  abline(h = median(daily_filt$transmission, na.rm = TRUE),
         col = "grey40", lty = 2)
  legend("topright", legend = paste0("median = ",
         round(median(daily_filt$transmission, na.rm = TRUE), 3)),
         bty = "n", text.col = "grey40")

  invisible(list(merged = merged, daily = daily))
}

###
# Memory efficient way of getting the dataframe with mean daily values - with max for a few variables
###
whole_weather_process <- function(save,
                                  raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/",
                                  output_directory = "/home/josimms/Documents/CASSIA/data") {

  variables <- c("RH672", "RH1250", "RHTd", "PAR", "CO2168", "T168", "T336", "Precip",
                 "tsoil_5", "tsoil_10", "wsoil_B1", "wsoil_B2", "Glob", "Glob67",
                 "Pamb336", "wpsoil_A", "wpsoil_B", "GPP")

  daily_result <- data.table::data.table(Date = character())
  hourly_result <- data.table::data.table(YMDH = character())

  for (variable in variables) {
    cat("Processing", variable, "\n")

    # Read and process files for each variable separately
    files <- list.files(raw.directory, pattern = variable, full.names = TRUE)

    variable_data <- data.table::rbindlist(lapply(files, function(file) {
      dt <- data.table::fread(file)

      # Simplify column names
      data.table::setnames(dt, gsub("HYY_META.|HYY_EDDY233.", "", names(dt)))

      # Add Date column
      dt[, Date := paste(Year, Month, Day, sep = "-")]
      dt[, YMDH := paste(Year, Month, Day, Hour, sep = "-")]

      return(dt)
    }))

    variable_data_YMDH <- variable_data[, .SD, .SDcols = c("YMDH", variable)]

    # Calculate daily statistics
    if (variable %in% c("Glob", "Glob67", "PAR")) {
      daily_stats <- variable_data[, .(
        mean = mean(get(variable), na.rm = TRUE),
        max = max(get(variable), na.rm = TRUE)
      ), by = Date]
      data.table::setnames(daily_stats, c("mean", "max"), paste0(variable, c("_mean", "_max")))
      hourly_stats <- variable_data[, .(
        mean = mean(get(variable), na.rm = TRUE),
        max = max(get(variable), na.rm = TRUE)
      ), by = YMDH]
      data.table::setnames(hourly_stats, c("mean", "max"), paste0(variable, c("_mean", "_max")))
    } else if (variable == "Precip") {
      daily_stats <- variable_data[, .(
        sum = sum(get(variable), na.rm = TRUE),
        mean = mean(get(variable), na.rm = TRUE)
      ), by = Date]
      data.table::setnames(daily_stats, c("sum", "mean"), paste0(variable, c("_sum", "_mean")))
      hourly_stats <- variable_data[, .(
        sum = sum(get(variable), na.rm = TRUE),
        mean = mean(get(variable), na.rm = TRUE)
      ), by = YMDH]
      data.table::setnames(hourly_stats, c("sum", "mean"), paste0(variable, c("_sum", "_mean")))
    } else {
      daily_stats <- variable_data[, .(mean = mean(get(variable), na.rm = TRUE)), by = Date]
      data.table::setnames(daily_stats, "mean", paste0(variable, "_mean"))
      hourly_stats <- variable_data[, .(mean = mean(get(variable), na.rm = TRUE)), by = YMDH]
      data.table::setnames(hourly_stats, "mean", paste0(variable, "_mean"))
    }

    # Merge with the result
    daily_result <- merge(daily_result, daily_stats, by = "Date", all = TRUE)
    hourly_result <- merge(hourly_result, hourly_stats, by = "YMDH", all = TRUE)

    # Clear memory
    rm(variable_data, daily_stats, hourly_stats)
    gc()
  }

  # Add Year, Month, and Day columns
  daily_result[, c("Year", "Month", "Day") := data.table::tstrsplit(Date, "-", fixed=TRUE)]
  hourly_result[, c("Year", "Month", "Day", "Hour") := data.table::tstrsplit(YMDH, "-", fixed=TRUE)]
  hourly_result[, Minute := 0]


  # Convert Year, Month, and Day to numeric for proper sorting
  daily_result[, ':='(Year = as.numeric(Year),
                      Month = as.numeric(Month),
                      Day = as.numeric(Day))]
  daily_result[, VPD := 10 * bigleaf::rH.to.VPD(0.01*RH672_mean, T336_mean)]

  hourly_result[, ':='(Year = as.numeric(Year),
                      Month = as.numeric(Month),
                      Day = as.numeric(Day))]
  hourly_result[, VPD := 10 * bigleaf::rH.to.VPD(0.01*RH672_mean, T336_mean)]

  # Add Monthly column
  daily_result[, Monthly := paste(Year, sprintf("%02d", Month), sep = "-")]

  # Sort the data by date
  data.table::setorder(daily_result, Year, Month, Day)

  # Reorder columns
  data.table::setcolorder(daily_result, c("Year", "Month", "Day", "Date", "Monthly"))

  # Save the result
  if (save) {
    data.table::fwrite(daily_result, file.path(output_directory, "daily_dataframe.csv"))
    data.table::fwrite(hourly_result, file.path(output_directory, "hourly_dataframe.csv"))

    cat("Processing complete. Results saved to", file.path(output_directory, "daily_dataframe.csv"), "\n")
  }

  return(daily_result)
}

###
# Gapfilling rows!
###

# Create a generalized gap-filling function
fill_missing_values <- function(data, pairs) {
  for (pair in pairs) {
    col1 <- pair[1]
    col2 <- pair[2]

    # Fill missing values based on mean-preserving substitution
    data <- data %>%
      mutate(
        !!col1 := ifelse(
          is.na(!!sym(col1)),
          !!sym(col2) - mean(!!sym(col2), na.rm = TRUE) + mean(!!sym(col1), na.rm = TRUE),
          !!sym(col1)
        ),
        !!col2 := ifelse(
          is.na(!!sym(col2)),
          !!sym(col1) - mean(!!sym(col1), na.rm = TRUE) + mean(!!sym(col2), na.rm = TRUE),
          !!sym(col2)
        )
      )
  }

  return(data)
}


###
# Final process!
###

raw_to_daily_monthly_hyytiala <- function(raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/",
                                          data.direct = "/home/josimms/Documents/CASSIA/data/") {
  library(dplyr)
  warning("This is originially ment for Joanna's own data processing. Thought that by making this function public it would help understadning of the weather processing. However it could be that paths used here don't work on your computer. If so contact Joanna!")

  ## Big files only run if really really necessary!f
  downloading_data()

  ### READING FILES AND MAKING MEAN, SUM AND MAX VALUES AS PER VARIABLE
  daily.dataframe <- whole_weather_process(T)
  # Note, as there are missing values in the day rather than create a multiplier for each day based on the missing data, just convert the mean

  plot(daily.dataframe$Precip_mean)

  ### saving dataframe raw!
  data.table::fwrite(daily.dataframe, paste0(data.direct, "daily.dataframe.hyytiala.csv"))

  ### Gapfill

  # Define column pairs for gap-filling
  column_pairs <- list(
    c("RH1250_mean", "RH672_mean"),
    c("RH1250_mean", "RHTd_mean"),
    c("T336_mean", "T168_mean"),
    c("tsoil_5_mean", "tsoil_10_mean"),
    c("wsoil_B1_mean", "wsoil_B2_mean"),
    c("wpsoil_A_mean", "wpsoil_B_mean"),
    c("PAR_mean", "Glob67_mean"),
    c("PAR_mean", "Glob_mean"),
    c("Glob_mean", "Glob67_mean"),
    c("Glob_max", "Glob67_max")
  )

  # Apply gap-filling
  all.gapfill <- fill_missing_values(daily.dataframe, column_pairs)
  # Only linear approx for variables that are used as at the moment doesn't work for every country
  variables.used <- zoo::na.approx(all.gapfill[,c("Year", "Month", "Day", "T336_mean", "VPD", "Glob_mean", "Glob_max", "wpsoil_A_mean", "PAR_mean", "Precip_sum", "Precip_mean", "tsoil_5_mean", "tsoil_10_mean", "wsoil_B1_mean")], maxgap = 7, na.rm = F)
  variables.used <- data.table::data.table(variables.used)

  # Make a water baseline temperarily

  NA_baselines <- variables.used %>%
    group_by(Year, Month, Day) %>%
    summarise(MB_mean = mean(wsoil_B1_mean, na.rm = T),
              Rain_mean = mean(Precip_mean, na.rm = T)) %>% # ??
    mutate_all(~replace(., is.infinite(.), NA)) %>%
    group_by(Day) %>%
    summarise(MB_mean = mean(MB_mean, na.rm = TRUE),
              Rain_mean = mean(Rain_mean, na.rm = TRUE))

  ### Different photosynthesis models

  ## phydro
  phydro <- variables.used[,c("Year", "Month", "T336_mean", "VPD", "PAR_mean", "Glob_max", "wpsoil_A_mean", "tsoil_5_mean", "tsoil_10_mean", "wsoil_B1_mean", "Precip_mean")]
  names(phydro)[c(3:7)] <- c("T", "VPD", "PAR", "PAR_max", "SWP")
  # Unit chnages
  # Temp = 'C, VPD = hPa to Pa, PPFD = umol m-2 s-1, SWP = kPa to -MPa
  phydro$VPD <- phydro$VPD * 100
  phydro$SWP <- 0.001 * phydro$SWP
  phydro$fAPAR <- 0.7
  phydro$Nitrogen <- 1
  phydro$P <- NA # Photosynthesis inputs!
  phydro$TSA = variables.used$tsoil_5_mean # Although not used in phydro
  phydro$TSB = variables.used$tsoil_10_mean # Although not used in phydro
  phydro$MB = variables.used$wsoil_B1_mean # Although not used in phydro
  phydro$Rain = variables.used$Precip_mean
  phydro$CO2 = 365
  phydro$PA = 101325

  phydro <- phydro[phydro$Year >= 2018]
  phydro$dates <- seq(as.Date("2018-01-01"), as.Date("2023-12-31"), by = "day")
  phydro$Decimal_year <- seq(phydro$Year[1],
                             phydro$Year[nrow(phydro)],
                             length.out = nrow(phydro))


  # Temporary MB baseline added for graphs
  MB_baseline_vector <- rep(NA_baselines$MB_mean, length.out = nrow(phydro))
  rain_baseline_vector <- rep(NA_baselines$Rain_mean, length.out = nrow(phydro))
  phydro$MB[is.na(phydro$MB)] <- NA_baseline_vector[is.na(phydro$MB)]
  phydro$Rain[is.na(phydro$Rain)] <- NA_baseline_vector[is.na(phydro$Rain)]

  # Plot the results
  par(mfrow = c(3, 3))
  plot(phydro$T, xlab = "", ylab = "", main = "Temperature, 'C")
  plot(phydro$VPD, xlab = "", ylab = "", main = "VPD, hPa")
  plot(phydro$PAR, xlab = "", ylab = "", main = "PPFD, umol m-2")
  plot(phydro$PAR_max, xlab = "", ylab = "", main = "max PPFD, max umol m-2")
  plot(phydro$SWP, xlab = "", ylab = "", main = "SWP, -MPa")
  plot(phydro$TSA, xlab = "", ylab = "", main = "Soil Temperature A, -MPa")
  plot(phydro$TSB, xlab = "", ylab = "", main = "Soil Temperature B, -MPa")
  plot(phydro$MB, xlab = "", ylab = "", main = "Soil Moisture, []")
  plot(phydro$Rain, xlab = "", ylab = "", main = "Precip mean, mm")

  # save file
  data.table::fwrite(phydro, paste0(data.direct, "phydro_smear_CASSIA_ready.csv"))

  ## preles
  preles <- variables.used[,c("Year", "PAR_mean", "T336_mean", "VPD")]
  # Temperature 'C, Precip sum mm,
  preles$PAR_mean <- 0.000001 * 24 * 60 * 60 * preles$PAR_mean # Sum of mmol m-2 s-1
  preles$VPD = 0.1 * preles$VPD # kPa
  preles$CO2 <- 380 # TODO: should I have real values for this
  preles$fAPAR <- 0.7
  preles$Nitrogen <- 1
  preles$P <- NA # Photosynthesis inputs!
  preles$TSA = variables.used$tsoil_5_mean # Although not used in preles
  preles$TSB = variables.used$tsoil_10_mean # Although not used in preles
  preles$MB = variables.used$wsoil_B1_mean # Although not used in preles
  preles$Rain = variables.used$Precip_sum # TODO: how much of a difference does the mean vs sum make for the precipitation response in the model?
  preles$CO2 = 365
  preles$PA = 101325
  names(preles)[2:3] <- c("PAR", "T")

  preles <- preles[preles$Year >= 2018]
  preles$dates <- seq(as.Date("2018-01-01"),
                      as.Date("2023-12-31"), by = "day")

  MB_baseline_vector <- rep(NA_baselines$MB_mean, length.out = nrow(preles))
  rain_baseline_vector <- rep(NA_baselines$Rain_mean, length.out = nrow(preles))
  preles$MB[is.na(preles$MB)] <- NA_baseline_vector[is.na(preles$MB)]
  preles$Rain[is.na(preles$Rain)] <- NA_baseline_vector[is.na(preles$Rain)]

  par(mfrow = c(3, 3))
  plot(preles$T, xlab = "", ylab = "", main = "Temperature, 'C")
  plot(preles$VPD, xlab = "", ylab = "", main = "VPD, kPa")
  plot(preles$PAR, xlab = "", ylab = "", main = "PPFD, sum mmol m-2")
  plot(preles$TSA, xlab = "", ylab = "", main = "Soil Temperature A, 'C")
  plot(preles$TSB, xlab = "", ylab = "", main = "Soil Temperature B, 'C")
  plot(preles$MB, xlab = "", ylab = "", main = "Soil Moisture, []")
  plot(preles$Rain, xlab = "", ylab = "", main = "Precip, sum mm")

  # save file
  data.table::fwrite(preles, paste0(data.direct, "preles_smear_CASSIA_ready.csv"))

  ### Import alternative data to check the similarity

  years <- c(2015, 2016, 2017, 2018)
  # Create a list to store the data frames
  weather_list <- lapply(years, function(year) {
    read.csv(file = paste0("/home/josimms/Documents/CASSIA/data/weather_original_", year, ".csv"),
             header = TRUE,
             sep = ",")
  })
  # Combine all data frames into a single data frame
  weather_original <- do.call(rbind, weather_list)

  ### Plot
  par(mfrow = c(2, 3))
  plot(weather_original$T, all.gapfill$T168_mean[all.gapfill$Year %in% years], main = "Temperature 'C", xlab = "", ylab = "")
  plot(weather_original$TSA, all.gapfill$tsoil_5_mean[all.gapfill$Year %in% years], main = "Temperature Soil A 'C", xlab = "", ylab = "")
  plot(weather_original$TSB, all.gapfill$tsoil_10_mean[all.gapfill$Year %in% years], main = "Temperature Soil A 'C", xlab = "", ylab = "")
  # TODO: soil moisture
  plot(weather_original$MB, all.gapfill$wsoil_B1_mean[all.gapfill$Year %in% years], main = "Soil Moisture ??", xlab = "", ylab = "")
  plot(weather_original$Rain, all.gapfill$Precip_mean[all.gapfill$Year %in% years], main = "Precipitation mm", xlab = "", ylab = "")

}
