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
# Memory efficient way of getting the dataframe with mean daily values - with max for a few variables
###
whole_weather_process <- function(save,
                                  raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/",
                                  output_directory = "/home/josimms/Documents/CASSIA/data") {

  variables <- c("RH672", "RH1250", "RHTd", "PAR", "CO2168", "T168", "T336", "Precip",
                 "tsoil_5", "tsoil_10", "wsoil_B1", "wsoil_B2", "Glob", "Glob67",
                 "Pamb336", "wpsoil_A", "wpsoil_B", "GPP")

  daily_result <- data.table::data.table(Date = character())

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

      return(dt)
    }))

    # Calculate daily statistics
    if (variable %in% c("Glob", "Glob67", "PAR")) {
      daily_stats <- variable_data[, .(
        mean = mean(get(variable), na.rm = TRUE),
        max = max(get(variable), na.rm = TRUE)
      ), by = Date]
      data.table::setnames(daily_stats, c("mean", "max"), paste0(variable, c("_mean", "_max")))
    } else if (variable == "Precip") {
      daily_stats <- variable_data[, .(
        sum = sum(get(variable), na.rm = TRUE),
        mean = mean(get(variable), na.rm = TRUE)
      ), by = Date]
      data.table::setnames(daily_stats, c("sum", "mean"), paste0(variable, c("_sum", "_mean")))
    } else {
      daily_stats <- variable_data[, .(mean = mean(get(variable), na.rm = TRUE)), by = Date]
      data.table::setnames(daily_stats, "mean", paste0(variable, "_mean"))
    }

    # Merge with the result
    daily_result <- merge(daily_result, daily_stats, by = "Date", all = TRUE)

    # Clear memory
    rm(variable_data, daily_stats)
    gc()
  }

  # Add Year, Month, and Day columns
  daily_result[, c("Year", "Month", "Day") := data.table::tstrsplit(Date, "-", fixed=TRUE)]

  # Convert Year, Month, and Day to numeric for proper sorting
  daily_result[, ':='(Year = as.numeric(Year),
                      Month = as.numeric(Month),
                      Day = as.numeric(Day))]
  daily_result[, VPD := 10 * bigleaf::rH.to.VPD(0.01*RH672_mean, T336_mean)]

  # Add Monthly column
  daily_result[, Monthly := paste(Year, sprintf("%02d", Month), sep = "-")]

  # Sort the data by date
  data.table::setorder(daily_result, Year, Month, Day)

  # Reorder columns
  data.table::setcolorder(daily_result, c("Year", "Month", "Day", "Date", "Monthly"))

  # Save the result
  if (save) {
    data.table::fwrite(daily_result, file.path(output_directory, "daily_dataframe.csv"))

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
                                          data.direct <- "/home/josimms/Documents/CASSIA/data/") {
  library(dplyr)
  warning("This is originially ment for Joanna's own data processing. Thought that by making this function public it would help understadning of the weather processing. However it could be that paths used here don't work on your computer. If so contact Joanna!")

  ## Big files only run if really really necessary!f
  downloading_data()

  ### READING FILES AND MAKING MEAN, SUM AND MAX VALUES AS PER VARIABLE
  daily.dataframe <- whole_weather_process(F)
  # Note, as there are missing values in the day rather than create a multiplier for each day based on the missing data, just convert the mean

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

  ### Different photosynthesis models

  ## phydro
  # TODO: Dates!
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
