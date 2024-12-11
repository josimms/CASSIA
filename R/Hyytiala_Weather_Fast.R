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
      daily_stats <- variable_data[, .(sum = sum(get(variable), na.rm = TRUE)), by = Date]
      data.table::setnames(daily_stats, "sum", paste0(variable, "_sum"))
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
# Final process!
###

raw_to_daily_monthly_hyytiala <- function(raw.directory = "/home/josimms/Documents/CASSIA_Calibration/Raw_Data/hyytiala_weather/",
                                          data.direct <- "/home/josimms/Documents/CASSIA/data/") {
  warning("This is originially ment for Joanna's own data processing. Thought that by making this function public it would help understadning of the weather processing. However it could be that paths used here don't work on your computer. If so contact Joanna!")

  ## Big files only run if really really necessary!
  downloading_data()

  ### READING FILES AND MAKING MEAN, SUM AND MAX VALUES AS PER VARIABLE
  daily.dataframe <- whole_weather_process(F)

  ### saving dataframe raw!
  fwrite(daily.dataframe, paste0(data.direct, "daily.dataframe.hyytiala.csv"))

  ### Gapfill TODO
  # daily.dataframe$RH672_mean <-


  ### TODO Create for different photosynthesis inputs

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
  plot(weather_original$T, daily.dataframe$T168_mean[daily.dataframe$Year %in% years], main = "Temperature 'C", xlab = "", ylab = "")
  plot(weather_original$TSA, daily.dataframe$tsoil_5_mean[daily.dataframe$Year %in% years], main = "Temperature Soil A 'C", xlab = "", ylab = "")
  plot(weather_original$TSB, daily.dataframe$tsoil_10_mean[daily.dataframe$Year %in% years], main = "Temperature Soil A 'C", xlab = "", ylab = "")
  # TODO: soil moisture
  plot(weather_original$MB, daily.dataframe$wsoil_B1_mean[daily.dataframe$Year %in% years], main = "Soil Moisture ??", xlab = "", ylab = "")
  plot(weather_original$Rain, daily.dataframe$Precip_sum[daily.dataframe$Year %in% years], main = "Precipitation mm", xlab = "", ylab = "")

}
