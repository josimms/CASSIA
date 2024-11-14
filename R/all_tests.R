###
# Ultimate test function
###

tests <- function() {
  # CASSIA extra parameters
  needle_mass_in = 4.467638

  # PRELES PARAMETERS
  # The values are taken from the Prebasso package!
  N_parameters = c(1/0.012, 0.0) # TODO: Fit the N parameters
  pPREL = c(413.0, 0.450, 0.118, 3.0, 0.748464, 12.74915, -3.566967, 18.4513, -0.136732,
            0.033942, 0.448975, 0.500, -0.364, 0.33271, 0.857291, 0.041781,
            0.474173, 0.278332, 1.5, 0.33, 4.824704, 0.0, 0.0, 180.0, 0.0, 0.0, 10.0,
            -999.9, -999.9, -999.9)

  ###
  # Xylem / cell enlargement plot
  ###
  n_rows = 100 # TODO: actual value
  max_ew_cells = 1.0 # TODO: actual value
  n_E_pot_old = 1.0 # TODO: actual value
  n_W_pot_old = 1.0 # TODO: actual value
  n_M_pot_old = 1.0 # TODO: actual value
  en_growth_vector = 1.0 # TODO: actual value
  tau_W_old = 1.0 # TODO: actual value
  carbon_daily_rate_ew = 1.0 # TODO: actual value
  carbon_daily_rate_lw = 1.0 # TODO: actual value
  g = 1.0 # TODO: actual value
  xylem_plot(t(parameters_p), common_p, t(sperling_p), sperling_extras,
             xylogensis_option, environmental_effect_xylogenesis, data_format,
             n_rows, max_ew_cells,
             n_E_pot_old, n_W_pot_old, n_M_pot_old, g, en_growth_vector,
             tau_W_old, carbon_daily_rate_ew, carbon_daily_rate_lw)
}

######
# Subfunctions for the ultermate test function
#####

# Function to read and combine weather data files for given years
read_and_combine_weather_data <- function(start_year, end_year, base_path = "./data/weather_original_") {
  weather_data_list <- lapply(start_year:end_year, function(year) {
    file_path <- paste0(base_path, year, ".csv")
    read.csv(file = file_path, header = TRUE, sep = ",")
  })
  do.call(rbind, weather_data_list)
}

process_weather_data <- function(using_spp_photosynthesis) {

  # Read and combine weather data from 2015 to 2018
  weather_original <- read_and_combine_weather_data(2015, 2018)
  weather_original$dates <- as.Date(strptime(paste(rep(2015:2018, times = c(365, 366, 365, 365)), weather_original$X), format = "%Y %j"))
  names(weather_original)[1] <- "date"

  # Add extra columns to the combined weather data
  extras <- data.frame(
    Nitrogen = rep(0.012, nrow(weather_original)),
    PAR = data_format[substring(data_format$Date, 1, 4) %in% 2015:2018, "PAR"],
    VPD = data_format[substring(data_format$Date, 1, 4) %in% 2015:2018, "VPD"],
    CO2 = data_format[substring(data_format$Date, 1, 4) %in% 2015:2018, "CO2"],
    fAPAR = rep(0.7, nrow(weather_original))
  )

  # Combine weather data with extra columns and remove specific rows
  weather_original <- cbind(weather_original, extras)
  weather_original <- weather_original

  # Set up dates
  dates_original <- seq(as.Date("2015-01-01"), as.Date("2018-12-31"), by = "day")
  dates_original <- dates_original
  dates <- dates_original
  start_year <- 2015
  end_year <- 2018

  # Adjust weather data and dates if not using spp photosynthesis
  if (!using_spp_photosynthesis) {
    weather_original <- data_format
    weather_original$Nitrogen <- rep(0.012, nrow(weather_original))
    weather_original$P <- rep(0.0, nrow(weather_original))

    adjust_temperature <- function(T_col, TSA_col, TSB_col) {
      TSA_col[is.na(TSA_col)] <- T_col[is.na(TSA_col)] - 7
      TSA_col[TSA_col < 0] <- 0
      TSB_col[is.na(TSB_col)] <- T_col[is.na(TSB_col)] - 7
      TSB_col[TSB_col < 0] <- 0
      list(TSA = TSA_col, TSB = TSB_col)
    }

    temp_adjusted <- adjust_temperature(weather_original$T, weather_original$TSA, weather_original$TSB)
    weather_original$TSA <- temp_adjusted$TSA
    weather_original$TSB <- temp_adjusted$TSB

    weather_original <- weather_original[-seq(from = 4, to = 19, by = 4) * 366, ]
    start_year <- 2005
    end_year <- 2023
    dates <- seq(as.Date("2005-01-01"), as.Date("2023-12-31"), by = "day")
    dates <- dates[-seq(from = 3, to = 18, by = 4) * 366 + 365]
  }

  # Return the processed weather data and dates
  return(list(weather_original = weather_original,
              dates = dates,
              dates_original = dates_original))
}

####
# Load data
####

load_data <- function(data_directory = "~/Documents/CASSIA_Calibration/Processed_Data/") {
  direct <- data_directory

  yu_data <- read.csv(paste0(direct, "yu.data.csv"))
  trenching_co2_fluxes <- read.delim(paste0(direct, "Trenching-CO2-fluxes-2013-2015.txt"), sep = "\t", header = FALSE)
  smearII_soil <- readxl::read_excel(paste0(direct, "smearII_soil.xlsx"))
  smearII_data <- readxl::read_excel(paste0(direct, "smearII_data.xlsx"))
  smearII <- read.csv(paste0(direct, "smearII.csv"))
  oyewole_2015_calibration_data <- read.csv(paste0(direct, "oyewole_2015_calibration_data.csv"))

  load(paste0(direct, "original.data.RData"))
  load(paste0(direct, "obs.vec.og.RData"))

  korhonen_mineral_N <- read.csv(paste0(direct, "korhonen_mineral_N.csv"))
  keskiarvot_kaittelyttain_vuot <- read.delim(paste0(direct, "keskiarvot_kasittelyttain_vuot_2013_2015.txt"))

  nitorgen_balance <- read.csv(paste0(direct, "ecosystem_balence_N.csv"))

  # Reading core drilling data files
  core_drilling_files <- list.files(direct, pattern = "Core Drilling")
  core_drilling_list <- lapply(paste0(direct, core_drilling_files), read.delim)

  # Print information about loaded data
  cat("Data loaded successfully.\n")
  cat("Loaded files:\n")
  cat(paste(core_drilling_files, collapse = "\n"), "\n")

  # Return a list of all loaded data frames and objects
  return(list(
    yu_data = yu_data,
    trenching_co2_fluxes = trenching_co2_fluxes,
    smearII_soil = smearII_soil,
    smearII_data = smearII_data,
    smearII = smearII,
    oyewole_2015_calibration_data = oyewole_2015_calibration_data,
    original_data = original.data,
    obs_vec = obs.vec.og,
    korhonen_mineral_N = korhonen_mineral_N,
    keskiarvot_kaittelyttain_vuot = keskiarvot_kaittelyttain_vuot,
    nitrogen_balance = nitorgen_balance,
    core_drilling_list = core_drilling_list
  ))
}

#####
# Parameters
#####

initialize_parameters <- function(calibration = FALSE, new_parameters = NULL) {
  # CASSIA extra parameters
  needle_mass_in <- 4.467638

  N_parameters <- c(1/0.012, 0.0) # TODO: Fit the N parameters
  pPREL <- c(413.0, 0.450, 0.118, 3.0, 0.748464, 12.74915, -3.566967, 18.4513, -0.136732,
             0.033942, 0.448975, 0.500, -0.364, 0.33271, 0.857291, 0.041781,
             0.474173, 0.278332, 1.5, 0.33, 4.824704, 0.0, 0.0, 180.0, 0.0, 0.0, 10.0,
             -999.9, -999.9, -999.9)

  parameters_test <- parameters_p
  parameters_test[c("lower_bound_needles", "lower_bound_phloem", "lower_bound_roots", "lower_bound_xylem_sh", "lower_bound_xylem_st"), 1] <- c(0.05, 0.13, 0.007, 0.009, 0.001)
  sperling_test <- sperling_p
  sperling_test[c("tau.s", "tau.t"), 1] <- c(3, 3)
  sperling_test[c("k_np", "k_pr", "k_pxsh", "k_pxst"), 1] <- c(100, 100, 100, 100)

  Throughfall <- 1

  if (calibration) {
    sperling_test <- sperling_p
    sperling_test[c(50:54, 35:39, 40:44, 45:49, 25, 26), 1] <- new_parameters[1:22]
    parameters_test <- parameters_p
    parameters_test[62:66, 1] <- new_parameters[23:27]
  }

  # TODO: All of the parameters should be checked
  parameters_R <- c(0.016906, # microbe_turnover (Preveen, 2013)
                    (0.001 * 2875) / 1881, # NC_in_root_opt (Heimisaari, 1995)
                    0.025, # NC_fungal_opt (Meyer, 2010)
                    1/28.73, # NC_microbe_opt (Heinonsalo, 2015)
                    0.5, # percentage_C_biomass (CASSIA)
                    # fungal
                    60, 40, 0.3, # N_limits: NH4, NO3, Norg
                    5, 5, 5, # N_k: NH4, NO3, Norg
                    0.3, 0.3, 0.3, # SWC_limits: NH4, NO3, Norg
                    # plant
                    60, 40, 0.3, # N_limits: NH4, NO3, Norg
                    5, 5, 5, # N_k: NH4, NO3, Norg
                    0.3, 0.3, 0.3, # SWC_limit: NH4, NO3, Norg
                    # microbes
                    60, 40, 0.3, 0.3, 10, # N_limits: NH4, NO3, Norg
                    5, 5, 5, 5, 5, # N_k: NH4, NO3, Norg
                    0.3, 0.3, 0.3, 0.3, 0.3, # SWC_limit: NH4, NO3, Norg
                    400, 400, 0.3, # C limits
                    10, # NH4_on_NO3
                    0.9, # optimal_root_fungal_biomass_ratio (TODO: Heinonsalo?)
                    1/625, # turnover_mantle, Meyer, 2010
                    1/50, # turnover_ERM, Meyer 2010
                    1/365, # turnover_roots, Meyer, 2010
                    1/625, # turnover_roots_mycorrhized Meyer, 2010
                    0.2, # turnover_fungal TODO: do I need this if the turnover is in Meyer?
                    1, # mantle_mass
                    1, # ERM_mass
                    0.2, # growth_C (Franklin, 2017) TODO: was as 0.3 or 0.03?
                    0.9, # growth_N (TODO)
                    1, # C_value_param_myco
                    1, # N_value_param_myco
                    1, # C_value_param_plant
                    1) # N_value_param_plant

  # Return parameters as a list
  return(list(
    needle_mass_in = needle_mass_in,
    N_parameters = N_parameters,
    pPREL = pPREL,
    parameters_test = parameters_test,
    sperling_test = sperling_test,
    Throughfall = Throughfall,
    parameters_R = parameters_R
  ))
}

#####
# Plot old CASSIA verison against new CASSIA version
#####

plot_comparison <- function(CASSIA_new_output, variables_new, Hyde_daily_original_plot, variables_original, soil_processes = FALSE) {
  par(mfrow = c(3, 3))

  dates = as.Date(strptime(paste(CASSIA_new_output$Growth$year, CASSIA_new_output$Growth$day), format = "%Y %j"))
  dates_original = as.Date(strptime(paste(Hyde_daily_original_plot$year, Hyde_daily_original_plot$day), format = "%Y %j"))

  for (var in 1:length(variables_new)) {
    if (var < length(variables_new)) {
      # Plot Outputs
      plot(dates, CASSIA_new_output$Growth[, variables_new[var]],
           main = "Outputs", xlab = "Date", ylab = gsub("_", " ", variables_new[var]), type = "l")
      lines(dates_original, Hyde_daily_original_plot[, variables_original[var]], col = "blue")

      # Plot New against Old
      plot(Hyde_daily_original_plot[, variables_original[var]],
           CASSIA_new_output$Growth[, variables_new[var]][-731],
           main = "New against old", xlab = "Original data", ylab = "New Data", col = "blue")
      abline(0, 1, col = "red")

      # Plot Residuals
      plot(dates_original, Hyde_daily_original_plot[, variables_original[var]] - CASSIA_new_output$Growth[, variables_new[var]][-731],
           main = "Residuals", xlab = "Date", ylab = "original - new output", col = "blue")

    } else {
      # Plot Outputs
      plot(dates, CASSIA_new_output$Preles[, variables_new[var]],
           main = "Outputs", xlab = "Date", ylab = variables_new[var], type = "l")
      lines(dates_original, Hyde_daily_original_plot$P * 1010 * 0.1, col = "blue")

      # Plot New against Old
      plot(Hyde_daily_original_plot$P * 1010 * 0.1,
           CASSIA_new_output$Preles[, variables_new[var]][-731],
           main = "New against old", xlab = "Original data", ylab = "New Data", col = "blue")
      abline(0, 1, col = "red")

      # Plot Residuals
      plot(dates_original, Hyde_daily_original_plot$P * 1010 * 0.1 - CASSIA_new_output$Preles[, variables_new[var]][-731],
           main = "Residuals", xlab = "Date", ylab = "original - new output", col = "blue")
      legend("topright", c("C++ Model", "Original R Model"), col = c("black", "blue"), bty = "n", lty = 1, cex = 0.75)
    }
  }
}

#####
# Weather variables
#####

plot_weather_variables <- function(weather_original, dates) {
  par(mfrow = c(3, 3))
  not_dates = !(names(weather_original) %in% c("dates", "Date", "date", "X"))
  for (clim in which(not_dates)) {
    plot(dates, weather_original[,clim],
         main = names(weather_original)[clim],
         ylab = names(weather_original)[clim],
         xlab = "Dates")
  }
}

#####
# Plot sugar start comparison
#####

plot_sugar_starch_comparison <- function(CASSIA_new_output, dates, original_data, yu_data) {
  par(mfrow = c(2, 1))

  # Plot Sugar
  sugar_total <- CASSIA_new_output$Sugar$sugar_needles +
    CASSIA_new_output$Sugar$sugar_phloem +
    CASSIA_new_output$Sugar$sugar_roots +
    CASSIA_new_output$Sugar$sugar_xylem_sh +
    CASSIA_new_output$Sugar$sugar_xylem_st

  plot(dates, sugar_total, main = "Sugar", col = "blue",
       type = "l", xlab = "Days of the Year", ylab = "Sugar, kg C",
       ylim = c(0, max(sugar_total, na.rm = TRUE)))
  abline(h = 0, lty = 2, col = "grey")
  lines(dates, sugar_total, col = "blue")
  lines(dates, CASSIA_new_output$Sugar$sugar_needles, col = "green")
  lines(dates, CASSIA_new_output$Sugar$sugar_phloem, col = "brown")
  lines(dates, CASSIA_new_output$Sugar$sugar_roots, col = "black")
  lines(dates, CASSIA_new_output$Sugar$sugar_xylem_sh, col = "red")
  lines(dates, CASSIA_new_output$Sugar$sugar_xylem_st, col = "orange")
  abline(h = 0.41, lty = 2, col = "blue")
  text(30, 0.43, "Expected Equilibrium", col = "blue", cex = 0.75)
  abline(h = sperling_p["SCb", "Hyde"], lty = 2, col = "pink")
  text(25, sperling_p["SCb", "Hyde"] + 0.02, "\"bloom\" threshold", col = "pink", cex = 0.75)
  points(as.Date(rownames(original_data)), original_data$needles_sugar, pch = "x", col = "green")
  points(as.Date(rownames(original_data)), original_data$phloem_sugar, pch = "x", col = "brown")
  points(as.Date(rownames(original_data)), original_data$roots_sugar, pch = "x", col = "black")
  points(as.Date(rownames(original_data)), original_data$xylem_sh_sugar, pch = "x", col = "red")
  points(as.Date(rownames(original_data)), original_data$xylem_st_sugar, pch = "x", col = "orange")
  points(as.Date(yu_data$Date), yu_data$sugar.needles, pch = "o", col = "green")
  points(as.Date(yu_data$Date), yu_data$sugar.phloem, pch = "o", col = "brown")
  points(as.Date(yu_data$Date), yu_data$sugar.root, pch = "o", col = "black")
  points(as.Date(yu_data$Date), yu_data$sugar.xylem.sh, pch = "o", col = "red")
  points(as.Date(yu_data$Date), yu_data$sugar.xylem.st, pch = "o", col = "orange")

  # Plot Starch
  starch_total <- CASSIA_new_output$Sugar$starch_needles +
    CASSIA_new_output$Sugar$starch_phloem +
    CASSIA_new_output$Sugar$starch_roots +
    CASSIA_new_output$Sugar$starch_xylem_sh +
    CASSIA_new_output$Sugar$starch_xylem_st

  plot(dates, starch_total, main = "Starch", col = "blue",
       type = "l", xlab = "Days of the Year", ylab = "Starch, kg C",
       ylim = c(0, max(starch_total, na.rm = TRUE)))
  abline(h = 0, lty = 2, col = "grey")
  lines(dates, starch_total, col = "blue")
  lines(dates, CASSIA_new_output$Sugar$starch_needles, col = "green")
  lines(dates, CASSIA_new_output$Sugar$starch_phloem, col = "brown")
  lines(dates, CASSIA_new_output$Sugar$starch_roots, col = "black")
  lines(dates, CASSIA_new_output$Sugar$starch_xylem_sh, col = "red")
  lines(dates, CASSIA_new_output$Sugar$starch_xylem_st, col = "orange")
  points(as.Date(rownames(original_data)), original_data$needles_starch, pch = "x", col = "green")
  points(as.Date(rownames(original_data)), original_data$phloem_starch, pch = "x", col = "brown")
  points(as.Date(rownames(original_data)), original_data$roots_starch, pch = "x", col = "black")
  points(as.Date(rownames(original_data)), original_data$xylem_sh_starch, pch = "x", col = "red")
  points(as.Date(rownames(original_data)), original_data$xylem_st_starch, pch = "x", col = "orange")
  points(as.Date(yu_data$Date), yu_data$starch.needles, pch = "o", col = "green")
  points(as.Date(yu_data$Date), yu_data$starch.phloem, pch = "o", col = "brown")
  points(as.Date(yu_data$Date), yu_data$starch.root, pch = "o", col = "black")
  points(as.Date(yu_data$Date), yu_data$starch.xylem.sh, pch = "o", col = "red")
  points(as.Date(yu_data$Date), yu_data$starch.xylem.st, pch = "o", col = "orange")

  par(mfrow = c(1, 1), new = TRUE)
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("left", c("Original Data", "Yu's Data", "Modelled Output"), bty = "n", pch = c("x", "o", "-"), title = "Data")
  legend("right", c("Needles", "Phloem", "Xylem, Shoot", "Xylem, Stem", "Roots"),
         bty = "n", col = c("green", "brown", "red", "orange", "black"), lty = 1, title = "Organs")
}

####
# Mycofon
####

plot_mycofon_data <- function(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, dates, loaded_data) {

  plot_fungal_data <- function(start, end, layout) {
    par(mfrow = layout)
    for (i in start:end) {
      # Calculate y-axis limits
      ylim <- if (any(!is.na(CASSIA_new_output_not_trenching$Fungal[,i]))) {
        range(c(CASSIA_new_output_not_trenching$Fungal[,i], CASSIA_new_output_trenching$Fungal[,i]), na.rm = TRUE)
      } else {
        c(0, 1)
      }

      # Plot data
      plot(dates, CASSIA_new_output_not_trenching$Fungal[,i],
           main = names(CASSIA_new_output_not_trenching$Fungal)[i],
           xlab = "Date", ylab = "", ylim = ylim, type = "l")
      lines(dates, CASSIA_new_output_trenching$Fungal[,i], col = "red")

      # Add reference data
      add_reference_data(names(CASSIA_new_output_not_trenching$Fungal)[i])

      # Check for all NAs
      if (all(is.na(CASSIA_new_output_not_trenching$Fungal[,i]))) {
        title(sub = paste0(names(CASSIA_new_output_not_trenching$Fungal)[i], " is all NAs!"))
        warning(paste0(names(CASSIA_new_output_not_trenching$Fungal)[i], " is all NAs!"))
      }

      # Add legend on last plot
      if (i == length(names(CASSIA_new_output_not_trenching$Fungal))) {
        legend("topleft", c("Control", "Trenched"),
               col = c("black", "red"), lty = 1, bty = "n")
      }
    }
  }

  add_reference_data <- function(var_name) {
    switch(var_name,
           "C_roots_NonStruct" = {
             points(as.Date(loaded_data$yu_data$Date), loaded_data$yu_data$sugar.root, pch = "o", col = "blue")
             points(as.Date(rownames(loaded_data$original_data)), loaded_data$original_data$roots_sugar, pch = "x", col = "blue")
           },
           "C_biomass" = {
             # TODO: check Ryhti data if not
             # TODO: Pauliina for the data that was originally used to calibrate CASSIA - root biomass
           },
           "C_fungal_NonStruct" = {
             # TODO: Add any specific plotting for this variable
           },
           "uptake_plant" = {
             abline(h = loaded_data$nitorgen_balance$Value_kg[loaded_data$nitorgen_balance$Item == "Uptake (aboveground)"], col = "blue", lty = 2)
             # TODO: why are the values different from the other file?
           },
           "uptake_fungal" = {
             abline(h = loaded_data$nitorgen_balance$Value_kg[loaded_data$nitorgen_balance$Item == "Uptake"], col = "blue", lty = 2)
             # TODO: why are the values different from the other file?
             # TODO: what is aboveground? What does this include?
             # TODO: what would be one tree?
           }
    )
  }

  # Plot first set of variables
  plot_fungal_data(1, 6, c(3, 2))

  # Plot second set of variables
  plot_fungal_data(7, length(names(CASSIA_new_output_not_trenching$Fungal)), c(3, 4))
}

#####
# Symphony
#####

plot_soil_data <- function(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, dates, loaded_data) {

  plot_soil_variables <- function(indices, mfrow = c(3, 3)) {
    par(mfrow = mfrow)
    for (i in indices) {
      ylim <- if (any(!is.na(CASSIA_new_output_not_trenching$Soil[,i]))) {
        range(c(CASSIA_new_output_not_trenching$Soil[,i], CASSIA_new_output_trenching$Soil[,i]), na.rm = TRUE)
      } else {
        c(0, 1)
      }

      plot(dates, CASSIA_new_output_not_trenching$Soil[,i],
           main = names(CASSIA_new_output_not_trenching$Soil)[i],
           xlab = "Date", ylab = "", ylim = ylim, type = "l")
      lines(dates, CASSIA_new_output_trenching$Soil[,i], col = "red")

      add_specific_data(names(CASSIA_new_output_not_trenching$Soil)[i])

      if (all(is.na(CASSIA_new_output_not_trenching$Soil[,i]))) {
        title(sub = paste0(names(CASSIA_new_output_not_trenching$Soil)[i], " is all NAs!"))
        warning(paste0(names(CASSIA_new_output_not_trenching$Soil)[i], " is all NAs!"))
      }

      if (i == tail(indices, 1)) {
        legend("topleft", c("Control", "Trenched"), col = c("black", "red"), lty = 1, bty = "n")
      }
    }
  }

  add_specific_data <- function(var_name) {
    switch(var_name,
           "C_decompose_SOM" = {
             abline(h = loaded_data$nitorgen_balance$Value_kg[loaded_data$nitorgen_balance$Item == "humus" & loaded_data$nitorgen_balance$Catagory == "Carbon"],
                    col = "blue", lty = 2)
           },
           "C_FOM_needles" = {
             points(loaded_data$karike_df_all$pvm, cumsum(replace(loaded_data$karike_df_all$neulanen, is.na(loaded_data$karike_df_all$neulanen), 0)),
                    col = "blue", pch = "+")
             title(sub = "Data: Needles in traps each year")
           },
           "C_FOM_woody" = {
             woody_sum <- rowSums(loaded_data$karike_df_all[c("oksa", "kuori", "käpy", "tikku")], na.rm = TRUE)
             points(loaded_data$karike_df_all$pvm, cumsum(woody_sum), col = "blue", pch = "+")
             title(sub = "Data: Branch + Bark + Cones + Sticks in traps each year")
           },
           "NH4" = {
             points(loaded_data$nitrogen$date, loaded_data$nitrogen$nh4, col = "blue") # TODO: units?
           },
           "NO3" = {
             points(loaded_data$nitrogen$date, loaded_data$nitrogen$no3, col = "blue") # TODO: units?
           }
    )
  }

  # Plot first set of variables
  plot_soil_variables(1:9)

  # Plot second set of variables
  plot_soil_variables(c(10:13, 19:20))
}

#####
# Uptake
#####

plot_nitrogen_uptake <- function(CASSIA_new_output_not_trenching, weather_original, oyewole_2015_calibration_data) {

  plot_uptake <- function(n_type) {
    par(mfrow = c(3, 4))

    ylim_ranges <- calculate_ylim(n_type)

    plot_against <- function(x_var, x_label) {
      plot_single(x_var, paste0("uptake_", n_type, "_fungal"), "fungal", x_label, ylim_ranges$Fungal)
      plot_single(x_var, paste0("uptake_", n_type, "_plant"), "plant", x_label, ylim_ranges$Plant)
      plot_single(x_var, paste0("uptake_", n_type, "_microbial_FOM"), "microbes_FOM", x_label, ylim_ranges$Plant)
      plot_single(x_var, paste0("uptake_", n_type, "_microbial_SOM"), "microbes_SOM", x_label, ylim_ranges$Plant)
    }

    plot_against(weather_original$TSB, "Temperature")
    plot_against(weather_original$MB, "Soil Moisture")

    n_range <- get_n_range(n_type)
    plot_against(n_range, n_type)

    if (n_type %in% c("NH4", "NO3", "Norg")) {
      add_oyewole_data(n_type)
    }
  }

  plot_single <- function(x, y, title_suffix, x_label, ylim) {
    y_data <- get_y_data(y)
    plot(x, y_data, main = paste0("uptake_", n_type, "_", title_suffix),
         ylab = "", xlab = x_label, ylim = ylim)
  }

  calculate_ylim <- function(n_type) {
    if (any(!is.na(CASSIA_new_output_not_trenching$Fungal[[paste0("uptake_", n_type, "_fungal")]]))) {
      list(
        Fungal = range(CASSIA_new_output_not_trenching$Fungal[[paste0("uptake_", n_type, "_fungal")]], na.rm = TRUE),
        Plant = range(CASSIA_new_output_not_trenching$Fungal[[paste0("uptake_", n_type, "_plant")]], na.rm = TRUE),
        Microbes = range(CASSIA_new_output_not_trenching$Preles[[paste0("uptake_", n_type, "_microbial")]], na.rm = TRUE)
      )
    } else {
      list(Fungal = c(0, 1), Plant = c(0, 1), Microbes = c(0, 1))
    }
  }

  get_n_range <- function(n_type) {
    if (any(!is.na(CASSIA_new_output_not_trenching$Soil[[n_type]]))) {
      CASSIA_new_output_not_trenching$Soil[[n_type]]
    } else {
      seq(0, 1, length.out = nrow(CASSIA_new_output_not_trenching$Soil))
    }
  }

  get_y_data <- function(y) {
    if (grepl("microbial", y)) {
      CASSIA_new_output_not_trenching$Preles[[y]]
    } else {
      CASSIA_new_output_not_trenching$Fungal[[y]]
    }
  }

  add_oyewole_data <- function(n_type) {
    n_comp <- switch(n_type,
                     "NH4" = "nh4",
                     "NO3" = "no3",
                     "Norg" = c("arginine", "glycine"))

    for (comp in n_comp) {
      points(oyewole_2015_calibration_data$concentration[oyewole_2015_calibration_data$n_comp == comp],
             oyewole_2015_calibration_data$control[oyewole_2015_calibration_data$n_comp == comp],
             pch = 17, col = "blue")
    }
  }

  # Plot for each nitrogen type
  for (n_type in c("NH4", "NO3", "Norg")) {
    plot_uptake(n_type)
  }
}

#####
# Respiration
#####

plot_respiration_data <- function(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, dates) {
  par(mfrow = c(3, 2))

  plot_respiration <- function(i) {
    range <- range(c(CASSIA_new_output_trenching$Respiration[,i],
                     CASSIA_new_output_not_trenching$Respiration[,i]), na.rm = TRUE)

    plot(dates, CASSIA_new_output_not_trenching$Respiration[,i],
         xlab = "Dates", ylab = "Respiration",
         ylim = range, xlim = range(dates),
         col = i, type = "l",
         main = clean_title(names(CASSIA_new_output_not_trenching$Respiration)[i]))

    if (i > 2) {
      lines(dates, CASSIA_new_output_trenching$Respiration[,i], col = 1, lty = 2)
    }
  }

  clean_title <- function(title) {
    gsub("_", " ", gsub("respiration_", "", title))
  }

  for (i in 1:ncol(CASSIA_new_output_not_trenching$Respiration)) {
    plot_respiration(i)
  }

  legend("topright", c("Control", "Trenching"), col = 1, lty = 1:2, bty = "n", cex = 0.75)
}

#####
# Total Ecosystem Respiration
#####

plot_total_ecosystem_respiration <- function(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, dates) {
  par(mfrow = c(1, 1))

  total_resp_control <- rowSums(CASSIA_new_output_not_trenching$Respiration)
  total_resp_trenching <- rowSums(CASSIA_new_output_not_trenching_trenching$Respiration)

  y_range <- range(c(total_resp_control, total_resp_trenching), na.rm = TRUE)

  plot(dates, total_resp_control,
       type = "l",
       main = "Total Ecosystem Respiration",
       xlab = "Dates",
       ylab = "Respiration",
       ylim = y_range,
       lwd = 2)

  lines(dates, total_resp_trenching, col = "red", lwd = 2)

  legend("topleft",
         c("Control", "Trenching"),
         col = c("black", "red"),
         lty = 1,
         lwd = 2,
         bty = "n")
}
