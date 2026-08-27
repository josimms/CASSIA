######
# Code for the functions in the ultimate test function
######

all_tests <- function(new_parameters, calibration, sperling_sugar_model, using_spp_photosynthesis, soil_processes, jussi_data) {

  ### WEATHER
  processed_data <- process_weather_data(settings_basic$photosynthesis_as_input)
  dates_original <- seq(as.Date("2015-01-01"), as.Date("2018-12-31"), by = "day")[-c(366+365)]

  ### VALIDATION DATA
  data_directory <- "~/Documents/CASSIA_Calibration/Processed_Data/"
  loaded_data <- load_data(data_directory)
  loaded_data[["karike_df_all"]] <- karike_df_all # TODO: where is this imported?

  ### PARAMETERS
  new_parameters <- rep(0.5, 27)  # Example new parameters
  calibration = F
  parameters_all <- initialize_parameters(calibration, new_parameters)

  ### PLOTTING INFORMATION
  variables_original <- c("bud", "wall_daily", "needle_daily", "root_daily", "height_daily", "Rg", "Rm", "P") # TODO: check if I want more, and that these are equivalent!
  variables_new <- c("bud_growth", "diameter_growth", "needle_growth", "root_growth", "height_growth", "respiration_growth", "respiration_maintenance", "GPP")

  Hyde_daily_original_plot <- Hyde_daily_original

  ###
  # Numerical tests of the functions in the model TODO: put in test that
  ###

  an.error.occured <- c()
  for (i in 1:nrow(processed_data$weather_original)) {
    weather_with_na_row = processed_data$weather_original
    weather_with_na_row[i,] = NA

    tryCatch({soil_processes = FALSE
    CASSIA_new_output_na_row = CASSIA_cpp(weather = weather_with_na_row,
                                          site = "Hyde",
                                          pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                          parameters = parameters_all$parameters_test,
                                          common = common_p,
                                          ratios = ratios_p,
                                          sperling = parameters_all$sperling_test,
                                          needle_mass_in = parameters_all$needle_mass_in,
                                          Throughfall = parameters_all$Throughfall)},
    error = function(e) {
      # Set error flag to TRUE if an error occurs
      an.error.occured[i] <<- TRUE
    }
    )
  }
  if (sum(an.error.occured) == nrow(processed_data$weather_original)) {
    print("Test passed NA in row")
  } else {
    stop("Test failed. Not all NAs in rows caught.")
  }

  an.error.occured <- c()
  for (j in c(2:7)) { ## Only checking for the original columns
    weather_with_na_column = processed_data$weather_original
    weather_with_na_column[,j] = NA

    tryCatch({CASSIA_new_output_na_col = CASSIA_cpp(weather = weather_with_na_column,
                                                    site = "Hyde",
                                                    pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                                    parameters = parameters_all$parameters_test,
                                                    common = common_p,
                                                    ratios = ratios_p,
                                                    sperling = parameters_all$sperling_test,
                                                    needle_mass_in = parameters_all$needle_mass_in,
                                                    Throughfall = parameters_all$Throughfall)},
             error = function(e) {
               # Set error flag to TRUE if an error occurs
               an.error.occured[j-1] <<- TRUE
             })
  }
  if (sum(an.error.occured) == 6) {
    print("Test passed NA in column")
  } else {
    stop("Test failed. Not all NAs in columns caught.")
  }

  photosynthesis_input_CASSIA_output <- CASSIA_cpp(weather = processed_data$weather_original,
                                                   site = "Hyde",
                                                   pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                                   parameters = parameters_all$parameters_test,
                                                   common = common_p,
                                                   ratios = ratios_p,
                                                   sperling = parameters_all$sperling_test,
                                                   needle_mass_in = parameters_all$needle_mass_in,
                                                   Throughfall = parameters_all$Throughfall,
                                                   photosynthesis_as_input = TRUE)

  if (sum(abs(processed_data$weather_original$P - photosynthesis_input_CASSIA_output$Preles$GPP) > rep(1e-17, nrow(photosynthesis_input_CASSIA_output$Preles))) > 1) {
    stop("Photosynthesis in not the same as photosynthesis out, when using_spp_photosynthesis is TRUE")
  }


  ###
  # Soil processes test
  ###

  ## PRELES

  ### Not Trenched
  CASSIA_new_output_no_tests = CASSIA_cpp(weather = processed_data$weather_original,
                                          site = "Hyde",
                                          pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                          parameters = parameters_all$parameters_test,
                                          common = common_p,
                                          ratios = ratios_p,
                                          sperling = parameters_all$sperling_test,
                                          parameters_R = parameters_all$parameters_R,
                                          needle_mass_in = parameters_all$needle_mass_in,
                                          Throughfall = parameters_all$Throughfall,
                                          soil = FALSE,
                                          photosynthesis_as_input = FALSE)

  ### Not Trenched
  CASSIA_new_output_not_trenching = CASSIA_cpp(weather = processed_data$weather_original,
                                               site = "Hyde",
                                               pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                               parameters = parameters_all$parameters_test,
                                               common = common_p,
                                               ratios = ratios_p,
                                               sperling = parameters_all$sperling_test,
                                               parameters_R = parameters_all$parameters_R,
                                               needle_mass_in = parameters_all$needle_mass_in,
                                               Throughfall = parameters_all$Throughfall,
                                               trenching_year = NA,
                                               soil = soil_processes,
                                               photosynthesis_as_input = FALSE)

  ### Trenched
  CASSIA_new_output_trenching = CASSIA_cpp(weather = processed_data$weather_original,
                                           site = "Hyde",
                                           pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                           parameters = parameters_all$parameters_test,
                                           common = common_p,
                                           ratios = ratios_p,
                                           sperling = parameters_all$sperling_test,
                                           parameters_R = parameters_all$parameters_R,
                                           needle_mass_in = parameters_all$needle_mass_in,
                                           Throughfall = parameters_all$Throughfall,
                                           trenching_year = 2015,
                                           soil = soil_processes,
                                           photosynthesis_as_input = FALSE)

  soil_against_yearly(CASSIA_new_output_no_tests, CASSIA_new_output_trenching, CASSIA_new_output_not_trenching)

  ## NO TESTS

  ### Not Trenched
  CASSIA_new_output_no_tests = CASSIA_cpp(weather = processed_data$weather_original,
                                          site = "Hyde",
                                          pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                          parameters = parameters_all$parameters_test,
                                          common = common_p,
                                          ratios = ratios_p,
                                          sperling = parameters_all$sperling_test,
                                          parameters_R = parameters_all$parameters_R,
                                          needle_mass_in = parameters_all$needle_mass_in,
                                          Throughfall = parameters_all$Throughfall,
                                          soil = FALSE,
                                          tests = TRUE)

  ### Not Trenched
  CASSIA_new_output_not_trenching = CASSIA_cpp(weather = processed_data$weather_original,
                                               site = "Hyde",
                                               pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                               parameters = parameters_all$parameters_test,
                                               common = common_p,
                                               ratios = ratios_p,
                                               sperling = parameters_all$sperling_test,
                                               parameters_R = parameters_all$parameters_R,
                                               needle_mass_in = parameters_all$needle_mass_in,
                                               Throughfall = parameters_all$Throughfall,
                                               trenching_year = NA,
                                               soil = soil_processes,
                                               tests = TRUE)

  ### Trenched
  CASSIA_new_output_trenching = CASSIA_cpp(weather = processed_data$weather_original,
                                           site = "Hyde",
                                           pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                           parameters = parameters_all$parameters_test,
                                           common = common_p,
                                           ratios = ratios_p,
                                           sperling = parameters_all$sperling_test,
                                           parameters_R = parameters_all$parameters_R,
                                           needle_mass_in = parameters_all$needle_mass_in,
                                           Throughfall = parameters_all$Throughfall,
                                           trenching_year = 2015,
                                           soil = soil_processes,
                                           tests = TRUE)

  soil_against_yearly(CASSIA_new_output_no_tests, CASSIA_new_output_trenching, CASSIA_new_output_not_trenching)

  ## TESTS

  ### Not Trenched
  CASSIA_new_output_no_tests = CASSIA_cpp(weather = processed_data$weather_original,
                                          site = "Hyde",
                                          pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                          parameters = parameters_all$parameters_test,
                                          common = common_p,
                                          ratios = ratios_p,
                                          sperling = parameters_all$sperling_test,
                                          parameters_R = parameters_all$parameters_R,
                                          needle_mass_in = parameters_all$needle_mass_in,
                                          Throughfall = parameters_all$Throughfall,
                                          soil = FALSE,
                                          tests = FALSE)

  ### Not Trenched
  CASSIA_new_output_not_trenching = CASSIA_cpp(weather = processed_data$weather_original,
                                               site = "Hyde",
                                               pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                               parameters = parameters_all$parameters_test,
                                               common = common_p,
                                               ratios = ratios_p,
                                               sperling = parameters_all$sperling_test,
                                               parameters_R = parameters_all$parameters_R,
                                               needle_mass_in = parameters_all$needle_mass_in,
                                               Throughfall = parameters_all$Throughfall,
                                               trenching_year = NA,
                                               soil = soil_processes,
                                               tests = FALSE)

  ### Trenched
  CASSIA_new_output_trenching = CASSIA_cpp(weather = processed_data$weather_original,
                                           site = "Hyde",
                                           pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                           parameters = parameters_all$parameters_test,
                                           common = common_p,
                                           ratios = ratios_p,
                                           sperling = parameters_all$sperling_test,
                                           parameters_R = parameters_all$parameters_R,
                                           needle_mass_in = parameters_all$needle_mass_in,
                                           Throughfall = parameters_all$Throughfall,
                                           trenching_year = 2015,
                                           soil = soil_processes,
                                           tests = FALSE)

  soil_against_yearly(CASSIA_new_output_no_tests, CASSIA_new_output_trenching, CASSIA_new_output_not_trenching)

  ### Extra Soil Processes

  if (soil_processes) {
    # Decision values
    # TODO: write the decision graphs

    # Mycofon - Date
    plot_mycofon_data(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, processed_data$dates, loaded_data)

    # Symphony - Date
    plot_soil_data(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, processed_data$weather_original$date, loaded_data)

    # Uptake - Environment
    plot_nitrogen_uptake(CASSIA_new_output_not_trenching, processed_data$weather_original, loaded_data$oyewole_2015_calibration_data)

    # Respiration soil breakdown
    plot_respiration_data(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, dates)
  }

  ###
  # Respiration
  ###

  plot_total_ecosystem_respiration(CASSIA_new_output_not_trenching, CASSIA_new_output_trenching, dates)

  # TODO: add a xylogenesis function here

  # TODO: add the water dependencies

  return(CASSIA_new_output_trenching)
}

