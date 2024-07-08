######
# Code for the functions in the ultimate test function
######

all_tests <- function(new_parameters, calibration, sperling_sugar_model, using_spp_photosynthesis, soil_processes, jussi_data) {
  ###
  # Settings and parameters
  ###
  storage_rest = F
  storage_grows = F
  LH_estim = T
  LN_estim = T
  mN_varies = T
  LD_estim = T
  sD_estim_T_count = F
  trees_grow = F
  growth_decreases = F
  needle_mass_grows = F
  mycorrhiza = T
  root_as_Ding = T
  xylogensis_option = F
  environmental_effect_xylogenesis = F
  temp_rise = F
  drought = F
  Rm_acclimation = F
  trenching = T
  etmodel = F
  LOGFLAG = F

  ### WEATHER
  processed_data <- process_weather_data(using_spp_photosynthesis)
  dates_original <- seq(as.Date("2015-01-01"), as.Date("2018-12-31"), by = "day")[-c(366+365)]

  ### VALIDATION DATA
  data_directory <- "~/Documents/CASSIA_Calibration/Processed_Data/"
  loaded_data <- load_data(data_directory)
  loaded_data[["karike_df_all"]] <- karike_df_all # TODO: where is this imported?

  ### PARAMETERS
  new_parameters <- rep(0.5, 27)  # Example new parameters
  parameters <- initialize_parameters(calibration, new_parameters)

  ###
  # Running the functions
  ###

  soil_processes = T
  if (soil_processes) {
    trenching = FALSE
    CASSIA_new_output = CASSIA_soil(processed_data$start_year, processed_data$end_year, processed_data$weather_original, GPP_ref,
                                    c(parameters$pPREL, parameters$N_parameters), t(parameters$parameters_test), common_p, t(ratios_p), t(parameters$sperling_test), parameters$parameters_R,
                                    parameters$needle_mass_in,
                                    parameters$Throughfall,
                                    storage_rest, storage_grows,
                                    LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                                    trees_grow, growth_decreases, needle_mass_grows,
                                    mycorrhiza, root_as_Ding, sperling_sugar_model,
                                    xylogensis_option, environmental_effect_xylogenesis,
                                    temp_rise, drought, Rm_acclimation,
                                    using_spp_photosynthesis, 2100, TRUE,
                                    etmodel, LOGFLAG)
    trenching = TRUE
    CASSIA_new_output_trenching = CASSIA_soil(processed_data$start_year, processed_data$end_year, processed_data$weather_original, GPP_ref,
                                              c(parameters$pPREL, parameters$N_parameters), t(parameters$parameters_test), common_p, t(ratios_p), t(parameters$sperling_test), parameters$parameters_R,
                                              parameters$needle_mass_in,
                                              parameters$Throughfall,
                                              storage_rest, storage_grows,
                                              LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                                              trees_grow, growth_decreases, needle_mass_grows,
                                              mycorrhiza, root_as_Ding, sperling_sugar_model,
                                              xylogensis_option, environmental_effect_xylogenesis,
                                              temp_rise, drought, Rm_acclimation,
                                              using_spp_photosynthesis, 2015, TRUE,
                                              etmodel, LOGFLAG)
  } else {
    CASSIA_new_output = CASSIA_yearly(processed_data$start_year, processed_data$end_year, processed_data$weather_original, GPP_ref,
                                      c(parameters$pPREL, parameters$N_parameters), t(parameters$parameters_test), common_p, t(ratios_p), t(parameters$sperling_test),
                                      parameters$needle_mass_in,
                                      parameters$Throughfall,
                                      storage_rest, storage_grows,
                                      LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                                      trees_grow, growth_decreases, needle_mass_grows,
                                      mycorrhiza, root_as_Ding, sperling_sugar_model,
                                      xylogensis_option, environmental_effect_xylogenesis,
                                      temp_rise, drought, Rm_acclimation,
                                      using_spp_photosynthesis, TRUE,
                                      etmodel, LOGFLAG)
  }

  ###
  # Weather Data Plots
  ###

  plot_weather_variables(processed_data$weather_original, processed_data$dates)

  ###
  # Previous CASSIA Outpoints against new outputs
  ###
  variables_original <- c("bud", "wall_daily", "needle_daily", "root_daily", "height_daily", "Rg", "Rm", "P") # TODO: check if I want more, and that these are equivalent!
  variables_new <- c("bud_growth", "diameter_growth", "needle_growth", "root_growth", "height_growth", "respiration_growth", "respiration_maintenance", "GPP")

  Hyde_daily_original_plot <- Hyde_daily_original[1:(365+365+365+365),]

  # TODO: debug
  plot_comparison(CASSIA_new_output, variables_new, processed_data$dates, processed_data$dates_original,
                  Hyde_daily_original_plot, variables_original, processed_data$Photosynthesis_Reference, soil_processes)

  ###
  # Sugar
  ###
  plot_sugar_starch_comparison(CASSIA_new_output, processed_data$dates, loaded_data$original_data, loaded_data$yu_data)

  ###
  # Soil Processes
  ###

  if (soil_processes) {
    # Mycofon - Date
    plot_mycofon_data(CASSIA_new_output, CASSIA_new_output_trenching, processed_data$dates, loaded_data)

    # Symphony - Date
    plot_soil_data(CASSIA_new_output, CASSIA_new_output_trenching, dates, loaded_data)

    # Uptake - Environment
    plot_nitrogen_uptake(CASSIA_new_output, processed_data$weather_original, loaded_data$oyewole_2015_calibration_data)

    # Respiration soil breakdown
    plot_respiration_data(CASSIA_new_output, CASSIA_new_output_trenching, dates)
  }

  ###
  # Respiration
  ###

  plot_total_ecosystem_respiration(CASSIA_new_output, CASSIA_new_output_trenching, dates)

  return(CASSIA_new_output_trenching)
}

