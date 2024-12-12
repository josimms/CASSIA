######
# Code for the functions in the ultimate test function
######

all_tests <- function(new_parameters, calibration, sperling_sugar_model, using_spp_photosynthesis, soil_processes, jussi_data) {
  ###
  # Settings and parameters
  ###
  settings_basic = list(
    storage_reset = TRUE,			# storage.reset<-TRUE=Same initial storage each year, storage.reset<-False, The storage on the last day of year X is  postponded to the first day of the year X+1
    storage_grows = FALSE,			# TRUE if the critical storage level increases with tree size.

    LN_estim = TRUE,				# LN depends on the GPP during previous july-august
    mN_varies = TRUE,				# needle mass (in maintenance respiration) is 2/3 of the total during period 1.10 - 31.5.

    LD_estim = TRUE,				# LD depends on the GPP during March-August
    sD_estim_T_count = FALSE,			# sD depends on the number of days when g in growing window - analogue to needles

    LH_estim = TRUE,

    trees_grow = FALSE,				# can be false if mature trees are modelled and not for a very long period
    growth_decreases = FALSE,			# the height and diameter growth (alfa_S and alfaD) decrease during the simulation
    needle_mass_grows = FALSE,		# Is needle mass dynamic i.e. the modelled growth is also respiring etc and following for some years? If true, note that root mass is related to needle mass

    phloem_trigger = FALSE,    # Phloem controls bud burst rather than whole tree sugar

    mycorrhiza = TRUE, 			# If allocation to mychorrhiza is taken into account
    root_as_Ding = TRUE,

    sperling_model = FALSE,       # Dynamic sugar model using Sperling's enzyme dynamics
    myco_model = FALSE,           # Joanna's mycomodel development!
    xylogenesis = FALSE,

    PRELES_GPP = FALSE,
    environment_effect_xylogenesis = FALSE,

    photosynthesis_as_input = TRUE,

    photoparameters = 3,
    temp_rise = FALSE,
    drought = FALSE,
    Rm_acclimation = TRUE,

    CASSIA_graphs = TRUE,
    tests = TRUE,

    etmodel = F,
    LOGFLAG = F
  )

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
  # Numerical tests of the functions in the model
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
  # Installing the model from github
  ###

  testing_remote = FALSE
  if (testing_remote) {
    devtools::install_github("josimms/CASSIA@adding_externals", force = TRUE)
  }

  ###
  # Weather Data Plots
  ###

  plot_weather_variables(processed_data$weather_original, processed_data$dates)

  ##################
  # PHOTOSYNTHESIS MODELS
  ##################

  ###
  # Preles CASSIA against PRELES original
  ###

  newwest_version = T
  if (!newwest_version) {
    devtools::install_github("josimms/Rprebasso", force = T)
  }

  processed_data <- process_weather_data(TRUE)

  preles_original <- Rprebasso::PRELES(processed_data$weather_original$PAR,
                                       processed_data$weather_original$T,
                                       processed_data$weather_original$VPD,
                                       processed_data$weather_original$Rain,
                                       processed_data$weather_original$CO2,
                                       processed_data$weather_original$fAPAR)

  preles_CASSIA <- preles_test(processed_data$weather_original)

  soil_processes = FALSE
  CASSIA_preles <- CASSIA_cpp(weather = processed_data$weather_original,
                              site = "Hyde",
                              pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                              parameters = parameters_all$parameters_test,
                              common = common_p,
                              ratios = ratios_p,
                              sperling = parameters_all$sperling_test,
                              needle_mass_in = parameters_all$needle_mass_in,
                              Throughfall = parameters_all$Throughfall,
                              photosynthesis_as_input = FALSE)

  CASSIA_preles_ecoevolutionary <- CASSIA_cpp(weather = processed_data$weather_original,
                                              site = "Hyde",
                                              pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                              parameters = parameters_all$parameters_test,
                                              common = common_p,
                                              ratios = ratios_p,
                                              sperling = parameters_all$sperling_test,
                                              needle_mass_in = parameters_all$needle_mass_in,
                                              Throughfall = parameters_all$Throughfall,
                                              photosynthesis_as_input = FALSE,
                                              ecoevolutionary = TRUE)

  CASSIA_preles_Xianglin <- CASSIA_cpp(weather = processed_data$weather_original,
                                       site = "Hyde",
                                       pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                       parameters = parameters_all$parameters_test,
                                       common = common_p,
                                       ratios = ratios_p,
                                       sperling = parameters_all$sperling_test,
                                       needle_mass_in = parameters_all$needle_mass_in,
                                       Throughfall = parameters_all$Throughfall,
                                       photosynthesis_as_input = FALSE,
                                       fAPAR_Tian = TRUE)

  CASSIA_preles_ecoevolutionary_Xianglin <- CASSIA_cpp(weather = processed_data$weather_original,
                                                       site = "Hyde",
                                                       pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                                       parameters = parameters_all$parameters_test,
                                                       common = common_p,
                                                       ratios = ratios_p,
                                                       sperling = parameters_all$sperling_test,
                                                       needle_mass_in = parameters_all$needle_mass_in,
                                                       Throughfall = parameters_all$Throughfall,
                                                       photosynthesis_as_input = FALSE,
                                                       ecoevolutionary = TRUE,
                                                       fAPAR_Tian = TRUE)

  par(mfrow = c(3, 2))
  ylim = c(0, max(c(CASSIA_preles_Xianglin$Preles$GPP, CASSIA_preles$Preles$GPP, preles_CASSIA$GPP)))
  plot(preles_original$GPP, ylab = "GPP", type = "l", ylim = ylim)
  lines(preles_CASSIA$GPP, col = "blue")
  lines(CASSIA_preles$Preles$GPP, col = "purple")
  lines(CASSIA_preles_ecoevolutionary$Preles$GPP, col = "orange")
  lines(CASSIA_preles_Xianglin$Preles$GPP, col = "green")
  lines(CASSIA_preles_ecoevolutionary_Xianglin$Preles$GPP, col = "yellow")

  ylim = c(0, max(c(preles_CASSIA$GPP, CASSIA_preles_Xianglin$Preles$GPP)))
  plot(preles_original$GPP, preles_CASSIA$GPP, xlab = "Original", ylab = "CASSIA", col = "blue", ylim = ylim)
  points(preles_original$GPP, CASSIA_preles$Preles$GPP, col = "purple")
  points(preles_original$GPP, CASSIA_preles_ecoevolutionary$Preles$GPP, col = "orange")
  points(preles_original$GPP, CASSIA_preles_Xianglin$Preles$GPP, col = "green")
  points(preles_original$GPP, CASSIA_preles_ecoevolutionary_Xianglin$Preles$GPP, col = "yellow")

  ylim = c(0, max(c(preles_CASSIA$ET, CASSIA_preles_Xianglin$Preles$ET)))
  plot(preles_original$ET, ylab = "ET", type = "l", ylim = ylim)
  lines(preles_CASSIA$ET, col = "blue")
  lines(CASSIA_preles$Preles$ET, col = "purple")
  lines(CASSIA_preles_ecoevolutionary$Preles$ET, col = "orange")
  lines(CASSIA_preles_Xianglin$Preles$ET, col = "green")
  lines(CASSIA_preles_ecoevolutionary_Xianglin$Preles$ET, col = "yellow")

  ylim = c(0, max(c(preles_CASSIA$ET, CASSIA_preles_Xianglin$Preles$ET)))
  plot(preles_original$ET, preles_CASSIA$ET, xlab = "Original", ylab = "CASSIA", col = "blue", ylim = ylim)
  points(preles_original$ET, CASSIA_preles$Preles$ET, col = "purple")
  points(preles_original$ET, CASSIA_preles_ecoevolutionary$Preles$ET, col = "orange")
  points(preles_original$ET, CASSIA_preles_Xianglin$Preles$ET, col = "green")
  points(preles_original$ET, CASSIA_preles_ecoevolutionary_Xianglin$Preles$ET, col = "yellow")

  plot(preles_original$SW, ylab = "Soil Water", type = "l")
  lines(preles_CASSIA$SoilWater, col = "blue")
  lines(CASSIA_preles$Preles$SoilWater, col = "purple")
  lines(CASSIA_preles_ecoevolutionary$Preles$SoilWater, col = "orange")
  lines(CASSIA_preles_Xianglin$Preles$SoilWater, col = "green")
  lines(CASSIA_preles_ecoevolutionary_Xianglin$Preles$SoilWater, col = "yellow")

  plot(preles_original$SW, preles_CASSIA$SoilWater, xlab = "Original", ylab = "CASSIA", col = "blue")
  points(preles_original$SW, CASSIA_preles$Preles$SoilWater, col = "purple")
  points(preles_original$SW, CASSIA_preles_ecoevolutionary$Preles$SoilWater, col = "orange")
  points(preles_original$SW, CASSIA_preles_Xianglin$Preles$SoilWater, col = "green")
  points(preles_original$SW, CASSIA_preles_ecoevolutionary_Xianglin$Preles$SoilWater, col = "yellow")
  legend("topleft", c("Original", "Preles Wrapper", "Preles in CASSIA", "Xianglin fAPAR"), col = c("black", "blue", "purple", "green"),
         pch = c(2, 2, 2), bty = "n")

  ##
  # Preles against original SPP CASSIA output
  ##

  processed_data <- process_weather_data(settings_basic$photosynthesis_as_input)
  soil_processes = FALSE

  CASSIA_new_output = CASSIA_cpp(weather = processed_data$weather_original,
                                 site = "Hyde",
                                 pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                 parameters = parameters_all$parameters_test,
                                 common = common_p,
                                 ratios = ratios_p,
                                 sperling = parameters_all$sperling_test,
                                 needle_mass_in = parameters_all$needle_mass_in,
                                 Throughfall = parameters_all$Throughfall,
                                 tests = settings_basic$tests)

  plot_comparison(CASSIA_new_output, variables_new,
                  Hyde_daily_original_plot, variables_original, soil_processes)

  plot_sugar_starch_comparison(CASSIA_new_output, processed_data$dates, loaded_data$original_data, loaded_data$yu_data)

  par(mfrow = c(2, 1))
  plot(CASSIA_new_output$Preles$GPP, col = "blue")
  points(CASSIA_preles$Preles$GPP, col = "black")

  plot_comparison(CASSIA_preles, variables_new,
                  Hyde_daily_original_plot, variables_original, soil_processes)

  ##
  # Phydro
  ##

  weather_Amazon <- read.delim("~/Documents/Austria/Plant-FATE/tests/data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE_new.csv", sep = ",")
  weather_Amazon$VPD <- 100 * weather_Amazon$VPD # Pa
  weather_Amazon$PAR <- weather_Amazon$PAR # umol m-2 s-1
  weather_Amazon$SWP <- - weather_Amazon$SWP
  # NOTE! Dates are absolutely not right! Just to see if the photosynthesis would work!
  weather_Amazon$dates <- seq(as.Date("1960-01-01"), as.Date("1960-07-10"), by = "day")
  names(weather_Amazon) <- c("Year", "Month", "Decimal_year", "T", "VPD", "PAR", "PAR_max", "SWP", "dates")
  weather_Amazon$fAPAR = 0.7
  weather_Amazon$Nitrogen = 1
  weather_Amazon$P = NA
  weather_Amazon$TSA = NA
  weather_Amazon$TSB = NA
  weather_Amazon$MB = NA
  weather_Amazon$Rain = NA
  weather_Amazon$CO2 = 365
  weather_Amazon$PA = 101325

  weather_Amazon_2013 <- weather_Amazon[weather_Amazon$Year > 2012,]

  # PAR [umol m-2 s-1], PAR Net radiation [W m-2], - MPa
  summary(weather_Amazon)

  CASSIA_phydro_amazon <- CASSIA_cpp(weather = weather_Amazon_2013,
                                     site = "Hyde",
                                     pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                     parameters = parameters_all$parameters_test,
                                     common = common_p,
                                     ratios = ratios_p,
                                     sperling = parameters_all$sperling_test,
                                     needle_mass_in = parameters_all$needle_mass_in,
                                     Throughfall = parameters_all$Throughfall,
                                     photosynthesis_as_input = FALSE,
                                     ecoevolutionary = TRUE)

  weather_ERAS_phydro <- read.delim("~/Documents/ERAS_dataset.csv", sep =",")
  weather_ERAS_phydro$VPD <- 0.1 * weather_ERAS_phydro$VPD # kPa
  weather_ERAS_phydro$dates <- seq(as.Date(paste(weather_ERAS_phydro$Year[1], weather_ERAS_phydro$Month[1], "01", sep = "-")),
                                   as.Date(paste(weather_ERAS_phydro$Year[nrow(weather_ERAS_phydro)], weather_ERAS_phydro$Month[nrow(weather_ERAS_phydro)], "31", sep = "-")),
                                   by = "day")
  names(weather_ERAS_phydro) <- c("Year", "Month", "Decimal_year", "T", "VPD", "PAR", "PAR_max", "SWP", "dates")
  weather_ERAS_phydro$fAPAR = 0.7
  weather_ERAS_phydro$Nitrogen = 1
  weather_ERAS_phydro$P = NA
  weather_ERAS_phydro$TSA = NA
  weather_ERAS_phydro$TSB = NA
  weather_ERAS_phydro$MB = NA
  weather_ERAS_phydro$Rain = NA
  weather_ERAS_phydro$CO2 = 421.38
  weather_ERAS_phydro$PA = 101325

  CASSIA_phydro_ecoevolutionary <- CASSIA_cpp(weather = weather_ERAS_phydro,
                                              site = "Hyde",
                                              pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                              parameters = parameters_all$parameters_test,
                                              common = common_p,
                                              ratios = ratios_p,
                                              sperling = parameters_all$sperling_test,
                                              needle_mass_in = parameters_all$needle_mass_in,
                                              Throughfall = parameters_all$Throughfall,
                                              photosynthesis_as_input = FALSE,
                                              ecoevolutionary = TRUE)

  plot(CASSIA_phydro_ecoevolutionary$Preles$GPP)

  ##
  # Running photosynthesis models with the same weather
  ##

  # TODO: redo with the ERA5 data!

  new_parameters <- rep(0.5, 27)  # Example new parameters
  calibration = F
  parameters_all <- initialize_parameters(calibration, new_parameters)

  data.direct = "./data/"
  phydro <- data.table::fread(paste0(data.direct, "phydro_smear_CASSIA_ready.csv"))
  preles <- data.table::fread(paste0(data.direct, "preles_smear_CASSIA_ready.csv"))

  smear_preles <- CASSIA_cpp(weather = preles,
                             site = "Hyde",
                             pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                             parameters = parameters_all$parameters_test,
                             common = common_p,
                             ratios = ratios_p,
                             sperling = parameters_all$sperling_test,
                             needle_mass_in = parameters_all$needle_mass_in,
                             Throughfall = parameters_all$Throughfall,
                             photosynthesis_as_input = FALSE)

  phydro_parameters_in = c(0.1008, 0.180496537959982, 5, 0.026263945805926, 0.011, 50,
                           0.5, -0.857817410110663, 4.1311874912949e17, 1, 2.45e-2, 2.0, 1.1, 0.1,
                           15, 10, 5, 0, 1, exp(-0.5 * 1.8), exp(-0.5 * 3.5), exp(-0.5 * 5.5))

  smear_phydro <- CASSIA_cpp(weather = phydro,
                             site = "Hyde",
                             pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                             parameters = parameters_all$parameters_test,
                             common = common_p,
                             ratios = ratios_p,
                             sperling = parameters_all$sperling_test,
                             needle_mass_in = parameters_all$needle_mass_in,
                             phydro = phydro_parameters_in,
                             Throughfall = parameters_all$Throughfall,
                             photosynthesis_as_input = FALSE,
                             ecoevolutionary = TRUE)

  par(mfrow = c(3, 1))
  plot(smear_preles$Preles$GPP, xlab = "Days since 2018-01-01", ylab = "Photosynthesis", col = "blue")
  points(smear_phydro$Preles$GPP, col = "green")

  plot(smear_preles$Growth$diameter_growth, xlab = "Days since 2018-01-01", ylab = "Diameter, kg C", col = "blue")
  points(smear_phydro$Growth$diameter_growth, col = "green")

  plot(smear_preles$Growth$root_growth, xlab = "Days since 2018-01-01", ylab = "Root, kg C", col = "blue")
  points(smear_phydro$Growth$root_growth, col = "green")

  ###
  # Sugar model
  ###


  ### SUGAR NO SOIL

  CASSIA_new_output = CASSIA_cpp(weather = processed_data$weather_original,
                                 site = "Hyde",
                                 pPREL = c(parameters_all$pPREL, parameters_all$N_parameters),
                                 parameters = parameters_all$parameters_test,
                                 common = common_p,
                                 ratios = ratios_p,
                                 sperling = parameters_all$sperling_test,
                                 needle_mass_in = parameters_all$needle_mass_in,
                                 Throughfall = parameters_all$Throughfall,
                                 sperling_model = T)

  plot_comparison(CASSIA_new_output, variables_new,
                  Hyde_daily_original_plot, variables_original, soil_processes)

  plot_sugar_starch_comparison(CASSIA_new_output, processed_data$dates, loaded_data$original_data, loaded_data$yu_data)

  ### SOIL

  ### Soil Processes
  soil_processes = TRUE
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
                                              soil = soil_processes)

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
                                            soil = soil_processes)

  plot_comparison(CASSIA_new_output_not_trenching, variables_new,
                  Hyde_daily_original_plot, variables_original, soil_processes)

  plot_comparison(CASSIA_new_output_trenching, variables_new,
                  Hyde_daily_original_plot, variables_original, soil_processes)

  plot_sugar_starch_comparison(CASSIA_new_output_not_trenching, processed_data$dates, loaded_data$original_data, loaded_data$yu_data)

  plot_sugar_starch_comparison(CASSIA_new_output_trenching, processed_data$dates, loaded_data$original_data, loaded_data$yu_data)

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

