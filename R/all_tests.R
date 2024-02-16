###
# Ultermate test function
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
# Code for the functions in the ultimate test function
######

all_tests <- function(new_parameters, calibration, sperling_sugar_model, using_spp_photosynthesis, soil_processes) {
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
  etmodel = F
  LOGFLAG = F

  # CASSIA extra parameters
  needle_mass_in = 4.467638

  N_parameters = c(1/0.012, 0.0) # TODO: Fit the N parameters
  pPREL = c(413.0, 0.450, 0.118, 3.0, 0.748464, 12.74915, -3.566967, 18.4513, -0.136732,
            0.033942, 0.448975, 0.500, -0.364, 0.33271, 0.857291, 0.041781,
            0.474173, 0.278332, 1.5, 0.33, 4.824704, 0.0, 0.0, 180.0, 0.0, 0.0, 10.0,
            -999.9, -999.9, -999.9)

  ###
  # Weather data processing
  ###

  weather_original_2015 = read.csv(file = "./data/weather_original_2015.csv", header = T, sep = ",")
  weather_original_2016 = read.csv(file = "./data/weather_original_2016.csv", header = T, sep = ",")
  weather_original_2017 = read.csv(file = "./data/weather_original_2017.csv", header = T, sep = ",")
  weather_original_2018 = read.csv(file = "./data/weather_original_2018.csv", header = T, sep = ",")
  # weather_original_2019 = data_format[,names(weather_original_2015)[-c(1, 3)]]
  weather_original = rbind(rbind(rbind(weather_original_2015, weather_original_2016), weather_original_2017), weather_original_2018)

  extras = data.frame(Nitrogen = rep(0.012, length = nrow(weather_original)),
                      PAR = data_format[substring(data_format$Date, 1, 4) %in% 2015:2018,c("PAR")],
                      VPD = data_format[substring(data_format$Date, 1, 4) %in% 2015:2018,c("VPD")],
                      CO2 = data_format[substring(data_format$Date, 1, 4) %in% 2015:2018,c("CO2")],
                      fAPAR = rep(0.7, length = nrow(weather_original)))
  weather_original <- cbind(weather_original, extras)

  start_year = 2015
  end_year = 2018
  dates_original = seq(as.Date("2015-01-01"), as.Date("2018-12-31"), by = "day")
  dates_original <- dates_original[-(365+366)]

  if (!using_spp_photosynthesis) {
    Photosynthesis_Reference <- weather_original$P
    weather_original <- data_format
    weather_original$Nitrogen <- rep(0.012, rep = nrow(weather_original))
    weather_original$P <- rep(0.0, rep = nrow(weather_original))
    weather_original$TSA[is.na(weather_original$TSA)] <- weather_original$T[is.na(weather_original$TSA)] - 5
    weather_original$TSA[weather_original$TSA < 0] <- 0
    weather_original$TSB[is.na(weather_original$TSB)] <- weather_original$T[is.na(weather_original$TSB)] - 5
    weather_original$TSB[weather_original$TSB < 0] <- 0
    weather_original <- weather_original[seq(from = 4, to = 19, by = 4)*-366,]
    start_year = 2005
    end_year = 2023
    dates = seq(as.Date("2005-01-01"), as.Date("2023-12-31"), by = "day")
    dates <- dates[-365+seq(from = 3, to = 18, by = 4)*-366]
  }

  direct <- "~/Documents/CASSIA_Calibration/"
  load(paste0(direct, "original.data.RData"))
  yu_data <- read.csv(paste0(direct, "yu.data.csv"))
  # TODO: is this the one with the right units?

  direct_raw_data <- "~/Documents/CASSIA_Calibration/Raw_Data/nitrogen_data/"
  nitorgen_balance <- read.csv(paste0(direct_raw_data, "ecosystem_balence_N.csv"))

  parameters_test <- parameters_p
  parameters_test[c("lower_bound_needles", "lower_bound_phloem", "lower_bound_roots", "lower_bound_xylem_sh", "lower_bound_xylem_st"),1] <- c(0.05, 0.13, 0.007, 0.009, 0.001)
  sperling_test <- sperling_p
  sperling_test[c("tau.s", "tau.t"),1] <- c(3, 3)
  sperling_test[c("k_np", "k_pr", "k_pxsh", "k_pxst"),1] <- c(100, 100, 100, 100)
  Throughfall = 1
  # SWinit, CWinit, SOGinit, Sinit

  if (calibration) {
    sperling_test <- sperling_p
    #sperling_test[c(50:54, 35:39, 40:44, 45:49, 25, 26),1] <- c(2.103433, 2.568670, 1.123682e+01, 7.150278e+00, 4.258494e-01, new_parameters[1:15], 2.567527, 5.467178)
    sperling_test[c(50:54, 35:39, 40:44, 45:49, 25, 26),1] <- new_parameters[1:22]
    parameters_test <- parameters_p
    parameters_test[62:66,1] <- new_parameters[23:27]
  }

  # TODO: All of the parameters should be checked
  parameters_R = c(0.016906, # microbe_turnover (Preveen, 2013)
                   (0.001*2875)/1881, # NC_in_root_opt (Heimisaari, 1995)
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
                   60, 40, 0.3, # N_limits: NH4, NO3, Norg
                   5, 5, 5, # N_k: NH4, NO3, Norg
                   0.3, 0.3, 0.3, # SWC_limit: NH4, NO3, Norg
                   60, 40, 0.3, # C limits
                   10, # NH4_on_NO3
                   0.9, # optimal_root_fungal_biomass_ratio (TODO: Heinonsalo?)
                   # TODO: 0.5 times just so I can get some good graphs for the presentation! Should be parametised
                   0, # 1/625, # turnover_mantle, Meyer, 2010
                   0, # 1/50, # turnover_ERM, Meyer 2010
                   0, # 1/365, # turnover_roots, Meyer, 2010
                   0, # 1/625, # turnover_roots_mycorrhized Meyer, 2010
                   0.2, # turnover_fungal TODO: do I need this if the turnover is in Meyer?
                   1, # mantle_mass
                   1, # ERM_mass
                   0.2, # growth_C (Franklin, 2017) TODO: was os 0.3 or 0.03?
                   0.9, # growth_N (TODO)
                   1, # C_value_param_myco
                   1, # N_value_param_myco
                   1, # C_value_param_plant
                   1) # N_value_param_plant

  if (soil_processes) {
    CASSIA_new_output = CASSIA_soil(start_year, end_year, weather_original, GPP_ref,
                                    c(pPREL, N_parameters), t(parameters_test), common_p, t(ratios_p), t(sperling_test), parameters_R,
                                    needle_mass_in,
                                    Throughfall,
                                    storage_rest, storage_grows,
                                    LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                                    trees_grow, growth_decreases, needle_mass_grows,
                                    mycorrhiza, root_as_Ding, sperling_sugar_model,
                                    xylogensis_option, environmental_effect_xylogenesis,
                                    temp_rise, drought, Rm_acclimation,
                                    using_spp_photosynthesis, TRUE,
                                    etmodel, LOGFLAG)
  } else {
    CASSIA_new_output = CASSIA_yearly(start_year, end_year, weather_original, GPP_ref,
                                      c(pPREL, N_parameters), t(parameters_test), common_p, t(ratios_p), t(sperling_test),
                                      needle_mass_in,
                                      Throughfall,
                                      storage_rest, storage_grows,
                                      LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                                      trees_grow, growth_decreases, needle_mass_grows,
                                      mycorrhiza, root_as_Ding, sperling_sugar_model,
                                      xylogensis_option, environmental_effect_xylogenesis,
                                      temp_rise, drought, Rm_acclimation,
                                      using_spp_photosynthesis, TRUE,
                                      etmodel, LOGFLAG)
  }


  variables_original <- c("bud", "wall_daily", "needle_daily", "root_daily", "height_daily", "Rg", "Rm", "P") # TODO: check if I want more, and that these are equivalent!
  variables_new <- c("bud_growth", "diameter_growth", "needle_growth", "root_growth", "height_growth", "respiration_growth", "respiration_maintenance", "GPP")

  Hyde_daily_original_plot <- Hyde_daily_original[1:(365+365+365+365),]

  ###
  # Previous Data
  ###

  par(mfrow = c(3, 3))
  for (var in 1:length(variables_new)) {
    if (var < length(variables_new)) {
      # TODO: Pauliina for the data that was originally used to calibrate CASSIA - root biomass in particular
      plot(dates, CASSIA_new_output[[1]][,c(variables_new[var])],
           main = "Outputs", xlab = "Date", ylab = variables_new[var], type = "l")
      lines(dates_original, Hyde_daily_original_plot[,c(variables_original[var])], col = "blue")
      plot(Hyde_daily_original_plot[,c(variables_original[var])][-731], CASSIA_new_output[[1]][,c(variables_new[var])][dates %in% dates_original], # TODO: the dates don't work here! Hence the errors
           main = "New against old", xlab = "Original data", ylab = "New Data", col = "blue")
      abline(0, 1, col = "red")
      plot(dates_original, Hyde_daily_original_plot[,c(variables_original[var])] - CASSIA_new_output[[1]][,c(variables_new[var])][dates %in% dates_original],
           main = "Residuals", xlab = "Date", ylab = "original - new output", col = "blue")
    } else {
      if (soil_processes) {
        photo_index = 5
      } else {
        photo_index = 3
      }
      plot(dates, CASSIA_new_output[[photo_index]][,c(variables_new[var])],
           main = "Outputs", xlab = "Date", ylab = variables_new[var], type = "l")
      lines(dates_original, Photosynthesis_Reference[-731], col = "blue")
      plot(Photosynthesis_Reference[-731][-length(Photosynthesis_Reference)+1], CASSIA_new_output[[photo_index]][,c(variables_new[var])][dates %in% dates_original], # TODO: the dates don't work here! Hence the errors
           main = "New against old", xlab = "Original data", ylab = "New Data")
      abline(0, 1, col = "red")
      plot(dates_original, Photosynthesis_Reference[-731] - CASSIA_new_output[[photo_index]][,c(variables_new[var])][dates %in% dates_original],
           main = "Residuals", xlab = "Date", ylab = "original - new output")
    }
  }

  ###
  # Sugar
  ###

  # TODO: make this into a function

  par(mfrow = c(2, 1))
  plot(dates, CASSIA_new_output$Sugar$sugar_needles +
         CASSIA_new_output$Sugar$sugar_phloem +
         CASSIA_new_output$Sugar$sugar_roots +
         CASSIA_new_output$Sugar$sugar_xylem_sh +
         CASSIA_new_output$Sugar$sugar_xylem_st, main = "Sugar", col = "blue",
       type = "l", xlab = "Days of the Year", ylab = "Sugar, kg C", ylim = c(0, max(CASSIA_new_output$Sugar$sugar_needles +
                                                                                      CASSIA_new_output$Sugar$sugar_phloem +
                                                                                      CASSIA_new_output$Sugar$sugar_roots +
                                                                                      CASSIA_new_output$Sugar$sugar_xylem_sh +
                                                                                      CASSIA_new_output$Sugar$sugar_xylem_st, na.rm = T)))
  abline(h = 0, lty = 2, col = "grey")
  lines(dates, CASSIA_new_output$Sugar$sugar_needles +
          CASSIA_new_output$Sugar$sugar_phloem +
          CASSIA_new_output$Sugar$sugar_roots +
          CASSIA_new_output$Sugar$sugar_xylem_sh +
          CASSIA_new_output$Sugar$sugar_xylem_st, col = "blue")
  lines(dates, CASSIA_new_output$Sugar$sugar_needles, col = "green")
  lines(dates, CASSIA_new_output$Sugar$sugar_phloem, col = "brown")
  lines(dates, CASSIA_new_output$Sugar$sugar_roots, col = "black")
  lines(dates, CASSIA_new_output$Sugar$sugar_xylem_sh, col = "red")
  lines(dates, CASSIA_new_output$Sugar$sugar_xylem_st, col = "orange")
  abline(h = 0.41, lty = 2, col = "blue")
  text(30, 0.43, "Expected Equilibrium", col = "blue", cex = 0.75)
  abline(h = sperling_p[c("SCb"),c("Hyde")], lty = 2, col = "pink")
  text(25, sperling_p[c("SCb"),c("Hyde")] + 0.02, "\"bloom\" threshold", col = "pink", cex = 0.75)
  points(as.Date(rownames(original.data)), original.data$needles_sugar, pch = "x", col = "green")
  points(as.Date(rownames(original.data)), original.data$phloem_sugar, pch = "x", col = "brown")
  points(as.Date(rownames(original.data)), original.data$roots_sugar, pch = "x", col = "black")
  points(as.Date(rownames(original.data)), original.data$xylem_sh_sugar, pch = "x", col = "red")
  points(as.Date(rownames(original.data)), original.data$xylem_st_sugar, pch = "x", col = "orange")
  points(as.Date(yu_data$Date), yu_data$sugar.needles, pch = "o", col = "green")
  points(as.Date(yu_data$Date), yu_data$sugar.phloem, pch = "o", col = "brown")
  points(as.Date(yu_data$Date), yu_data$sugar.root, pch = "o", col = "black")
  points(as.Date(yu_data$Date), yu_data$sugar.xylem.sh, pch = "o", col = "red")
  points(as.Date(yu_data$Date), yu_data$sugar.xylem.st, pch = "o", col = "orange")

  plot(dates, CASSIA_new_output$Sugar$starch_needles +
         CASSIA_new_output$Sugar$starch_phloem +
         CASSIA_new_output$Sugar$starch_roots +
         CASSIA_new_output$Sugar$starch_xylem_sh +
         CASSIA_new_output$Sugar$starch_xylem_st, main = "Starch", col = "blue",
       type = "l", xlab = "Days of the Year", ylab = "Starch, kg C", ylim = c(0, max(CASSIA_new_output$Sugar$starch_needles +
                                                                                CASSIA_new_output$Sugar$starch_phloem +
                                                                                CASSIA_new_output$Sugar$starch_roots +
                                                                                CASSIA_new_output$Sugar$starch_xylem_sh +
                                                                                CASSIA_new_output$Sugar$starch_xylem_st, na.rm = T)))
  abline(h = 0, lty = 2, col = "grey")
  lines(dates, CASSIA_new_output$Sugar$starch_needles +
          CASSIA_new_output$Sugar$starch_phloem +
          CASSIA_new_output$Sugar$starch_roots +
          CASSIA_new_output$Sugar$starch_xylem_sh +
          CASSIA_new_output$Sugar$starch_xylem_st, col = "blue")
  lines(dates, CASSIA_new_output$Sugar$starch_needles, col = "green")
  lines(dates, CASSIA_new_output$Sugar$starch_phloem, col = "brown")
  lines(dates, CASSIA_new_output$Sugar$starch_roots, col = "black")
  lines(dates, CASSIA_new_output$Sugar$starch_xylem_sh, col = "red")
  lines(dates, CASSIA_new_output$Sugar$starch_xylem_st, col = "orange")
  points(as.Date(rownames(original.data)), original.data$needles_starch, pch = "x", col = "green")
  points(as.Date(rownames(original.data)), original.data$phloem_starch, pch = "x", col = "brown")
  points(as.Date(rownames(original.data)), original.data$roots_starch, pch = "x", col = "black")
  points(as.Date(rownames(original.data)), original.data$xylem_sh_starch, pch = "x", col = "red")
  points(as.Date(rownames(original.data)), original.data$xylem_st_starch, pch = "x", col = "orange")
  points(as.Date(yu_data$Date), yu_data$starch.needles, pch = "o", col = "green")
  points(as.Date(yu_data$Date), yu_data$starch.phloem, pch = "o", col = "brown")
  points(as.Date(yu_data$Date), yu_data$starch.root, pch = "o", col = "black")
  points(as.Date(yu_data$Date), yu_data$starch.xylem.sh, pch = "o", col = "red")
  points(as.Date(yu_data$Date), yu_data$starch.xylem.st, pch = "o", col = "orange")

  par(mfrow = c(1, 1), new=TRUE)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("left", c("Original Data", "Yu's Data", "Modelled Output"), bty = "n", pch = c("x", "o", "-"), title = "Data")
  legend("right", c("Needles", "Phloem", "Xylem, Shoot", "Xylem, Stem", "Roots"),
         bty = "n", col = c("green", "brown", "red", "orange", "black"), lty = 1, title = "Organs")

  ###
  # Soil Processes
  ###
  if (soil_processes) {
    # Mycofon
    par(mfrow = c(3, 2))
    for (i in 1:6) {
      if (sum(is.na(CASSIA_new_output$Fungal[,i])) < nrow(CASSIA_new_output$Fungal)) {
        ylim = c(min(CASSIA_new_output$Fungal[,i], na.rm = T), max(CASSIA_new_output$Fungal[,i], na.rm = T))
      } else {
        ylim = c(0, 1)
      }
        # Plot
      plot(dates, CASSIA_new_output$Fungal[,i], main = names(CASSIA_new_output$Fungal)[i],
           xlab = "Date", ylab = "", ylim = ylim, type = "l")

      # Reference Data
      if (names(CASSIA_new_output$Fungal)[i] == "C_roots_NonStruct") {
        points(as.Date(yu_data$Date), yu_data$sugar.root, pch = "o", col = "blue")
        points(as.Date(rownames(original.data)), original.data$roots_sugar, pch = "x", col = "blue")
      } else if (names(CASSIA_new_output$Fungal)[i] == "C_biomass") {
        # TODO: check Ryhti data if not
        # TODO: Pauliina for the data that was originally used to calibrate CASSIA - root biomass
      }

      if (sum(is.na(CASSIA_new_output$Fungal[,i])) == nrow(CASSIA_new_output$Fungal)) {
        title(sub = paste0(names(CASSIA_new_output$Fungal)[i], " is all NAs!"))
        warning(paste0(names(CASSIA_new_output$Fungal)[i], " is all NAs!"))
      }
      # TODO: rest of the data should come from the Italian study
    }

    par(mfrow = c(3, 4))
    for (i in 7:length(names(CASSIA_new_output$Fungal))) {
      if (sum(is.na(CASSIA_new_output$Fungal[,i])) < nrow(CASSIA_new_output$Fungal)) {
        ylim = c(min(CASSIA_new_output$Fungal[,i], na.rm = T), max(CASSIA_new_output$Fungal[,i], na.rm = T))
      } else {
        ylim = c(0, 1)
      }

      # Plot
      plot(dates, CASSIA_new_output$Fungal[,i], main = names(CASSIA_new_output$Fungal)[i],
           xlab = "Date", ylab = "", ylim = ylim, type = "l")

      # Test data
      if (names(CASSIA_new_output$Fungal)[i] == "uptake_plant") {
        abline(h = nitorgen_balance$Value_kg[nitorgen_balance$Item == "Uptake (aboveground)"], col = "blue", lty = 2)
        # TODO: why are the values different from the other file?
      } else if (names(CASSIA_new_output$Fungal)[i] == "uptake_fungal") {
        abline(h = nitorgen_balance$Value_kg[nitorgen_balance$Item == "Uptake"], col = "blue", lty = 2)
        # TODO: why are the values different from the other file?
        # TODO: what is aboveground? What does this inlcude?
        # TODO: what would be one tree?
      }

      if (sum(is.na(CASSIA_new_output$Fungal[,i])) == nrow(CASSIA_new_output$Fungal)) {
        title(sub = paste0(names(CASSIA_new_output$Fungal)[i], " is all NAs!"))
        warning(paste0(names(CASSIA_new_output$Fungal)[i], " is all NAs!"))
      }
    }

    # Symphony
    par(mfrow = c(3, 3))
    for (i in 1:8) {
      if (sum(is.na(CASSIA_new_output$Soil[,i])) < nrow(CASSIA_new_output$Soil)) {
        ylim = c(min(CASSIA_new_output$Soil[,i], na.rm = T), max(CASSIA_new_output$Soil[,i], na.rm = T))
      } else {
        ylim = c(0, 1)
      }

      # TOO: units
      plot(dates, CASSIA_new_output$Soil[,i], main = names(CASSIA_new_output$Soil)[i],
           xlab = "Date", ylab = "", ylim = ylim, type = "l")

      if (names(CASSIA_new_output$Soil)[i] == "C_decompose_FOM") {
        # TODO: is this in the respiration paper?
      } else if (names(CASSIA_new_output$Soil)[i] == "C_decompose_SOM") {
        # TODO: is this in the respiration paper?
        abline(h = nitorgen_balance$Value_kg[nitorgen_balance$Item == "humus" & nitorgen_balance$Catagory == "Carbon"], col = "blue", lty = 2)
      } else if (names(CASSIA_new_output$Soil)[i] == "C_FOM_needles") {
        points(karike_df_all$pvm, cumsum(replace(karike_df_all$neulanen, is.na(karike_df_all$neulanen), 0)), col = "blue", pch = "+")
        title(sub = "Data: Needles in traps each year")
      } else if (names(CASSIA_new_output$Soil)[i] == "C_FOM_roots") {
        # TODO:
      } else if (names(CASSIA_new_output$Soil)[i] == "C_FOM_woody") {
        # TODO: this doesn't really make sense...?
        points(karike_df_all$pvm, cumsum(replace(karike_df_all$oksa + karike_df_all$kuori + karike_df_all$käpy + karike_df_all$tikku,
                                                 is.na(karike_df_all$oksa + karike_df_all$kuori + karike_df_all$käpy + karike_df_all$tikku),
                                                 0)),
               col = "blue", pch = "+")
        title(sub = "Data: Branch + Bark + Cones + Sticks in traps each year")
      } else if (names(CASSIA_new_output$Soil)[i] == "C_SOM") {
        # TODO:
      }

      if (sum(is.na(CASSIA_new_output$Soil[,i])) == nrow(CASSIA_new_output$Soil)) {
        title(sub = paste0(names(CASSIA_new_output$Soil)[i], " is all NAs!"))
        warning(paste0(names(CASSIA_new_output$Soil)[i], " is all NAs!"))
      }
    }

    par(mfrow = c(3, 3))
    for (i in c(9:12, 18:20)) {
      if (sum(is.na(CASSIA_new_output$Soil[,i])) < nrow(CASSIA_new_output$Soil)) {
        ylim = c(min(CASSIA_new_output$Soil[,i], na.rm = T), max(CASSIA_new_output$Soil[,i], na.rm = T))
      } else {
        ylim = c(0, 1)
      }

      plot(dates, CASSIA_new_output$Soil[,i], main = names(CASSIA_new_output$Soil)[i],
           xlab = "Date", ylab = "", ylim = ylim, type = "l")

      if (names(CASSIA_new_output$Soil)[i] == "NH4") {
        points(nitrogen$date, nitrogen$nh4, col = "blue") # TODO: units?
      } else if (names(CASSIA_new_output$Soil)[i] == "NO3") {
        points(nitrogen$date, nitrogen$no3, col = "blue") # TODO: units?
      }

      if (sum(is.na(CASSIA_new_output$Soil[,i])) == nrow(CASSIA_new_output$Soil)) {
        title(sub = paste0(names(CASSIA_new_output$Soil)[i], " is all NAs!"))
        warning(paste0(names(CASSIA_new_output$Soil)[i], " is all NAs!"))
      }
    }


    par(mfrow = c(3, 4))
    if (sum(is.na(CASSIA_new_output$Fungal$uptake_NH4_fungal)) < nrow(CASSIA_new_output$Fungal)) {
      ylim_Fungal = c(min(CASSIA_new_output$Fungal$uptake_NH4_fungal, na.rm = T), max(CASSIA_new_output$Fungal$uptake_NH4_fungal, na.rm = T))
      ylim_Plant = c(min(CASSIA_new_output$Fungal$uptake_NH4_plant, na.rm = T), max(CASSIA_new_output$Fungal$uptake_NH4_plant, na.rm = T))
      ylim_Microbes = c(min(CASSIA_new_output$Preles$uptake_NH4_microbial, na.rm = T), max(CASSIA_new_output$Fungal$uptake_NH4_microbial, na.rm = T))
    } else {
      ylim_Fungal = c(0, 1)
      ylim_Plant = c(0, 1)
      ylim_Microbes = c(0, 1)
    }

    plot(weather_original$TSB, CASSIA_new_output$Fungal$uptake_NH4_fungal, main = "uptake_NH4_fungal",
         ylab = "", xlab = "Temperature", ylim = ylim_Fungal)
    plot(weather_original$TSB, CASSIA_new_output$Fungal$uptake_NH4_plant, main = "uptake_NH4_plant",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$TSB, CASSIA_new_output$Preles$uptake_NH4_microbial_FOM, main = "uptake_NH4_microbes_FOM",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$TSB, CASSIA_new_output$Preles$uptake_NH4_microbial_SOM, main = "uptake_NH4_microbes_SOM",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Fungal$uptake_NH4_fungal, main = "uptake_NH4_fungal",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Fungal)
    plot(weather_original$MB, CASSIA_new_output$Fungal$uptake_NH4_plant, main = "uptake_NH4_plant",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Preles$uptake_NH4_microbial_FOM, main = "uptake_NH4_microbes_FOM",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Preles$uptake_NH4_microbial_SOM, main = "uptake_NH4_microbes_SOM",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    if (sum(is.na(CASSIA_new_output$Soil$NH4)) < nrow(CASSIA_new_output$Soil)) {
      NH4_range = CASSIA_new_output$Soil$NH4
    } else {
      NH4_range = seq(0, 1, length.out = nrow(CASSIA_new_output$Soil))
    }

    plot(NH4_range, CASSIA_new_output$Fungal$uptake_NH4_fungal, main = "uptake_NH4_fungal",
         ylab = "", xlab = "NH4", ylim = ylim_Fungal)
    plot(NH4_range, CASSIA_new_output$Fungal$uptake_NH4_plant, main = "uptake_NH4_plant",
         ylab = "", xlab = "NH4", ylim = ylim_Plant)
    points(oyewole_2015_calibration_data$concentration[oyewole_2015_calibration_data$n_comp == "nh4"],
           oyewole_2015_calibration_data$control[oyewole_2015_calibration_data$n_comp == "nh4"], pch = 17, col = "blue")
    plot(NH4_range, CASSIA_new_output$Preles$uptake_NH4_microbial_FOM, main = "uptake_NH4_microbes_FOM",
         ylab = "", xlab = "NH4", ylim = ylim_Plant)
    plot(NH4_range, CASSIA_new_output$Preles$uptake_NH4_microbial_SOM, main = "uptake_NH4_microbes_SOM",
         ylab = "", xlab = "NH4", ylim = ylim_Plant)

    # TODO: check this paper, was it tree or fungal uptake?
    if (sum(is.na(CASSIA_new_output$Fungal$uptake_NO3_fungal)) < nrow(CASSIA_new_output$Fungal)) {
      ylim_Fungal = c(min(CASSIA_new_output$Fungal$uptake_NO3_fungal, na.rm = T), max(CASSIA_new_output$Fungal$uptake_NO3_fungal, na.rm = T))
      ylim_Plant = c(min(CASSIA_new_output$Fungal$uptake_NO3_plant, na.rm = T), max(CASSIA_new_output$Fungal$uptake_NO3_plant, na.rm = T))
      ylim_Microbes = c(min(CASSIA_new_output$Preles$uptake_NO3_microbial, na.rm = T), max(CASSIA_new_output$Fungal$uptake_NO3_microbial, na.rm = T))
    } else {
      ylim_Fungal = c(0, 1)
      ylim_Plant = c(0, 1)
      ylim_Microbes = c(0, 1)
    }

    # TODO: plots

    plot(weather_original$TSB, CASSIA_new_output$Fungal$uptake_NO3_fungal, main = "uptake_NO3_fungal",
         ylab = "", xlab = "Temperature", ylim = ylim_Fungal)
    plot(weather_original$TSB, CASSIA_new_output$Fungal$uptake_NO3_plant, main = "uptake_NO3_plant",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$TSB, CASSIA_new_output$Preles$uptake_NO3_microbial_FOM, main = "uptake_NO3_microbes_FOM",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$TSB, CASSIA_new_output$Preles$uptake_NO3_microbial_SOM, main = "uptake_NO3_microbes_SOM",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Fungal$uptake_NO3_fungal, main = "uptake_NO3_fungal",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Fungal)
    plot(weather_original$MB, CASSIA_new_output$Fungal$uptake_NO3_plant, main = "uptake_NO3_plant",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Preles$uptake_NO3_microbial_FOM, main = "uptake_NO3_microbes_FOM",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Preles$uptake_NO3_microbial_SOM, main = "uptake_NO3_microbes_SOM",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    if (sum(is.na(CASSIA_new_output$Soil$NO3)) < nrow(CASSIA_new_output$Soil)) {
      NO3_range = CASSIA_new_output$Soil$NO3
    } else {
      NO3_range = seq(0, 1, length.out =  nrow(CASSIA_new_output$Soil))
    }
    plot(NO3_range, CASSIA_new_output$Fungal$uptake_NO3_fungal, main = "uptake_NO3_fungal",
         ylab = "", xlab = "NO3", ylim = ylim_Fungal)
    plot(NO3_range, CASSIA_new_output$Fungal$uptake_NO3_plant, main = "uptake_NO3_plant",
         ylab = "", xlab = "NO3", ylim = ylim_Plant)
    points(oyewole_2015_calibration_data$concentration[oyewole_2015_calibration_data$n_comp == "no3"],
           oyewole_2015_calibration_data$control[oyewole_2015_calibration_data$n_comp == "no3"], pch = 17, col = "blue")
    plot(NO3_range, CASSIA_new_output$Preles$uptake_NO3_microbial_FOM, main = "uptake_NO3_microbes_FOM",
         ylab = "", xlab = "NO3", ylim = ylim_Plant)
    plot(NO3_range, CASSIA_new_output$Preles$uptake_NO3_microbial_SOM, main = "uptake_NO3_microbes_SOM",
         ylab = "", xlab = "NO3", ylim = ylim_Plant)

    if (sum(is.na(CASSIA_new_output$Fungal$uptake_Norg_fungal)) < nrow(CASSIA_new_output$Fungal)) {
      ylim_Fungal = c(min(CASSIA_new_output$Fungal$uptake_Norg_fungal, na.rm = T), max(CASSIA_new_output$Fungal$uptake_Norg_fungal, na.rm = T))
      ylim_Plant = c(min(CASSIA_new_output$Fungal$uptake_Norg_plant, na.rm = T), max(CASSIA_new_output$Fungal$uptake_Norg_plant, na.rm = T))
      ylim_Microbes = c(min(CASSIA_new_output$Preles$uptake_Norg_microbial, na.rm = T), max(CASSIA_new_output$Fungal$uptake_Norg_microbial, na.rm = T))
    } else {
      ylim_Fungal = c(0, 1)
      ylim_Plant = c(0, 1)
      ylim_Microbes = c(0, 1)
    }

    plot(weather_original$TSB, CASSIA_new_output$Fungal$uptake_Norg_fungal, main = "uptake_Norg_fungal",
         ylab = "", xlab = "Temperature", ylim = ylim_Fungal)
    plot(weather_original$TSB, CASSIA_new_output$Fungal$uptake_Norg_plant, main = "uptake_Norg_plant",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$TSB, CASSIA_new_output$Preles$uptake_Norg_microbial_FOM, main = "uptake_Norg_microbes_FOM",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$TSB, CASSIA_new_output$Preles$uptake_Norg_microbial_SOM, main = "uptake_Norg_microbes_SOM",
         ylab = "", xlab = "Temperature", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Fungal$uptake_Norg_fungal, main = "uptake_Norg_fungal",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Fungal)
    plot(weather_original$MB, CASSIA_new_output$Fungal$uptake_Norg_plant, main = "uptake_Norg_plant",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Preles$uptake_Norg_microbial_FOM, main = "uptake_Norg_microbes_FOM",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    plot(weather_original$MB, CASSIA_new_output$Preles$uptake_Norg_microbial_SOM, main = "uptake_Norg_microbes_SOM",
         ylab = "", xlab = "Soil Moisture", ylim = ylim_Plant)
    if (sum(is.na(CASSIA_new_output$Soil$N_FOM)) < nrow(CASSIA_new_output$Soil)) {
      Norg_range = CASSIA_new_output$Soil$N_FOM
    } else {
      Norg_range = seq(0, 1, length.out =  nrow(CASSIA_new_output$Soil))
    }
    plot(Norg_range, CASSIA_new_output$Fungal$uptake_Norg_fungal, main = "uptake_Norg_fungal",
         ylab = "", xlab = "Norg", ylim = ylim_Fungal)
    plot(Norg_range, CASSIA_new_output$Fungal$uptake_Norg_plant, main = "uptake_Norg_plant",
         ylab = "", xlab = "Norg", ylim = ylim_Plant)
    points(oyewole_2015_calibration_data$concentration[oyewole_2015_calibration_data$n_comp == "arginine"],
           oyewole_2015_calibration_data$control[oyewole_2015_calibration_data$n_comp == "arginine"], pch = 17, col = "blue")
    points(oyewole_2015_calibration_data$concentration[oyewole_2015_calibration_data$n_comp == "glycine"],
           oyewole_2015_calibration_data$control[oyewole_2015_calibration_data$n_comp == "glycine"], pch = 17, col = "blue")
    plot(NO3_range, CASSIA_new_output$Preles$uptake_Norg_microbial_FOM, main = "uptake_Norg_microbes_FOM",
         ylab = "", xlab = "Norg", ylim = ylim_Plant)
    plot(NO3_range, CASSIA_new_output$Preles$uptake_Norg_microbial_SOM, main = "uptake_Norg_microbes_SOM",
         ylab = "", xlab = "Norg", ylim = ylim_Plant)
  }

  return(CASSIA_new_output)
}
