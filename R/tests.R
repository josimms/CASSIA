###
# Ultermate test function
###

#' @export
tests <- function() {
  start_year = 2010
  end_year = 2019

  # INITIAL SETTINGS OF THE ORIGINAL MODEL
  storage_rest = T
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
  sperling_sugar_model = F
  xylogensis_option = F
  environmental_effect_xylogenesis = F
  temp_rise = F
  drought = F
  Rm_acclimation = F
  etmodel = F
  LOGFLAG = F

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
  # Respiration tests
  ###
  respiration_plot()

  ###
  # Repola test
  ###
  repola_plot()

  ###
  # PRELES test
  ###
  PRELES_plot(data_format, N_parameters)

  ###
  # Growth plot - TODO: do this!
  ###
  # TODO: make this!
  growth_plot()

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

  ###
  # Outputs between R code and the cpp code WITHOUT the sugar model
  ###

  # TODO: add nitrogen here!

  test_against_original_data()

  ###
  # Outputs between R code and the cpp code WITH the sugar model
  # TODO: needs to be calibrated first
  ###

  sugar_plot(CASSIA_out, site, SCb = sperling_p[c("SCb"), c("Hyde")], cpp = F)

}

######
# Code for the functions in the ultimate test function
######

### Weather tests

tests <- function(weather) {
  ###
  # Test to see if the weather files are okay
  ###
  if (sum(names(weather) %in% c("date", "Date", "T", "P", "TSA", "TSB", "MB", "Rain", "PAR", "VPD", "fAPAR"))) {stop("Incomplete weather data - incorrect variables, or named incorrectly")}
}

### Test against original data

test_against_original_data <- function(new_parameters, calibration, sperling_sugar_model, using_spp_photosynthesis) {
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

  weather_original_2015 = read.csv(file = "./data/weather_original_2015.csv", header = T, sep = ",")
  weather_original_2016 = read.csv(file = "./data/weather_original_2016.csv", header = T, sep = ",")
  weather_original_2017 = read.csv(file = "./data/weather_original_2017.csv", header = T, sep = ",")
  weather_original_2018 = read.csv(file = "./data/weather_original_2018.csv", header = T, sep = ",")
  weather_original_2019 = data_format[,names(weather_original_2015)[-c(1, 3)]]
  weather_original = rbind(rbind(rbind(weather_original_2015, weather_original_2016), weather_original_2017), weather_original_2018)


  extras = data.frame(Nitrogen = rep(0.012, length = nrow(weather_original)),
                      PAR = data_format[substring(data_format$Date, 1, 4) %in% 2015:2019,c("PAR")],
                      VPD = data_format[substring(data_format$Date, 1, 4) %in% 2015:2019,c("VPD")],
                      CO2 = data_format[substring(data_format$Date, 1, 4) %in% 2015:2019,c("CO2")],
                      fAPAR = rep(0.7, length = nrow(weather_original)))
  weather_original <- cbind(weather_original, extras)
  weather_original <- weather_original[-c(365+365),]

  direct <- "~/Documents/CASSIA_Calibration/"
  load(paste0(direct, "original.data.RData"))

  yu_data <- read.csv(paste0(direct, "yu.data.csv"))

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

  CASSIA_new_output = CASSIA_yearly(2015, 2018, weather_original, GPP_ref,
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

  variables_original <- c("bud", "wall_daily", "needle_daily", "root_daily", "height_daily", "Rg", "Rm", "P") # TODO: check if I want more, and that these are equivalent!
  variables_new <- c("bud_growth", "diameter_growth", "needle_growth", "root_growth", "height_growth", "respiration_growth", "respiration_maintenance", "GPP")

  dates = seq(as.Date("2015-01-01"), as.Date("2019-12-31"), by = "day")
  dates = dates[-c(365+366)]

  Hyde_daily_original_plot <- Hyde_daily_original[1:(365+365+365+365+365),]

  par(mfrow = c(3, 3))
  for (var in 1:length(variables_new)) {
    if (var < length(variables_new)) {
      plot(dates, Hyde_daily_original_plot[,c(variables_original[var])],
           main = "Outputs", xlab = "Date", ylab = variables_new[var], type = "l")
      lines(dates, CASSIA_new_output[[1]][,c(variables_new[var])], col = "blue")
      plot(Hyde_daily_original_plot[,c(variables_original[var])], CASSIA_new_output[[1]][,c(variables_new[var])],
           main = "New against old", xlab = "Original data", ylab = "New Data")
      abline(0, 1, col = "red")
      plot(dates, Hyde_daily_original_plot[,c(variables_original[var])] - CASSIA_new_output[[1]][,c(variables_new[var])],
           main = "Residuals", xlab = "Date", ylab = "original - new output")
    } else {
      plot(dates, weather_original$P,
           main = "Outputs", xlab = "Date", ylab = variables_new[var], type = "l")
      lines(dates, CASSIA_new_output[[3]][,c(variables_new[var])], col = "blue")
      plot(weather_original$P, CASSIA_new_output[[3]][,c(variables_new[var])],
           main = "New against old", xlab = "Original data", ylab = "New Data")
      abline(0, 1, col = "red")
      plot(dates, weather_original$P - CASSIA_new_output[[3]][,c(variables_new[var])],
           main = "Residuals", xlab = "Date", ylab = "original - new output")
    }
  }

  ###
  # Sugar
  ###

  par(mfrow = c(2, 1), xpd=TRUE)
  plot(dates, CASSIA_new_output$Sugar$sugar_needles +
         CASSIA_new_output$Sugar$sugar_phloem +
         CASSIA_new_output$Sugar$sugar_roots +
         CASSIA_new_output$Sugar$sugar_xylem_sh +
         CASSIA_new_output$Sugar$sugar_xylem_st, main = "Sugar", col = "blue",
       type = "l", xlab = "Days of the Year", ylab = "Sugar, kg C", ylim = c(0, max(CASSIA_new_output$Sugar$sugar_needles +
                                                                                      CASSIA_new_output$Sugar$sugar_phloem +
                                                                                      CASSIA_new_output$Sugar$sugar_roots +
                                                                                      CASSIA_new_output$Sugar$sugar_xylem_sh +
                                                                                      CASSIA_new_output$Sugar$sugar_xylem_st)))
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
                                                                                CASSIA_new_output$Sugar$starch_xylem_st)))
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

}

