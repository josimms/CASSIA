symphony_plus_vs_original_plot <- function() {
  symphony_parameters = c(0.0317917, # A
                          0.016906, # s
                          0.0368857, # r
                          0.00505757, # m_p
                          3.98799e-4, # e_p
                          0.00369772, # r_p
                          0.0121216, # k
                          0.0909091, # alpha
                          0.0142857, # beta
                          0.0110068, # i
                          0.00262647, # l
                          4.22868e-4, # y
                          0.00929094, # u
                          0.0289652, # e
                          400, # Ca
                          0.0627704) # psi i

  Photosynthesis = PRELES(data_format$PAR,
                          data_format$T,
                          data_format$VPD,
                          data_format$Rain,
                          data_format$CO2,
                          data_format$fAPAR,
                          rep(1, length = nrow(data_format)))$GPP

  symphony_results = symphony(symphony_parameters)
  symphony_plus_results = symphony_plus(symphony_parameters, Photosynthesis)

  par(mfrow = c(3, 2))
  names = c("Carbon Decomposers FOM", "Carbon Decomposers SOM", "Carbon Plant", "Carbon FOM", "Mineral N in Soil", "Carbon in SOM")
  for (i in 1:length(symphony_plus_results)) {
    axis_min = min(symphony_plus_results[[i]], symphony_results[[i]])
    axis_max = max(symphony_plus_results[[i]], symphony_results[[i]])
    plot(symphony_plus_results[[i]], col = "blue",
         main = names[i],
         xlab = "Days", ylim = c(axis_min, axis_max),
         ylab = names(symphony_plus_results)[i], type = "l")
    lines(symphony_results[[i]], col = "black")
  }
}

symphony_multiple_FOM_daily_plot <- function(Tsb, SoilWater,
                                             C_FOM_needles, C_FOM_woody, C_FOM_roots, C_FOM_mantle, C_FOM_ERM,
                                             C_SOM, N_SOM,
                                             C_decompose_FOM, C_decompose_SOM, N_decompose_FOM, N_decompose_SOM,
                                             imobilisation,  assimilation,
                                             Litter_needles, Litter_woody, Litter_roots, Litter_mantle, Litter_ERM,
                                             NH4, NO3, N_FOM_needles, N_FOM_woody, N_FOM_roots, N_FOM_mantle, N_FOM_ERM,
                                             NH4_used_Plant, NH4_used_Fungal, NO3_used_Plant, NO3_used_Fungal, FOM_Norg_used_Plant, FOM_Norg_used_Fungal,
                                             SOM_Norg_used,
                                             Rm, Q10, N_limits_microbes, N_k_microbes, SWC_k_microbes,
                                             NC_microbe_opt, microbe_turnover) {

  ####
  # Simulations
  ####

  ### Initialisation

  symphony_temperature <- symphony_soilwater <-
    symphony_C_FOM_needles <- symphony_C_FOM_woody <- symphony_C_FOM_roots <- symphony_C_FOM_mycelium <-
    symphony_C_SOM <- symphony_N_SOM <- symphony_C_decompose_FOM <- symphony_C_decompose_SOM <- symphony_N_decompose_FOM <- symphony_N_decompose_SOM <-
    symphony_Litter_needles <- symphony_Litter_woody <- symphony_Litter_roots <- symphony_Litter_mycelium <-
    symphony_NH4 <- symphony_NO3 <-  symphony_N_FOM_needles <-  symphony_N_FOM_woody <-  symphony_N_FOM_roots <-  symphony_N_FOM_mycelium <-
    symphony_NH4_used_Plant <-  symphony_NH4_used_Fungal <-  symphony_NO3_used_Plant <-  symphony_NO3_used_Fungal <-
    symphony_FOM_Norg_used_Plant <-  symphony_FOM_Norg_used_Fungal <- symphony_SOM_Norg_used <- symphony_NC_microbe_opt <- symphony_microbe_turnover <-
    NULL

  symphony_multiple_FOM_daily(Tsb, SoilWater,
                              C_FOM_needles, C_FOM_woody, C_FOM_roots, C_FOM_mantle, C_FOM_ERM,
                              C_SOM, N_SOM,
                              C_decompose_FOM, C_decompose_SOM, N_decompose_FOM, N_decompose_SOM,
                              Litter_needles, Litter_woody, Litter_roots, Litter_mantle, Litter_ERM,
                              imobilisation,  assimilation,
                              NH4, NO3, N_FOM_needles, N_FOM_woody, N_FOM_roots, N_FOM_mantle, N_FOM_ERM,
                              NH4_used_Plant, NH4_used_Fungal, NO3_used_Plant, NO3_used_Fungal, FOM_Norg_used_Plant, FOM_Norg_used_Fungal,
                              SOM_Norg_used,
                              c(Rm, Q10), N_limits_microbes, N_k_microbes, SWC_k_microbes,
                              NC_microbe_opt, microbe_turnover)

  C_range = seq(0, 25, by = 0.01)
  N_range = seq(0, 25, by = 0.01)
  rate_range = seq(0, 1, by = 0.01)

  # TODO: would like to have all of these outputs as a graph, so how should I go about this? Is there a way of doing it where the output could be a list rather than having to do 8 for loops for everything?
  #for (C in C_range) {
  #
  #}

  #for (N in N_range) {
  #
  #}

  #for (rate in rate_range) {
  #
  #}


  ####
  # Plots
  ####

  par(2, 2)
  plot(C_range, symphony_temperature,
       xlab = "Temperature range",
       ylab = "?")

  plot(C_range, symphony_soilwater,
       xlab = "Soil water range",
       ylab = "?")

  plot(C_range, symphony_C_FOM_needles,
       xlab = "C range",
       ylab = "?",
       ylim = range(c(symphony_C_FOM_needles, symphony_C_FOM_woody, symphony_C_FOM_roots, symphony_C_FOM_mycelium)))
  points(C_range, symphony_C_FOM_woody)
  points(C_range, symphony_C_FOM_roots)
  points(C_range, symphony_C_FOM_mycelium)
  points(C_range, symphony_C_SOM)
  legend(C_range[1], max(symphony_C_FOM_woody), c("Needles", "Woody", "Roots", "Mycelim"), title = "FOM Compartment")

  plot(N_range, symphony_N_FOM_needles,
       xlab = "N range",
       ylab = "?",
       ylim = range(c(symphony_N_FOM_needles, symphony_N_FOM_woody, symphony_N_FOM_roots, symphony_N_FOM_mycelium,
                      symphony_NH4, symphony_NO3, symphony_N_SOM)))
  points(N_range, symphony_N_FOM_woody)
  points(N_range, symphony_N_FOM_roots)
  points(N_range, symphony_N_FOM_mycelium)
  points(N_range, symphony_NH4)
  points(N_range, symphony_NO3)
  points(N_range, symphony_N_SOM)
  legend(N_range[1], max(symphony_N_FOM_woody), c("Needles", "Woody", "Roots", "Mycelim"), title = "FOM Compartment")

  plot(rate_range, symphony_Litter_needles,
       xlab = "Litter decomposition range",
       ylab = "?",
       ylim = range(c(symphony_Litter_needles, symphony_Litter_woody, symphony_Litter_roots, symphony_Litter_mycelium)))
  points(rate_range, symphony_Litter_woody)
  points(rate_range, symphony_Litter_roots)
  points(rate_range, symphony_Litter_mycelium)
  legend(rate_range[1], max(symphony_Litter_woody), c("Needles", "Woody", "Roots", "Mycelim"), title = "FOM Compartment")

  plot(rate_range, symphony_NC_microbe_opt,
       xlab = "Rate range",
       ylab = "？",
       ylim = range(c(symphony_NC_microbe_opt, symphony_microbe_turnover)))
  points(rate_range, symphony_microbe_turnover)

}

nitrogen_graphs_plant <- function(Temperature,
                                  SoilWater,
                                  micorization,
                                  NH4,
                                  NO3,
                                  Norg,
                                  N_limits_R,
                                  N_k_R,
                                  SWC_k_R,
                                  parameters,
                                  demand) {
  NH4_high = 20
  NH4_medium = 0.5*NH4_high
  NH4_low = 0.01*NH4_high

  # Indervidual functions should follow this pattern
  temp_change_N_Uptake <- SWC_change_N_Uptake <- conc_change_N_Uptake <- NULL
  temp_change_N_Uptake_Fungal <- SWC_change_N_Uptake_Fungal <- conc_change_N_Uptake_Fungal <- conc_change_Norg_Uptake_Fungal <- conc_change_NO3_Uptake_Fungal <- conc_change_NH4_Uptake_Fungal <- NULL
  temp_change_N_Uptake_Plant <- SWC_change_N_Uptake_Plant <- conc_change_NH4_Uptake_Plant_NH4_out <- conc_change_NH4_Uptake_Plant_NO3_out <- conc_change_NH4_Uptake_Plant_Norg_out <- NULL
  conc_change_Norg_Uptake_Plant <- conc_change_NO3_Uptake_Plant_NH4_20 <- conc_change_NO3_Uptake_Plant_NH4_10 <- conc_change_NO3_Uptake_Plant_NH4_5 <- conc_change_NO3_Uptake_Plant_NH4_2 <- conc_change_NH4_Uptake_Plant <- NULL
  conc_change_N_Uptake_Plant_NH4_out <- conc_change_N_Uptake_Plant_NO3_out <- conc_change_N_Uptake_Plant_Norg_out <- NULL
  conc_change_NH4_Uptake_Plant_NH4_out <- conc_change_NH4_Uptake_Plant_NO3_out <- conc_change_NH4_Uptake_Plant_Norg_out <- NULL
  conc_change_NO3_Uptake_Plant_NH4_out <- conc_change_NO3_Uptake_Plant_NO3_out <- conc_change_NO3_Uptake_Plant_Norg_out <- NULL
  conc_change_NO3_Uptake_Plant_NH4_10_NH4_out <- conc_change_NO3_Uptake_Plant_NH4_10_NO3_out <- conc_change_NO3_Uptake_Plant_NH4_10_Norg_out <- NULL
  conc_change_Norg_Uptake_Plant_NH4_out <- conc_change_Norg_Uptake_Plant_NO3_out <- conc_change_Norg_Uptake_Plant_Norg_out <- conc_change_N_Uptake_Plant <- NULL
  temperature_range = -20:20
  count = 1
  for (Temp in temperature_range) {
    temp_change_N_Uptake[count] = uptake_N(10, Temp, 10, 0.2, 0.8, 0.5)
    temp_change_N_Uptake_Plant[count] = Plant_N_Uptake(Temp, SoilWater, micorization, NH4, NO3, Norg, N_limits_R, N_k_R, SWC_k_R, parameters, demand)$N_to_plant
    count = count + 1
  }
  concentration_range = seq(0, 600, by = 0.01)
  count = 1
  for (conc in concentration_range) {
    conc_change_N_Uptake[count] <- uptake_N(conc, Temperature, 10, 0.2, SoilWater, 0.5)
    conc_change_N_Uptake_Plant[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, conc, conc, N_limits_R, N_k_R, SWC_k_R, parameters, demand)$N_to_plant
    conc_change_N_Uptake_Plant_NH4_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, conc, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NH4_used
    conc_change_N_Uptake_Plant_NO3_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, conc, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NO3_used
    conc_change_N_Uptake_Plant_Norg_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, conc, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$Norg_used
    conc_change_Norg_Uptake_Plant[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, NO3, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$N_to_plant
    conc_change_Norg_Uptake_Plant_NH4_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, NO3, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NH4_used
    conc_change_Norg_Uptake_Plant_NO3_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, NO3, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NO3_used
    conc_change_Norg_Uptake_Plant_Norg_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, NO3, conc, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$Norg_used
    conc_change_NO3_Uptake_Plant_NH4_10[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$Norg_used
    conc_change_NO3_Uptake_Plant_NH4_10_NH4_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NH4_used
    conc_change_NO3_Uptake_Plant_NH4_10_NO3_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NO3_used
    conc_change_NO3_Uptake_Plant_NH4_10_Norg_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$Norg_used
    conc_change_NO3_Uptake_Plant_NH4_20[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4_high, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$N_to_plant
    conc_change_NO3_Uptake_Plant_NH4_5[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4_medium, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$N_to_plant
    conc_change_NO3_Uptake_Plant_NH4_2[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, NH4_low, conc, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$N_to_plant
    conc_change_NH4_Uptake_Plant[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, NO3, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$N_to_plant
    conc_change_NH4_Uptake_Plant_NH4_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, NO3, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NH4_used
    conc_change_NH4_Uptake_Plant_NO3_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, NO3, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$NO3_used
    conc_change_NH4_Uptake_Plant_Norg_out[count] = Plant_N_Uptake(Temperature, SoilWater, micorization, conc, NO3, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$Norg_used
    count = count + 1
  }
  SWC_range = seq(0, 1, by = 0.01)
  count = 1
  for (SWC in SWC_range) {
    SWC_change_N_Uptake[count] <- uptake_N(10, 15, 10, 0.2, SWC, 0.5)
    SWC_change_N_Uptake_Plant[count] = Plant_N_Uptake(Temperature, SWC, micorization, NH4, NO3, Norg, N_limits_R, N_k_R, SWC_k_R,parameters, demand)$N_to_plant
    count = count + 1
  }

  par(mfrow = c(3, 2))
  plot(concentration_range, conc_change_NH4_Uptake_Plant, xlab = "NH4", ylab = "total N to plant",
       main = "NH4 response (Plant) \n with temperature 15, SWC 0.8, other N 10",
       ylim = c(0, max(conc_change_NH4_Uptake_Plant)))
  lines(concentration_range, conc_change_NH4_Uptake_Plant_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NH4_Uptake_Plant_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NH4_Uptake_Plant_Norg_out, col = "green", lty = 2, lwd = 2)
  legend(0, max(unlist(c(conc_change_NH4_Uptake_Plant, conc_change_NH4_Uptake_Plant_NH4_out, conc_change_NH4_Uptake_Plant_NO3_out, conc_change_NH4_Uptake_Plant_Norg_out))),
         c("NH4", "No3", "FOM"), col = c("black", "blue", "green"), lty = 2, title = "N Type", bty = "n")

  plot(concentration_range, conc_change_Norg_Uptake_Plant, xlab = "Norg", ylab = "total N to plant",
       main = "Norg response (Plant) \n with temperature 15, SWC 0.8, other N 10",
       ylim = c(0, max(conc_change_Norg_Uptake_Plant)))
  lines(concentration_range, conc_change_Norg_Uptake_Plant_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_Norg_Uptake_Plant_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_Norg_Uptake_Plant_Norg_out, col = "green", lty = 2, lwd = 2)

  plot(concentration_range, conc_change_NO3_Uptake_Plant_NH4_10, xlab = "NO3", ylab = "total N to plant",
       main = "NO3 response (Plant) \n with temperature 15, SWC 0.8, other 10",
       ylim = c(0,  max(c(conc_change_NO3_Uptake_Plant_NH4_10))))
  lines(concentration_range, conc_change_NO3_Uptake_Plant_NH4_10_NH4_out, lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NO3_Uptake_Plant_NH4_10_NO3_out, lty = 2, lwd = 2, col = "blue")
  lines(concentration_range, conc_change_NO3_Uptake_Plant_NH4_10_Norg_out, lty = 2, lwd = 2, col = "green")

  plot(concentration_range, conc_change_NO3_Uptake_Plant_NH4_20, xlab = "NO3", ylab = "total N to plant",
       main = "NO3 response (Plant) \n with temperature 15, SWC 0.8, Organic 10, NH4 variable",
       ylim = c(0,  max(c(conc_change_NO3_Uptake_Plant_NH4_20, conc_change_NO3_Uptake_Plant_NH4_5, conc_change_NO3_Uptake_Plant_NH4_2))))
  points(concentration_range, conc_change_NO3_Uptake_Plant_NH4_5, col = "pink")
  points(concentration_range, conc_change_NO3_Uptake_Plant_NH4_2, col = "red")
  legend(0, max(unlist(c(conc_change_NO3_Uptake_Plant_NH4_10, conc_change_NO3_Uptake_Plant_NH4_5, conc_change_NO3_Uptake_Plant_NH4_2))),
         c("High", "Medium", "Low"), col = c("black", "pink", "red"), pch = 1, title = "NH4", bty = "n")

  plot(temperature_range, temp_change_N_Uptake_Plant, xlab = "Temperature", ylab = "total N to plant", main = "Temperature response (Plant) \n with concentation 10, SWC 0.8")

  plot(SWC_range, SWC_change_N_Uptake, xlab = "SWC", ylab = "total N to plant", main = "SWC response (Plant) \n with concentation 10, temperature 15")

}

nitrogen_graphs_fungal <- function(Temperature,
                                   SoilWater,
                                   NH4, NO3, Norg,
                                   N_limits_fungal,
                                   N_k_fungal,
                                   SWC_k_fungal,
                                   demand) {

  NH4_high = 20
  NH4_medium = 0.5*NH4_high
  NH4_low = 0.01*NH4_high

  #####
  ## Generating the values
  #####
  temp_change_N_Uptake_Fungal <- SWC_change_N_Uptake_Fungal <- conc_change_N_Uptake_Fungal <- conc_change_Norg_Uptake_Fungal <- conc_change_NO3_Uptake_Fungal <- conc_change_NH4_Uptake_Fungal <- NULL
  conc_change_N_Uptake_Fungal_NH4_out <- conc_change_N_Uptake_Fungal_NO3_out <- conc_change_N_Uptake_Fungal_Norg_out <- NULL
  conc_change_NO3_Uptake_Fungal_NH4_10 <- conc_change_NO3_Uptake_Fungal_NH4_10_NH4_out <- conc_change_NO3_Uptake_Fungal_NH4_10_NO3_out <- conc_change_NO3_Uptake_Fungal_NH4_10_Norg_out <- NULL
  conc_change_NO3_Uptake_Fungal_NH4_5 <- conc_change_NO3_Uptake_Fungal_NH4_2 <- SWC_change_N_Uptake <- NULL
  conc_change_NH4_Uptake_Fungal_NH4_out <- conc_change_NH4_Uptake_Fungal_NO3_out <- conc_change_NH4_Uptake_Fungal_Norg_out <- NULL
  conc_change_Norg_Uptake_Fungal_NH4_out <- conc_change_Norg_Uptake_Fungal_NO3_out <- conc_change_Norg_Uptake_Fungal_Norg_out <- NULL

  temperature_range = -20:20
  count = 1
  for (Temp in temperature_range) {
    temp_change_N_Uptake_Fungal[count] = Fungal_N_Uptake(Temp, SoilWater, NH4, NO3, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    count = count + 1
  }
  concentration_range = seq(0, 20, by = 0.01)
  count = 1
  for (conc in concentration_range) {
    conc_change_N_Uptake_Fungal[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                         conc, conc, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    conc_change_N_Uptake_Fungal_NH4_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                 conc, conc, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NH4_used
    conc_change_N_Uptake_Fungal_NO3_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                 conc, conc, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NO3_used
    conc_change_N_Uptake_Fungal_Norg_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                  conc, conc, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$Norg_used
    conc_change_NO3_Uptake_Fungal_NH4_10[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                  NH4, conc, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    conc_change_NO3_Uptake_Fungal_NH4_10_NH4_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                          NH4, conc, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NH4_used
    conc_change_NO3_Uptake_Fungal_NH4_10_NO3_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                          NH4, conc, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NO3_used
    conc_change_NO3_Uptake_Fungal_NH4_10_Norg_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                           NH4, conc, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$Norg_used
    conc_change_NO3_Uptake_Fungal_NH4_5[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                 NH4_medium, conc, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    conc_change_NO3_Uptake_Fungal_NH4_2[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                 NH4_low, conc, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    conc_change_NH4_Uptake_Fungal[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                           conc, NO3, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    conc_change_NH4_Uptake_Fungal_NH4_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                   conc, NO3, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NH4_used
    conc_change_NH4_Uptake_Fungal_NO3_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                   conc, NO3, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NO3_used
    conc_change_NH4_Uptake_Fungal_Norg_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                    conc, NO3, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$Norg_used
    conc_change_Norg_Uptake_Fungal[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                            NH4, NO3, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    conc_change_Norg_Uptake_Fungal_NH4_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                    NH4, NO3, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NH4_used
    conc_change_Norg_Uptake_Fungal_NO3_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                    NH4, NO3, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$NO3_used
    conc_change_Norg_Uptake_Fungal_Norg_out[count] = Fungal_N_Uptake(Temperature, SoilWater,
                                                                     NH4, NO3, conc, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$Norg_used
    count = count + 1
  }

  SWC_range = seq(0, 1, by = 0.01)
  count = 1
  for (SWC in SWC_range) {
    SWC_change_N_Uptake_Fungal[count] = Fungal_N_Uptake(Temperature, SWC,
                                                        NH4, NO3, Norg, N_limits_fungal, N_k_fungal, SWC_k_fungal, demand)$N_to_fungal
    count = count + 1
  }

  #####
  ## Plots
  #####
  par(mfrow = c(3, 2))
  plot(concentration_range, conc_change_NH4_Uptake_Fungal, xlab = "NH4", ylab = "total N to Fungal", main = "NH4 response (Fungal) \n with temperature 15, SWC 0.8, other N 10",
       ylim = c(0, max(conc_change_NH4_Uptake_Fungal)))
  lines(concentration_range, conc_change_NH4_Uptake_Fungal_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NH4_Uptake_Fungal_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NH4_Uptake_Fungal_Norg_out, col = "green", lty = 2, lwd = 2)
  legend(0, max(unlist(c(conc_change_NH4_Uptake_Fungal, conc_change_NH4_Uptake_Fungal_NH4_out, conc_change_NH4_Uptake_Fungal_NO3_out, conc_change_NH4_Uptake_Fungal_Norg_out))),
         c("NH4", "NO3", "FOM"), col = c("black", "blue", "green"), lty = 2, title = "N Type", bty = "n")

  plot(concentration_range, conc_change_Norg_Uptake_Fungal, xlab = "Norg", ylab = "total N to Fungal",
       main = "Norg response (Fungal) \n with temperature 15, SWC 0.8, other N 10",
       ylim = c(0, max(conc_change_Norg_Uptake_Fungal)))
  lines(concentration_range, conc_change_Norg_Uptake_Fungal_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_Norg_Uptake_Fungal_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_Norg_Uptake_Fungal_Norg_out, col = "green", lty = 2, lwd = 2)

  plot(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_10, xlab = "NO3", ylab = "total N to Fungal",
       main = "NO3 response (Fungal) \n with temperature 15, SWC 0.8, other 10",
       ylim = c(0,  max(c(conc_change_NO3_Uptake_Fungal_NH4_10))))
  lines(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_10_NH4_out, lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_10_NO3_out, lty = 2, lwd = 2, col = "blue")
  lines(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_10_Norg_out, lty = 2, lwd = 2, col = "green")

  plot(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_10, xlab = "NO3", ylab = "total N to Fungal",
       main = "NO3 response (Fungal) \n with temperature 15, SWC 0.8, Organic 10, NH4 variable",
       ylim = c(0,  max(c(conc_change_NO3_Uptake_Fungal_NH4_10, conc_change_NO3_Uptake_Fungal_NH4_5, conc_change_NO3_Uptake_Fungal_NH4_2))))
  points(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_5, col = "pink")
  points(concentration_range, conc_change_NO3_Uptake_Fungal_NH4_2, col = "red")
  legend(0, max(unlist(c(conc_change_NO3_Uptake_Fungal_NH4_10, conc_change_NO3_Uptake_Fungal_NH4_5, conc_change_NO3_Uptake_Fungal_NH4_2))),
         c("High", "Medium", "Low"), col = c("black", "pink", "red"), pch = 1, title = "NH4", bty = "n")

  plot(temperature_range, temp_change_N_Uptake_Fungal, xlab = "Temperature", ylab = "total N to Fungal", main = "Temperature response (Fungal) \n with concentation 10, SWC 0.8")

  plot(SWC_range, SWC_change_N_Uptake, xlab = "SWC", ylab = "total N to Fungal", main = "SWC response (Fungal) \n with concentation 10, temperature 15")
}

nitrogen_graphs <- function(N_type, Temperature, N_limit, k, SoilWater, SoilWaterk) {
  temp_change_N_Uptake <- SWC_change_N_Uptake <- conc_change_N_Uptake <- NULL

  temperature_range = -20:20
  count = 1
  for (Temp in temperature_range) {
    temp_change_N_Uptake[count] = uptake_N(N_type, Temp, N_limit, k, SoilWater, SoilWaterk)
    count = count + 1
  }

  concentration_range = seq(0, 600, by = 0.01)
  count = 1
  for (conc in concentration_range) {
    conc_change_N_Uptake[count] <- uptake_N(conc, Temperature, N_limit, k, SoilWater, SoilWaterk)
    count = count + 1
  }

  SWC_range = seq(0, 1, by = 0.01)
  count = 1
  for (SWC in SWC_range) {
    SWC_change_N_Uptake[count] <- uptake_N(N_type, Temperature, N_limit, k, SWC, SoilWaterk)
    count = count + 1
  }

  data_path = "~/Documents/CASSIA_Calibration/Processed_Data/"
  data = read.delim(paste0(data_path, "oyewole_2015_calibration_data.csv"), sep = ",", dec = ",")

  par(mfrow = c(3, 1))
  plot(concentration_range, conc_change_N_Uptake, xlab = "Concentration", ylab = "Uptake of N",
       main = "Concentration response \n with temperature 15, SWC 0.8", type = "l",
       ylim = c(0, max(data$control, conc_change_N_Uptake)))
  points(data$concentration, data$control, col = as.factor(data$n_comp))
  arrows(data$concentration, data$control - data$error,
         data$concentration, data$control + data$error, length=0.05, angle=90, code=3,
         col = as.factor(data$n_comp))
  legend(0, max(conc_change_N_Uptake, data$control), c("Agrinine", "Glycine", "NH4", "NO3"), col = c(1, 2, 3, 4), bty = "n", pch = 1)
  plot(temperature_range, temp_change_N_Uptake, xlab = "Temperature", ylab = "Uptake of N", main = "Temperature response \n with concentation 10, SWC 0.8", type = "l")
  plot(SWC_range, SWC_change_N_Uptake, xlab = "SWC", ylab = "Uptake of N", main = "SWC response \n with concentation 10, temperature 15", type = "l")
}

nitrogen_graphs_microbe <- function(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, Rm, Q10) {

  #####
  ## Generating the values
  #####
  temp_change_N_Uptake_Microbes_C_out <- SWC_change_N_Uptake_Microbes_C_out <- conc_change_N_Uptake_Microbes_C_out <- conc_change_Norg_Uptake_Microbes_C_out <- conc_change_NO3_Uptake_Microbes_C_out <- conc_change_NH4_Uptake_Microbes_C_out <- NULL
  conc_change_N_Uptake_Microbes_NH4_out <- conc_change_N_Uptake_Microbes_NO3_out <- conc_change_N_Uptake_Microbes_Norg_out <- NULL
  temp_change_N_Uptake_Microbes_NH4_out <- temp_change_N_Uptake_Microbes_NO3_out <- temp_change_N_Uptake_Microbes_Norg_out <- NULL
  conc_change_NO3_Uptake_Microbes_NH4_1000_C_out <- conc_change_NO3_Uptake_Microbes_NH4_1000_NH4_out <- conc_change_NO3_Uptake_Microbes_NH4_1000_NO3_out <- conc_change_NO3_Uptake_Microbes_NH4_1000_Norg_out <- NULL
  conc_change_NO3_Uptake_Microbes_NH4_500_C_out <- conc_change_NO3_Uptake_Microbes_NH4_200_C_out <- NULL
  conc_change_NH4_Uptake_Microbes_NH4_out <- conc_change_NH4_Uptake_Microbes_NO3_out <- conc_change_NH4_Uptake_Microbes_Norg_out <- NULL
  conc_change_Norg_Uptake_Microbes_NH4_out <- conc_change_Norg_Uptake_Microbes_NO3_out <- conc_change_Norg_Uptake_Microbes_Norg_out <- NULL
  conc_change_Norg_FOM_Uptake_Microbes <- NULL

  temperature_range = -20:20
  count = 1
  for (Temp in temperature_range) {
    temp_change_N_Uptake_Microbes_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, Norg, Temp, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    temp_change_N_Uptake_Microbes_NH4_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, Norg, Temp, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NH4_uptaken
    temp_change_N_Uptake_Microbes_NO3_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, Norg, Temp, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NO3_uptaken
    temp_change_N_Uptake_Microbes_Norg_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, Norg, Temp, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$Norg_uptaken
    count = count + 1
  }

  concentration_range = seq(0, 30, by = 1)
  count = 1
  for (conc in concentration_range) {
    conc_change_N_Uptake_Microbes_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, conc, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    conc_change_N_Uptake_Microbes_NH4_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, conc, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NH4_uptaken
    conc_change_N_Uptake_Microbes_NO3_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, conc, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NO3_uptaken
    conc_change_N_Uptake_Microbes_Norg_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, conc, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$Norg_uptaken
    conc_change_NO3_Uptake_Microbes_NH4_1000_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, conc, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    conc_change_NO3_Uptake_Microbes_NH4_1000_NH4_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, conc, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NH4_uptaken
    conc_change_NO3_Uptake_Microbes_NH4_1000_NO3_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, conc, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NO3_uptaken
    conc_change_NO3_Uptake_Microbes_NH4_1000_Norg_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, conc, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$Norg_uptaken
    conc_change_NO3_Uptake_Microbes_NH4_500_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, conc, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    conc_change_NO3_Uptake_Microbes_NH4_200_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, conc, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    conc_change_NH4_Uptake_Microbes_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, NO3, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    conc_change_NH4_Uptake_Microbes_NH4_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, NO3, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NH4_uptaken
    conc_change_NH4_Uptake_Microbes_NO3_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, NO3, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NO3_uptaken
    conc_change_NH4_Uptake_Microbes_Norg_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, conc, NO3, Norg, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$Norg_uptaken
    conc_change_Norg_Uptake_Microbes_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    conc_change_Norg_Uptake_Microbes_NH4_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NH4_uptaken
    conc_change_Norg_Uptake_Microbes_NO3_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$NO3_uptaken
    conc_change_Norg_Uptake_Microbes_Norg_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, conc, Temperature, SoilWater, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$Norg_uptaken
    count = count + 1
  }

  SWC_range = seq(0, 1, by = 0.01)
  count = 1
  for (SWC in SWC_range) {
    SWC_change_N_Uptake_Microbes_C_out[count] = Microbe_Uptake(C_microbe, N_microbe, C_soil_component, NC_microbe_opt, NH4, NO3, Norg, Temperature, SWC, NC_Litter, imnolisation, assimilation, N_limits_microbes, N_k_microbes, SWC_k_microbes, SOM_decomposers, c(Rm, Q10))$C_uptaken
    count = count + 1
  }

  #####
  ## Plots
  #####

  par(mfrow = c(3, 2))
  plot(concentration_range, conc_change_NH4_Uptake_Microbes_C_out, xlab = "NH4", ylab = "total N to Microbes",
       main = "NH4 response (Microbes) \n with temperature 15, SWC 0.8, other N 1000", col = "red",
       ylim = range(c(conc_change_NH4_Uptake_Microbes_C_out, conc_change_NH4_Uptake_Microbes_NH4_out,
                      conc_change_NH4_Uptake_Microbes_NO3_out, conc_change_NH4_Uptake_Microbes_Norg_out)))
  lines(concentration_range, conc_change_NH4_Uptake_Microbes_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NH4_Uptake_Microbes_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NH4_Uptake_Microbes_Norg_out, col = "green", lty = 2, lwd = 2)
  legend(0, max(unlist(c(conc_change_NH4_Uptake_Microbes_C_out, conc_change_NH4_Uptake_Microbes_NH4_out, conc_change_NH4_Uptake_Microbes_NO3_out, conc_change_NH4_Uptake_Microbes_Norg_out))),
         c("NH4", "NO3", "FOM", "C"), col = c("black", "blue", "green", "red"), lty = 2, title = "N Type", bty = "n")

  plot(concentration_range, conc_change_NO3_Uptake_Microbes_NH4_1000_C_out, xlab = "NO3", ylab = "total N to Microbes",
       main = "NO3 response (Microbes) \n with temperature 15, SWC 0.8, other 1000", col = "red",
       ylim = range(c(conc_change_NO3_Uptake_Microbes_NH4_1000_C_out, conc_change_NO3_Uptake_Microbes_NH4_1000_NH4_out,
                      conc_change_NO3_Uptake_Microbes_NH4_1000_NO3_out, conc_change_NO3_Uptake_Microbes_NH4_1000_Norg_out)))
  lines(concentration_range, conc_change_NO3_Uptake_Microbes_NH4_1000_NH4_out, lty = 2, lwd = 2)
  lines(concentration_range, conc_change_NO3_Uptake_Microbes_NH4_1000_NO3_out, lty = 2, lwd = 2, col = "blue")
  lines(concentration_range, conc_change_NO3_Uptake_Microbes_NH4_1000_Norg_out, lty = 2, lwd = 2, col = "green")

  plot(concentration_range, conc_change_Norg_Uptake_Microbes_C_out, xlab = "Norg", ylab = "total N to Microbes",
       main = "Norg response (Microbes) \n with temperature 15, SWC 0.8, other N 1000", col = "red",
       ylim = range(c(conc_change_Norg_Uptake_Microbes_C_out, conc_change_Norg_Uptake_Microbes_NH4_out,
                      conc_change_Norg_Uptake_Microbes_NO3_out, conc_change_Norg_Uptake_Microbes_Norg_out)))
  lines(concentration_range, conc_change_Norg_Uptake_Microbes_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_Norg_Uptake_Microbes_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_Norg_Uptake_Microbes_Norg_out, col = "green", lty = 2, lwd = 2)

  plot(concentration_range, conc_change_N_Uptake_Microbes_C_out, xlab = "N", ylab = "total N to Microbes",
       main = "Norg response (Microbes) \n with temperature 15, SWC 0.8, other N 1000", col = "red",
       ylim = range(c(conc_change_N_Uptake_Microbes_C_out, conc_change_N_Uptake_Microbes_NH4_out,
                      conc_change_N_Uptake_Microbes_NO3_out, conc_change_N_Uptake_Microbes_Norg_out), na.rm = T))
  lines(concentration_range, conc_change_N_Uptake_Microbes_NH4_out, col = "black", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_N_Uptake_Microbes_NO3_out, col = "blue", lty = 2, lwd = 2)
  lines(concentration_range, conc_change_N_Uptake_Microbes_Norg_out, col = "green", lty = 2, lwd = 2)

  plot(temperature_range, temp_change_N_Uptake_Microbes_C_out, xlab = "Temperature", ylab = "total N to Microbes",
       main = "Temperature response (Microbes) \n with concentation 10, SWC 0.8", col = "black")

  plot(SWC_range, SWC_change_N_Uptake_Microbes_C_out, xlab = "SWC", ylab = "total N to Microbes",
       main = "SWC response (Microbes) \n with concentation 10, temperature 15", col = "black")
}



decision_graphs <- function(C_roots, N_roots, C_fungal, N_fungal,
                            optimal_root_fungal_biomass_ratio,
                            N_allo, max_C_allocation_CASSIA,
                            NC_fungal_opt, growth_C, growth_N) {

  ####
  # Initialisation
  ####
  plant_decision_C_roots_mycofon <- plant_decision_C_roots_franklin <- plant_decision_C_fungal_mycofon <- plant_decision_C_fungal_franklin <-
    plant_transfer_C_roots_mycofon <- plant_transfer_C_roots_franklin <- plant_transfer_C_fungal_mycofon <- plant_transfer_C_fungal_franklin <-
    myco_decision_N_roots_mycofon <- myco_decision_N_roots_franklin <- plant_decision_N_roots_mycofon <- plant_decision_N_roots_franklin <-
    myco_transfer_N_roots_mycofon <- myco_transfer_N_roots_mycofon <- plant_transfer_N_roots_mycofon <- plant_transfer_N_roots_franklin <-
    myco_decision_N_fungal_mycofon <- myco_decision_N_fungal_franklin <- myco_transfer_N_fungal_mycofon <- myco_transfer_N_fungal_franklin <-
    plant_decision_op_mycofon <- plant_decision_op_franklin <- plant_transfer_op_mycofon <- plant_transfer_op_franklin <-
    myco_decision_op_fungal_mycofon <- myco_decision_op_fungal_franklin <- myco_transfer_op_fungal_mycofon <- myco_transfer_op_fungal_franklin <-
    myco_decision_growth_C_fungal_mycofon <- myco_decision_growth_C_fungal_franklin <- myco_transfer_growth_C_fungal_mycofon <- myco_transfer_growth_C_fungal_franklin <-
    myco_decision_growth_N_fungal_mycofon <- myco_decision_growth_N_fungal_franklin <- myco_transfer_growth_N_fungal_mycofon <- myco_transfer_growth_N_fungal_franklin <-
    plant_decision_N_allo_mycofon <- plant_decision_N_allo_franklin <- plant_transfer_N_allo_mycofon <- plant_transfer_N_allo_franklin <-
    plant_decision_CASSIA_mycofon <- plant_decision_CASSIA_franklin <- plant_transfer_CASSIA_mycofon <- plant_transfer_CASSIA_franklin <-
    myco_transfer_N_roots_franklin <-
    NULL

  ####
  # State tests
  ####
  C_range <- seq(100, 300)
  count = 1
  for (C in C_range) {
    plant_decision_C_roots_mycofon[count] <- plant_decision(C, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Mycofon_demand
    plant_decision_C_roots_franklin[count] <- plant_decision(C, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Franklin_demand

    plant_decision_C_fungal_mycofon[count] <- plant_decision(C_roots, N_roots, C, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Mycofon_demand
    plant_decision_C_fungal_franklin[count] <- plant_decision(C_roots, N_roots, C, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Franklin_demand

    plant_transfer_C_roots_mycofon[count] <- plant_decision(C, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Mycofon_allocation
    plant_transfer_C_roots_franklin[count] <- plant_decision(C, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Franklin_allocation

    plant_transfer_C_fungal_mycofon[count] <- plant_decision(C_roots, N_roots, C, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Mycofon_allocation
    plant_transfer_C_fungal_franklin[count] <- plant_decision(C_roots, N_roots, C, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Franklin_allocation

    count = count + 1
  }
  N_range <- seq(50, 250)
  count = 1
  for (N in N_range) {
    myco_decision_N_roots_mycofon[count] <- myco_decision(C_fungal, N, C_roots, N_roots, NC_fungal_opt, growth_C, growth_N)$Mycofon_demand
    myco_decision_N_roots_franklin[count] <- myco_decision(C_fungal, N, C_roots, N_roots, NC_fungal_opt, growth_C, growth_N)$Franklin_demand

    plant_decision_N_roots_mycofon[count] <- plant_decision(C_roots, N, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Mycofon_demand
    plant_decision_N_roots_franklin[count] <- plant_decision(C_roots, N, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Franklin_demand

    myco_transfer_N_roots_mycofon[count] <- myco_decision(C_fungal, N, C_roots, N_roots, NC_fungal_opt, growth_C, growth_N)$Mycofon_allocation
    myco_transfer_N_roots_franklin[count] <- myco_decision(C_fungal, N, C_roots, N_roots, NC_fungal_opt, growth_C, growth_N)$Franklin_allocation

    plant_transfer_N_roots_mycofon[count] <- plant_decision(C_roots, N, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Mycofon_allocation
    plant_transfer_N_roots_franklin[count] <- plant_decision(C_roots, N, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, max_C_allocation_CASSIA)$Franklin_allocation

    myco_decision_N_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N, NC_fungal_opt, growth_C, growth_N)$Mycofon_demand
    myco_decision_N_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N, NC_fungal_opt, growth_C, growth_N)$Franklin_demand

    myco_transfer_N_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N, NC_fungal_opt, growth_C, growth_N)$Mycofon_allocation
    myco_transfer_N_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N, NC_fungal_opt, growth_C, growth_N)$Franklin_allocation

    count = count + 1
  }
  optimum_point <- seq(0.2, 0.8, by = 0.01)
  count = 1
  for (op in optimum_point) {
    plant_decision_op_mycofon[count] <- plant_decision(C_roots, N_roots, C_fungal, op, N_allo, max_C_allocation_CASSIA)$Mycofon_demand
    plant_decision_op_franklin[count] <- plant_decision(C_roots, N_roots, C_fungal, op, N_allo, max_C_allocation_CASSIA)$Franklin_demand

    plant_transfer_op_mycofon[count] <- plant_decision(C, N_roots, C_fungal, op, N_allo, max_C_allocation_CASSIA)$Mycofon_allocation
    plant_transfer_op_franklin[count] <- plant_decision(C, N_roots, C_fungal, op, N_allo, max_C_allocation_CASSIA)$Franklin_allocation

    myco_decision_op_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, op, growth_C, growth_N)$Mycofon_demand
    myco_decision_op_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, op, growth_C, growth_N)$Franklin_demand

    myco_transfer_op_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, op, growth_C, growth_N)$Mycofon_allocation
    myco_transfer_op_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, op, growth_C, growth_N)$Franklin_allocation

    myco_decision_growth_C_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, op, growth_N)$Mycofon_demand
    myco_decision_growth_C_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, op, growth_N)$Franklin_demand

    myco_transfer_growth_C_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, op, growth_N)$Mycofon_allocation
    myco_transfer_growth_C_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, op, growth_N)$Franklin_allocation

    myco_decision_growth_N_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, growth_C, op)$Mycofon_demand
    myco_decision_growth_N_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, growth_C, op)$Franklin_demand

    myco_transfer_growth_N_fungal_mycofon[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, growth_C, op)$Mycofon_allocation
    myco_transfer_growth_N_fungal_franklin[count] <- myco_decision(C_fungal, N_fungal, C_roots, N_roots, NC_fungal_opt, growth_C, op)$Franklin_allocation

    plant_decision_N_allo_mycofon[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, op, max_C_allocation_CASSIA)$Mycofon_demand
    plant_decision_N_allo_franklin[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, op, max_C_allocation_CASSIA)$Franklin_demand

    plant_transfer_N_allo_mycofon[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, op, max_C_allocation_CASSIA)$Mycofon_allocation
    plant_transfer_N_allo_franklin[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, op, max_C_allocation_CASSIA)$Franklin_allocation

    plant_decision_CASSIA_mycofon[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, op)$Mycofon_demand
    plant_decision_CASSIA_franklin[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, op)$Franklin_demand

    plant_transfer_CASSIA_mycofon[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, op)$Mycofon_allocation
    plant_transfer_CASSIA_franklin[count] <- plant_decision(C_roots, N_roots, C_fungal, optimal_root_fungal_biomass_ratio, N_allo, op)$Franklin_allocation

    count = count + 1
  }


  ####
  # Plots
  ####

  par(mfrow = c(3, 2))
  plot(C_range, plant_transfer_C_roots_mycofon, col = "black",
       ylab = "Plant: Transfer of C (decision value assumed 1)", xlab = "C concentration",
       main = "C roots range test",
       ylim = range(c(plant_transfer_C_roots_mycofon, plant_transfer_C_roots_franklin)))
  points(C_range, plant_transfer_C_roots_franklin, col = "red")
  legend(C_range[1], max(plant_transfer_C_roots_franklin), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(C_range, plant_transfer_C_fungal_mycofon, col = "black",
       ylab = "Plant: Transfer of C (decision value assumed 1)", xlab = "C concentration",
       main = "C fungal range test",
       ylim = range(c(plant_transfer_C_fungal_mycofon, plant_transfer_C_fungal_franklin)))
  points(C_range, plant_transfer_C_fungal_franklin, col = "red")
  legend(C_range[1], max(plant_transfer_C_fungal_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(N_range, plant_transfer_N_roots_mycofon, col = "black",
       ylab = "Plant: Transfer of C (decision value assumed 1)", xlab = "C concentration",
       main = "N roots range test",
       ylim = range(c(plant_transfer_N_roots_mycofon, plant_transfer_N_roots_franklin)))
  points(N_range, plant_transfer_N_roots_franklin, col = "red")
  legend(N_range[1], max(plant_transfer_N_roots_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(N_range, myco_transfer_N_fungal_mycofon, col = "black",
       ylab = "Fungi: Transfer of N (decision value assumed 1)", xlab = "C concentration",
       main = "N fungal range test",
       ylim = range(c(myco_transfer_N_fungal_mycofon, myco_transfer_N_fungal_franklin)))
  points(N_range, myco_transfer_N_fungal_franklin, col = "red")
  legend(N_range[1], max(myco_transfer_N_fungal_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(N_range, myco_transfer_N_roots_mycofon, col = "black",
       ylab = "Fungi: Transfer of N (decision value assumed 1)", xlab = "C concentration",
       main = "N roots range test",
       ylim = range(c(myco_transfer_N_roots_mycofon, myco_transfer_N_roots_franklin)))
  points(N_range, myco_transfer_N_roots_franklin, col = "red")
  legend(N_range[1], max(myco_transfer_N_roots_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(optimum_point, plant_transfer_op_mycofon, col = "black",
       ylab = "Plant: Transfer of C (decision value assumed 1)", xlab = "C concentration",
       main = "Optimal value range test",
       ylim = range(c(plant_transfer_op_mycofon, plant_transfer_op_franklin)))
  points(optimum_point, plant_transfer_op_franklin, col = "red")
  legend(optimum_point[1], max(plant_transfer_op_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(optimum_point, myco_transfer_op_fungal_mycofon, col = "black",
       ylab = "Fungi: Transfer of C (decision value assumed 1)", xlab = "C concentration",
       main = "Optimal value range test",
       ylim = range(c(myco_transfer_op_fungal_mycofon, myco_transfer_op_fungal_franklin)))
  points(optimum_point, myco_transfer_op_fungal_franklin, col = "red")
  legend(optimum_point[1], max(myco_transfer_op_fungal_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(optimum_point, plant_transfer_N_allo_mycofon, col = "black",
       ylab = "Plant: Transfer of C (decision value assumed 1)", xlab = "C concentration",
       main = "N uptake range test",
       ylim = range(c(plant_transfer_N_allo_mycofon, plant_transfer_N_allo_franklin)))
  points(optimum_point, plant_transfer_N_allo_franklin, col = "red")
  legend(optimum_point[1], max(plant_transfer_N_allo_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(optimum_point, plant_transfer_CASSIA_mycofon, col = "black",
         ylab = "Plant: Transfer of C (decision value assumed 1)", xlab = "C concentration",
         main = "Maximum C for allocation from CASSIA range test",
         ylim = range(c(plant_transfer_CASSIA_mycofon, plant_transfer_CASSIA_franklin)))
  points(optimum_point, plant_transfer_CASSIA_franklin, col = "red")
  legend(optimum_point[1], max(plant_transfer_CASSIA_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(optimum_point, myco_transfer_growth_C_fungal_mycofon, col = "black",
       ylab = "Fungi: Transfer of N (decision value assumed 1)", xlab = "C concentration",
       main = "Growth C dependency parameter range test",
       ylim = range(c(myco_transfer_growth_C_fungal_mycofon, myco_transfer_growth_C_fungal_franklin)))
  points(optimum_point, myco_transfer_growth_C_fungal_franklin, col = "red")
  legend(optimum_point[1], max(myco_transfer_growth_C_fungal_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

  plot(optimum_point, myco_transfer_growth_N_fungal_mycofon, col = "black",
       ylab = "Fungi: Transfer of N (decision value assumed 1)", xlab = "C concentration",
       main = "Growth N dependency parameter range test",
       ylim = range(c(myco_transfer_growth_N_fungal_mycofon, myco_transfer_growth_N_fungal_franklin)))
  points(optimum_point, myco_transfer_growth_N_fungal_franklin, col = "red")
  legend(optimum_point[1], max(myco_transfer_growth_N_fungal_mycofon), c("Mycofon", "Franklin"), col = c("black", "red"), pch = c(1, 1),
         title = "Method", bty = "n")

}

respiration_graphs <- function(Rm, Q10) {
  temp_respiration <- NULL
  temp = -30:30
  count = 1
  for (tem in temp) {
    temp_respiration[count] <- respiration(tem, Rm, Q10)
    count = count + 1
  }

  par(mfrow = c(1, 1))
  plot(temp, temp_respiration, main = "Respiration",
       xlab = "Temperature range", ylab = "Respiration, kg C")
}

myco_growth_graphs <- function(C_fungal_in, N_fungal_in, growth_C, growth_N) {
  myco_growth_C_fungal <- myco_growth_N_fungal <- myco_growth_C_and_N_fungal <- NULL
  range <- -10:50
  count = 1
  for (ele in range) {
    myco_growth_C_and_N_fungal[count] <- myco_growth(ele, ele, growth_C, growth_N)$growth
    myco_growth_C_fungal[count] <- myco_growth(ele, N_fungal_in, growth_C, growth_N)$growth
    myco_growth_N_fungal[count] <- myco_growth(C_fungal_in, ele, growth_C, growth_N)$growth
    count = count + 1
  }

  myco_growth_yearly_N_5 <- myco_growth_yearly_N_10 <- myco_growth_yearly_N_20 <- NULL
  myco_growth_yearly_N_5[1] <- myco_growth_yearly_N_10[1] <- myco_growth_yearly_N_20[1] <- C_fungal_in
  for (t in 2:50) {
    myco_growth_yearly_N_5[t] <- myco_growth_yearly_N_5[t-1] + myco_growth(myco_growth_yearly_N_5[t-1], 5, growth_C, growth_N)$growth
    myco_growth_yearly_N_10[t] <- myco_growth_yearly_N_10[t-1] + myco_growth(myco_growth_yearly_N_10[t-1], 10, growth_C, growth_N)$growth
    myco_growth_yearly_N_20[t] <- myco_growth_yearly_N_20[t-1] + myco_growth(myco_growth_yearly_N_20[t-1], 20, growth_C, growth_N)$growth
  }

  par(mfrow = c(2, 2))
  plot(range, myco_growth_C_fungal, main = "Growth\n C concentration changed", xlab = "C concentration", ylab = "Growth")
  plot(range, myco_growth_N_fungal, main = "Growth\n N concentration changed", xlab = "N concentration", ylab = "Growth")
  plot(range, myco_growth_C_and_N_fungal, main = "Growth\n with N and C concentration changed", xlab = "C and N concentration", ylab = "Growth")
  plot(1:50, myco_growth_yearly_N_20, main = "Growth over time\n assuming N and C are 20 for each time interval", xlab = "Time", ylab = "Growth")
  points(1:50, myco_growth_yearly_N_10, col = "blue")
  points(1:50, myco_growth_yearly_N_5, col = "red")
  legend(0, max(myco_growth_yearly_N_20), c(5, 10, 20), title = "N = ", col = c("red", "blue", 1), pch = 1, bty = "n")
}

mycofon_balence_graphs <- function(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                   C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                   Rm, Q10, NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                   N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                   mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                   max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy) {
  ####
  # Initalisation
  ####

  mycofon_balence_C_roots_roots <- mycofon_balence_C_roots_N_roots <- mycofon_balence_C_roots_fungal <- mycofon_balence_C_roots_N_fungal <- NULL
  mycofon_balence_C_fungal_roots <- mycofon_balence_C_fungal_N_roots <- mycofon_balence_C_fungal_fungal <- mycofon_balence_C_fungal_N_fungal <- NULL
  mycofon_balence_N_roots_roots <- mycofon_balence_N_roots_N_roots <- mycofon_balence_N_roots_fungal <- mycofon_balence_N_roots_N_fungal <- NULL
  mycofon_balence_N_fungal_roots <- mycofon_balence_N_fungal_N_roots <- mycofon_balence_N_fungal_fungal <- mycofon_balence_N_fungal_N_fungal <- NULL
  mycofon_balence_optimal_root_fungal_biomass_ratio_roots <- mycofon_balence_optimal_root_fungal_biomass_ratio_N_roots <-
    mycofon_balence_optimal_root_fungal_biomass_ratio_fungal <- mycofon_balence_optimal_root_fungal_biomass_ratio_N_fungal <- NULL
  mycofon_balence_NH4_roots <- mycofon_balence_NH4_N_roots <- mycofon_balence_NH4_fungal <- mycofon_balence_NH4_N_fungal <- NULL
  mycofon_balence_NO3_roots <- mycofon_balence_NO3_N_roots <- mycofon_balence_NO3_fungal <- mycofon_balence_NO3_N_fungal <- NULL
  mycofon_balence_FOM_Norg_roots <- mycofon_balence_FOM_Norg_N_roots <- mycofon_balence_FOM_Norg_fungal <- mycofon_balence_FOM_Norg_N_fungal <- NULL
  mycofon_balence_T_roots <- mycofon_balence_Tsb_roots <- mycofon_balence_T_N_fungal <- mycofon_balence_Tsb_N_fungal <- NULL
  mycofon_balence_SWC_roots <- mycofon_balence_SWC_N_roots <- mycofon_balence_SWC_fungal <- mycofon_balence_SWC_N_fungal <- NULL
  mycofon_balence_T_fungal <- mycofon_balence_Tsb_fungal <- mycofon_balence_T_N_roots <- mycofon_balence_Tsb_N_roots <- NULL
  mycofon_balence_mm_roots <- mycofon_balence_mm_N_roots <- mycofon_balence_mm_fungal <- mycofon_balence_mm_N_fungal <- NULL
  mycofon_balence_ERM_roots <- mycofon_balence_ERM_N_roots <- mycofon_balence_ERM_fungal <- mycofon_balence_ERM_N_fungal <- NULL
  mycofon_balence_NC_in_fungai_opt_roots <- mycofon_balence_NC_in_fungai_opt_N_roots <- mycofon_balence_NC_in_fungai_opt_fungal <- mycofon_balence_NC_in_fungai_opt_N_fungal <- NULL
  mycofon_balence_turnover_roots_roots <- mycofon_balence_turnover_roots_N_roots <- mycofon_balence_turnover_roots_fungal <- mycofon_balence_turnover_roots_N_fungal <- NULL
  mycofon_balence_turnover_mycorrhized_roots <- mycofon_balence_turnover_mycorrhized_N_roots <- mycofon_balence_turnover_mycorrhized_fungal <- mycofon_balence_turnover_mycorrhized_N_fungal <- NULL
  mycofon_balence_turnover_mantle_roots <- mycofon_balence_turnover_mantle_N_roots <- mycofon_balence_turnover_mantle_fungal <- mycofon_balence_turnover_mantle_N_fungal <- NULL
  mycofon_balence_turnover_ERM_roots <- mycofon_balence_turnover_ERM_N_roots <- mycofon_balence_turnover_ERM_fungal <- mycofon_balence_turnover_ERM_N_fungal <- NULL

  ####
  # vector generation
  ####

  temp_range = -20:20
  count = 1
  for (tem in temp_range) {
    mycofon_balence_T_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                      C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                      c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, tem, Tsb, SWC,
                                                      N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                      mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                      max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_Tsb_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, tem, SWC,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots

    mycofon_balence_T_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, tem, Tsb, SWC,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_Tsb_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                          C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                          c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, tem, SWC,
                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                          mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots

    mycofon_balence_T_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                       C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                       c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, tem, Tsb, SWC,
                                                       N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                       mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                       max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_Tsb_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, tem, SWC,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal

    mycofon_balence_T_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, tem, Tsb, SWC,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal
    mycofon_balence_Tsb_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                           C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                           c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, tem, SWC,
                                                           N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                           mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                           max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal
    count = count + 1
  }
  concentration_range = seq(0.01, 1500, by = 1)
  count = 1
  for (conc in concentration_range) {
    mycofon_balence_C_roots_roots[count] <- mycofon_balence(conc, N_roots, optimal_root_fungal_biomass_ratio,
                                                            C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                            c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                            N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                            mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                            max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_C_roots_N_roots[count] <- mycofon_balence(conc, N_roots, optimal_root_fungal_biomass_ratio,
                                                              C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                              c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                              N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                              mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                              max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_C_roots_fungal[count] <- mycofon_balence(conc, N_roots, optimal_root_fungal_biomass_ratio,
                                                             C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                             c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                             N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                             mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                             max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_C_roots_N_fungal[count] <- mycofon_balence(conc, N_roots, optimal_root_fungal_biomass_ratio,
                                                               C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                               c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                               N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                               mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                               max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_C_fungal_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                             conc, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                             c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                             N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                             mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                             max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_C_fungal_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                               conc, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                               c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                               N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                               mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                               max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_C_fungal_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                              conc, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                              c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                              N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                              mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                              max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_C_fungal_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                conc, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_N_roots_roots[count] <- mycofon_balence(C_roots, conc, optimal_root_fungal_biomass_ratio,
                                                            C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                            c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                            N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                            mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                            max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_N_roots_N_roots[count] <- mycofon_balence(C_roots, conc, optimal_root_fungal_biomass_ratio,
                                                              C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                              c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                              N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                              mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                              max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_N_roots_fungal[count] <- mycofon_balence(C_roots, conc, optimal_root_fungal_biomass_ratio,
                                                             C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                             c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                             N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                             mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                             max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_N_roots_N_fungal[count] <- mycofon_balence(C_roots, conc, optimal_root_fungal_biomass_ratio,
                                                               C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                               c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                               N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                               mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                               max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_N_fungal_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                             C_fungal, conc, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                             c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                             N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                             mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                             max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_N_fungal_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                               C_fungal, conc, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                               c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                               N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                               mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                               max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_N_fungal_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                              C_fungal, conc, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                              c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                              N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                              mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                              max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_N_fungal_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                C_fungal, conc, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_NH4_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), conc, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_NH4_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                          C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                          c(Rm, Q10), conc, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                          mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_NH4_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), conc, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_NH4_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                           C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                           c(Rm, Q10), conc, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                           N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                           mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                           max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_NO3_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), NH4, conc, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_NO3_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                          C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                          c(Rm, Q10), NH4, conc, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                          mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_NO3_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), NH4, conc, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_NO3_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                           C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                           c(Rm, Q10), NH4, conc, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                           N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                           mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                           max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_FOM_Norg_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                             C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                             c(Rm, Q10), NH4, NO3, conc, NC_in_fungai_opt, temp, Tsb, SWC,
                                                             N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                             mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                             max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_FOM_Norg_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                               C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                               c(Rm, Q10), NH4, NO3, conc, NC_in_fungai_opt, temp, Tsb, SWC,
                                                               N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                               mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                               max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_FOM_Norg_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                              C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                              c(Rm, Q10), NH4, NO3, conc, NC_in_fungai_opt, temp, Tsb, SWC,
                                                              N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                              mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                              max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_FOM_Norg_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                c(Rm, Q10), NH4, NO3, conc, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_mm_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                       C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                       c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                       N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                       conc, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                       max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_mm_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         conc, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_mm_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        conc, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_mm_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                          C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                          c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                          conc, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_ERM_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        mantle_mass, conc, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_ERM_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                          C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                          c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                          mantle_mass, conc, parameters_NH4_on_NO3, growth_C, growth_N,
                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_ERM_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         mantle_mass, conc, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_ERM_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                           C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                           c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                           N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                           mantle_mass, conc, parameters_NH4_on_NO3, growth_C, growth_N,
                                                           max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal
    count = count + 1
  }
  SWC_range = seq(0.01, 1, by = 0.01)
  count = 1
  for (SoilWater in SWC_range) {
    mycofon_balence_SWC_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                        c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SoilWater,
                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                        mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_SWC_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                          C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                          c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SoilWater,
                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                          mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_SWC_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SoilWater,
                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_SWC_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                           C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                           c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SoilWater,
                                                           N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                           mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                           max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_optimal_root_fungal_biomass_ratio_roots[count] <- mycofon_balence(C_roots, N_roots, SoilWater,
                                                                                      C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                                      c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                                      N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                                      mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                                      max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_optimal_root_fungal_biomass_ratio_N_roots[count] <- mycofon_balence(C_roots, N_roots, SoilWater,
                                                                                        C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                                        c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                                        N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                                        mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                                        max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_optimal_root_fungal_biomass_ratio_fungal[count] <- mycofon_balence(C_roots, N_roots, SoilWater,
                                                                                       C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                                       c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                                       N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                                       mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                                       max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_optimal_root_fungal_biomass_ratio_N_fungal[count] <- mycofon_balence(C_roots, N_roots, SoilWater,
                                                                                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_turnover_roots_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                   C_fungal, N_fungal, SoilWater, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                   c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                   N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                   mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                   max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_turnover_roots_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                     C_fungal, N_fungal, SoilWater, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                     c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                     N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                     mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                     max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_turnover_roots_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                    C_fungal, N_fungal, SoilWater, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                    c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                    N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                    mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                    max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_turnover_roots_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                      C_fungal, N_fungal, SoilWater, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                                                                      c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                      N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                      mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                      max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_turnover_mycorrhized_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                         C_fungal, N_fungal, turnover_roots, SoilWater, turnover_mantle, turnover_ERM,
                                                                         c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                         N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                         mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                         max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_turnover_mycorrhized_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                           C_fungal, N_fungal, turnover_roots, SoilWater, turnover_mantle, turnover_ERM,
                                                                           c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                           N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                           mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                           max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_turnover_mycorrhized_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                          C_fungal, N_fungal, turnover_roots, SoilWater, turnover_mantle, turnover_ERM,
                                                                          c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                          N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                          mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                          max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_turnover_mycorrhized_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                            C_fungal, N_fungal, turnover_roots, SoilWater, turnover_mantle, turnover_ERM,
                                                                            c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                            N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                            mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                            max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_turnover_mantle_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                    C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, SoilWater, turnover_ERM,
                                                                    c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                    N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                    mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                    max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_turnover_mantle_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                      C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, SoilWater, turnover_ERM,
                                                                      c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                      N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                      mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                      max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_turnover_mantle_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                     C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, SoilWater, turnover_ERM,
                                                                     c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                     N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                     mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                     max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_turnover_mantle_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                       C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, SoilWater, turnover_ERM,
                                                                       c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                       N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                       mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                       max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    mycofon_balence_turnover_ERM_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                    C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, SoilWater,
                                                                    c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                    N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                    mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                    max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_roots
    mycofon_balence_turnover_ERM_N_roots[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                   C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, SoilWater,
                                                                   c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                   N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                   mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                   max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_roots
    mycofon_balence_turnover_ERM_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                  C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, SoilWater,
                                                                  c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                  N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                  mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                  max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$C_fungal
    mycofon_balence_turnover_ERM_N_fungal[count] <- mycofon_balence(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                                                                    C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, SoilWater,
                                                                    c(Rm, Q10), NH4, NO3, FOM, NC_in_fungai_opt, temp, Tsb, SWC,
                                                                    N_limits_Plant, N_k_Plant, SWC_k_Plant, N_limits_Fungal, N_k_Fungal, SWC_k_Fungal,
                                                                    mantle_mass, ERM_mass, parameters_NH4_on_NO3, growth_C, growth_N,
                                                                    max_C_allocation_CASSIA, to_CASSIA, mycofon_stratergy)$N_fungal

    count = count + 1
  }

  ####
  # Plots!
  ####

  par(mfrow = c(3, 2))
  axis_1 = c(min(mycofon_balence_T_roots, mycofon_balence_T_N_roots, mycofon_balence_T_fungal, mycofon_balence_T_N_fungal),
             max(mycofon_balence_T_roots, mycofon_balence_T_N_roots, mycofon_balence_T_fungal, mycofon_balence_T_N_fungal))
  plot(temp_range, mycofon_balence_T_roots, main = "kg C", xlab = "Temperature Range: 'C", ylab = "C kg eq", ylim = axis_1)
  points(temp_range, mycofon_balence_T_N_roots, col = "blue")
  points(temp_range, mycofon_balence_T_fungal, pch = "+")
  points(temp_range, mycofon_balence_T_N_fungal, pch = "+", col = "blue")
  legend(temp_range[1], axis_1[2], c("Roots C", "Fungal C", "Roots N", "Fungal N"), col = c("black", "black", "blue", "blue"),  pch = c("o", "+", "o", "+"), bty = "n")

  axis_2 = c(min(mycofon_balence_Tsb_roots, mycofon_balence_Tsb_N_roots, mycofon_balence_Tsb_fungal, mycofon_balence_Tsb_N_fungal),
             max(mycofon_balence_Tsb_roots, mycofon_balence_Tsb_N_roots, mycofon_balence_Tsb_fungal, mycofon_balence_Tsb_N_fungal))
  plot(temp_range, mycofon_balence_Tsb_roots, main = "kg C eq", xlab = "B horizon temperature Range: 'C", ylab = "C kg eq", ylim = axis_2)
  points(temp_range, mycofon_balence_Tsb_N_roots, col = "blue")
  points(temp_range, mycofon_balence_Tsb_fungal, pch = "+")
  points(temp_range, mycofon_balence_Tsb_N_fungal, pch = "+", col = "blue")

  axis_3 = c(min(mycofon_balence_SWC_roots, mycofon_balence_SWC_N_roots, mycofon_balence_SWC_fungal, mycofon_balence_SWC_N_fungal),
             max(mycofon_balence_SWC_roots, mycofon_balence_SWC_N_roots, mycofon_balence_SWC_fungal, mycofon_balence_SWC_N_fungal))
  plot(SWC_range, mycofon_balence_SWC_roots, main = "kg C eq", xlab = "Soil water concentration range: %", ylab = "C kg", ylim = axis_3)
  points(SWC_range, mycofon_balence_SWC_N_roots, col = "blue")
  points(SWC_range, mycofon_balence_SWC_fungal, pch = "+")
  points(SWC_range, mycofon_balence_SWC_N_fungal, pch = "+", col = "blue")

  axis_4 = c(min(mycofon_balence_NH4_roots, mycofon_balence_NH4_N_roots, mycofon_balence_NH4_fungal, mycofon_balence_NH4_N_fungal),
             max(mycofon_balence_NH4_roots, mycofon_balence_NH4_N_roots, mycofon_balence_NH4_fungal, mycofon_balence_NH4_N_fungal))
  plot(concentration_range, mycofon_balence_NH4_roots, main = "kg C", xlab = "Concentration range: NH4", ylab = "C kg", ylim = axis_4)
  points(concentration_range, mycofon_balence_NH4_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_NH4_fungal, pch = "+")
  points(concentration_range, mycofon_balence_NH4_N_fungal, pch = "+", col = "blue")

  axis_5 = c(min(mycofon_balence_NO3_roots, mycofon_balence_NO3_N_roots, mycofon_balence_NO3_fungal, mycofon_balence_NO3_N_fungal),
             max(mycofon_balence_NO3_roots, mycofon_balence_NO3_N_roots, mycofon_balence_NO3_fungal, mycofon_balence_NO3_N_fungal))
  plot(concentration_range, mycofon_balence_NO3_roots, main = "kg C", xlab = "Concentration range: NO3", ylab = "C kg", ylim = axis_5)
  points(concentration_range, mycofon_balence_NO3_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_NO3_fungal, pch = "+")
  points(concentration_range, mycofon_balence_NO3_N_fungal, pch = "+", col = "blue")

  axis_6 = c(min(mycofon_balence_FOM_Norg_roots, mycofon_balence_FOM_Norg_N_roots, mycofon_balence_FOM_Norg_fungal, mycofon_balence_FOM_Norg_N_fungal),
             max(mycofon_balence_FOM_Norg_roots, mycofon_balence_FOM_Norg_N_roots, mycofon_balence_FOM_Norg_fungal, mycofon_balence_FOM_Norg_N_fungal))
  plot(concentration_range, mycofon_balence_FOM_Norg_roots, main = "kg C", xlab = "Concentration range: FOM_Norg", ylab = "C kg", ylim = axis_6)
  points(concentration_range, mycofon_balence_FOM_Norg_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_FOM_Norg_fungal, pch = "+")
  points(concentration_range, mycofon_balence_FOM_Norg_N_fungal, pch = "+", col = "blue")

  par(mfrow = c(2, 2))
  axis_7 = c(min(mycofon_balence_turnover_roots_roots, mycofon_balence_turnover_roots_N_roots, mycofon_balence_turnover_roots_fungal, mycofon_balence_turnover_roots_N_fungal),
             max(mycofon_balence_turnover_roots_roots, mycofon_balence_turnover_roots_N_roots, mycofon_balence_turnover_roots_fungal, mycofon_balence_turnover_roots_N_fungal))
  plot(SWC_range, mycofon_balence_turnover_roots_roots, main = "kg C", xlab = "Root turnover rate: %", ylab = "C kg", ylim = axis_7)
  points(SWC_range, mycofon_balence_turnover_roots_N_roots, col = "blue")
  points(SWC_range, mycofon_balence_turnover_roots_fungal, pch = "+")
  points(SWC_range, mycofon_balence_turnover_roots_N_fungal, pch = "+", col = "blue")

  axis_8 = c(min(mycofon_balence_turnover_mycorrhized_roots, mycofon_balence_turnover_mycorrhized_N_roots, mycofon_balence_turnover_mycorrhized_fungal, mycofon_balence_turnover_mycorrhized_N_fungal),
             max(mycofon_balence_turnover_mycorrhized_roots, mycofon_balence_turnover_mycorrhized_N_roots, mycofon_balence_turnover_mycorrhized_fungal, mycofon_balence_turnover_mycorrhized_N_fungal))
  plot(SWC_range, mycofon_balence_turnover_mycorrhized_roots, main = "kg C", xlab = "Mycorrhized root turnover rate: %", ylab = "C kg", ylim = axis_8)
  points(SWC_range, mycofon_balence_turnover_mycorrhized_N_roots, col = "blue")
  points(SWC_range, mycofon_balence_turnover_mycorrhized_fungal, pch = "+")
  points(SWC_range, mycofon_balence_turnover_mycorrhized_N_fungal, pch = "+", col = "blue")

  axis_9 = c(min(mycofon_balence_turnover_ERM_roots, mycofon_balence_turnover_ERM_N_roots, mycofon_balence_turnover_ERM_fungal, mycofon_balence_turnover_ERM_N_fungal),
             max(mycofon_balence_turnover_ERM_roots, mycofon_balence_turnover_ERM_N_roots, mycofon_balence_turnover_ERM_fungal, mycofon_balence_turnover_ERM_N_fungal))
  plot(SWC_range, mycofon_balence_turnover_ERM_roots, main = "kg C", xlab = "ERM turnover rate: %", ylab = "C kg", ylim = axis_9)
  points(SWC_range, mycofon_balence_turnover_ERM_N_roots, col = "blue")
  points(SWC_range, mycofon_balence_turnover_ERM_fungal, pch = "+")
  points(SWC_range, mycofon_balence_turnover_ERM_N_fungal, pch = "+", col = "blue")

  axis_9_5 = c(min(mycofon_balence_turnover_mantle_roots, mycofon_balence_turnover_mantle_N_roots, mycofon_balence_turnover_mantle_fungal, mycofon_balence_turnover_mantle_N_fungal),
             max(mycofon_balence_turnover_mantle_roots, mycofon_balence_turnover_mantle_N_roots, mycofon_balence_turnover_mantle_fungal, mycofon_balence_turnover_mantle_N_fungal))
  plot(SWC_range, mycofon_balence_turnover_mantle_roots, main = "kg C", xlab = "Mantle turnover rate: %", ylab = "C kg", ylim = axis_9_5)
  points(SWC_range, mycofon_balence_turnover_mantle_N_roots, col = "blue")
  points(SWC_range, mycofon_balence_turnover_mantle_fungal, pch = "+")
  points(SWC_range, mycofon_balence_turnover_mantle_N_fungal, pch = "+", col = "blue")

  par(mfrow = c(2, 2))
  axis_10 = c(min(mycofon_balence_C_roots_roots, mycofon_balence_C_roots_N_roots, mycofon_balence_C_roots_fungal, mycofon_balence_C_roots_N_fungal, na.rm = T),
              max(mycofon_balence_C_roots_roots, mycofon_balence_C_roots_N_roots, mycofon_balence_C_roots_fungal, mycofon_balence_C_roots_N_fungal, na.rm = T))
  plot(concentration_range, mycofon_balence_C_roots_roots, main = "kg C", xlab = "Concentration range: C roots", ylab = "C kg", ylim = axis_10)
  points(concentration_range, mycofon_balence_C_roots_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_C_roots_fungal, pch = "+")
  points(concentration_range, mycofon_balence_C_roots_N_fungal, pch = "+", col = "blue")
  legend(temp_range[1], axis_10[2], c("Roots C", "Fungal C", "Roots N", "Fungal N"), col = c("black", "black", "blue", "blue"),  pch = c("o", "+", "o", "+"), bty = "n")

  axis_11 = c(min(mycofon_balence_N_roots_roots, mycofon_balence_N_roots_N_roots, mycofon_balence_N_roots_fungal, mycofon_balence_N_roots_N_fungal, na.rm = T),
              max(mycofon_balence_N_roots_roots, mycofon_balence_N_roots_N_roots, mycofon_balence_N_roots_fungal, mycofon_balence_N_roots_N_fungal, na.rm = T))
  plot(concentration_range, mycofon_balence_N_roots_roots, main = "kg C", xlab = "Concentration range: N roots", ylab = "C kg", ylim = axis_11)
  points(concentration_range, mycofon_balence_N_roots_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_N_roots_fungal, pch = "+")
  points(concentration_range, mycofon_balence_N_roots_N_fungal, pch = "+", col = "blue")

  axis_12 = c(min(mycofon_balence_C_fungal_roots, mycofon_balence_C_fungal_N_roots, mycofon_balence_C_fungal_fungal, mycofon_balence_C_fungal_N_fungal, na.rm = T),
              max(mycofon_balence_C_fungal_roots, mycofon_balence_C_fungal_N_roots, mycofon_balence_C_fungal_fungal, mycofon_balence_C_fungal_N_fungal, na.rm = T))
  plot(concentration_range, mycofon_balence_C_fungal_roots, main = "kg C", xlab = "Concentration range: C fungal", ylab = "C kg", ylim = axis_12)
  points(concentration_range, mycofon_balence_C_fungal_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_C_fungal_fungal, pch = "+")
  points(concentration_range, mycofon_balence_C_fungal_N_fungal, pch = "+", col = "blue")

  axis_13 = c(min(mycofon_balence_N_fungal_roots, mycofon_balence_N_fungal_N_roots, mycofon_balence_N_fungal_fungal, mycofon_balence_N_fungal_N_fungal, na.rm = T),
              max(mycofon_balence_N_fungal_roots, mycofon_balence_N_fungal_N_roots, mycofon_balence_N_fungal_fungal, mycofon_balence_N_fungal_N_fungal, na.rm = T))
  plot(concentration_range, mycofon_balence_N_fungal_roots, main = "kg C", xlab = "Concentration range: N fungal", ylab = "C kg", ylim = axis_13)
  points(concentration_range, mycofon_balence_N_fungal_N_roots, col = "blue")
  points(concentration_range, mycofon_balence_N_fungal_fungal, pch = "+")
  points(concentration_range, mycofon_balence_N_fungal_N_fungal, pch = "+", col = "blue")

}
