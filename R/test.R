tests <- function() {
  ###
  # Check that the outputs are providing the right behaviour
  ###
  # TODO: make the graph titles match the input variables

  respiration_graphs(Rm, Q10)

  ## NH4 version
  NH4 = 0.3
  NO3 = 0.03
  Norg = 0.3
  Temp = 15
  Tsb = 5
  SoilWater = 0.5
  N_limits_plant = c(0.3, 0.3, 0.3)
  N_k_plant = c(10, 10, 10)
  SWC_sat_plant = c(0.2, 0.2, 0.2)
  nitrogen_graphs(NH4, Temp, SoilWater, N_limits_plant[1], N_k_plant[1], SWC_sat_plant[1])

  C_roots = 0.03
  N_roots = 0.03
  C_fungal = 0.03
  N_fungal = 0.03
  optimal_root_fungal_biomass_ratio = 0.5
  NC_fungal_opt = 0.5
  decision_graphs(C_roots, N_roots, C_fungal, N_fungal,
                  optimal_root_fungal_biomass_ratio,
                  NC_fungal_opt)

  growth_C = 0.3
  growth_N = 0.3
  C_fungal_in = 5
  N_fungal_in = 5
  myco_growth_graphs(C_fungal_in, N_fungal_in, growth_C, growth_N)

  micorization = 0.9
  demand = 1
  NH4_on_NO3 = 0.3
  nitrogen_graphs_plant(Temp, SoilWater, micorization, NH4, NO3, Norg,
                        N_limits_plant, N_k_plant, SWC_sat_plant, NH4_on_NO3, demand)

  N_limits_fungal = c(0.3, 0.3, 0.3)
  N_k_fungal = c(10, 10, 10)
  SWC_sat_fungal = c(0.2, 0.2, 0.2)
  nitrogen_graphs_fungal(Temp, SoilWater,
                         NH4, NO3, Norg,
                         N_limits_fungal, N_k_fungal,
                         SWC_sat_fungal, demand)
  # TODO: not sure that this is doing the right thing! Could be the graphs that are being printed!

  C_microbe = 26
  N_microbe = 26
  C_SOM = 50
  C_FOM = 26
  NC_microbe_opt = 0.5
  NC_Litter = 0.8
  imobilisation = 0.2
  assimilation = 0.2
  N_limits_microbes = c(0.3, 0.3, 0.3)
  N_k_microbes = c(10, 10, 10)
  SWC_sat_microbes = c(0.2, 0.2, 0.2)
  nitrogen_graphs_microbe(C_microbe, N_microbe, C_SOM, NC_microbe_opt, NH4, NO3, Norg,
                          Temp, SoilWater, NC_Litter, imobilisation, assimilation,
                          N_limits_microbes, N_k_microbes, SWC_sat_microbes,
                          TRUE)
  nitrogen_graphs_microbe(C_microbe, N_microbe, C_FOM, NC_microbe_opt, NH4, NO3, Norg,
                          Temp, SoilWater, NC_Litter, imobilisation, assimilation,
                          N_limits_microbes, N_k_microbes, SWC_sat_microbes,
                          FALSE)


  parameters_R = c(0.016906, # microbe_turnover (Preveen, 2013)
                   (0.001*2875)/1881, # NC_in_root_opt (Heimisaari, 1995)
                   0.025, # NC_fungal_opt (Meyer, 2010)
                   1/28.73, # NC_microbe_opt (Heinonsalo, 2015)
                   0.5, # percentage_C_biomass (CASSIA)
                   # fungal
                   0.3, 0.3, 0.3, # N_limits: NH4, NO3, Norg
                   50, 50, 50, # N_k: NH4, NO3, Norg
                   0.2, 0.2, 0.2, # SWC_limits: NH4, NO3, Norg
                   # plant
                   0.3, 0.3, 0.3, # N_limits: NH4, NO3, Norg
                   50, 50, 50, # N_k: NH4, NO3, Norg
                   0.2, 0.2, 0.2, # SWC_limit: NH4, NO3, Norg
                   # microbes
                   0.3, 0.3, 0.3, # N_limits: NH4, NO3, Norg
                   50, 50, 50, # N_k: NH4, NO3, Norg
                   0.2, 0.2, 0.2, # SWC_limit: NH4, NO3, Norg
                   0.3, 50, 0.2, # C limits
                   10, # NH4_on_NO3
                   0.9, # optimal_root_fungal_biomass_ratio (TODO: Heinonsalo?)
                   1/625, # turnover_mantle, Meyer, 2010
                   1/50, # turnover_ERM, Meyer 2010
                   1/365, # turnover_roots, Meyer, 2010
                   1/625, # turnover_roots_mycorrhized Meyer, 2010
                   0.2, # turnover_fungal TODO: do I need this if the turnover is in Meyer?
                   1, # mantle_mass
                   1, # ERM_mass
                   0.03, # growth_C (Franklin, 2017)
                   0.03, # growth_N (TODO)
                   1, # C_value_param_myco
                   1, # N_value_param_myco
                   1, # C_value_param_plant
                   1) # N_value_param_plant
  C_roots_NonStruct = 0.03
  N_roots_NonStruct = 0.03
  C_fungal_NonStruct = 0.03
  N_fungal_NonStruct = 0.03
  max_C_allocation_CASSIA = 0.3
  mycofon_balence_graphs(C_roots, N_roots, C_fungal, N_fungal,
                         C_roots_NonStruct, N_roots_NonStruct, C_fungal_NonStruct, N_fungal_NonStruct,
                         max_C_allocation_CASSIA,
                         parameters_R,
                         NH4, NO3, Norg,
                         Temp, Tsb, SoilWater,
                         TRUE)
  mycofon_balence_graphs(C_roots, N_roots, C_fungal, N_fungal,
                         C_roots_NonStruct, N_roots_NonStruct, C_fungal_NonStruct, N_fungal_NonStruct,
                         max_C_allocation_CASSIA,
                         parameters_R,
                         NH4, NO3, Norg,
                         Temp, Tsb, SoilWater,
                         FALSE)

  Litter_mantle = 1
  Litter_ERM = 1
  N_FOM_mantle = 1
  N_FOM_ERM = 1
  NC_needles = 1 # TODO: find a real value for this!
  NC_woody = 1 # TODO: find a real value for this!
  NC_roots = 1 # TODO: find a real value for this!
  NC_mantle = 1 # TODO: find a real value for this!
  NC_ERM = 1 # TODO: find a real value for this!

  symphony_multiple_FOM_daily_plot(Tsb, SoilWater,
                                   C_FOM_needles, C_FOM_woody, C_FOM_roots, C_FOM_mantle, C_FOM_ERM,
                                   C_SOM, N_SOM,
                                   C_decompose_FOM, C_decompose_SOM, N_decompose_FOM, N_decompose_SOM,
                                   Litter_needles, Litter_woody, Litter_roots, Litter_mantle, Litter_ERM,
                                   imobilisation,  assimilation,
                                   NH4, NO3, N_FOM_needles, N_FOM_woody, N_FOM_roots, N_FOM_mantle, N_FOM_ERM,
                                   NH4_used_Plant, NH4_used_Fungal, NO3_used_Plant, NO3_used_Fungal, FOM_Norg_used_Plant, FOM_Norg_used_Fungal,
                                   SOM_Norg_used,
                                   Rm, Q10, N_limits_microbes, N_k_microbes, SWC_sat_microbes,
                                   NC_microbe_opt, microbe_turnover)

}
