tests <- function() {
  ###
  # Check that the outputs are providing the right behaviour
  ###
  # TODO: make the graph titles match the input variables

  respiration_graphs(Rm, Q10)

  ## NH4 version
  nitrogen_graphs(NH4, Temp, N_limits_plant[1], N_k_plant[1], SoilWater, SWC_sat_plant[1])

  N_allo = 0.1
  max_C_allocation_CASSIA = 0.1
  decision_graphs(C_roots, N_roots, C_fungal, N_fungal,
                  optimal_root_fungal_biomass_ratio,
                  N_allo, max_C_allocation_CASSIA,
                  NC_fungal_opt, growth_C, growth_N)

  myco_growth_graphs(C_fungal_in, N_fungal_in, growth_C, growth_N)

  nitrogen_graphs_plant(Temp, SoilWater, micorization, NH4, NO3, Norg,
                        N_limits_plant, N_k_plant, SWC_sat_plant, c(Rm, Q10), 1)
  nitrogen_graphs_fungal(Temp, SoilWater,
                         NH4, NO3, Norg,
                         N_limits_fungal, N_k_fungal,
                         SWC_sat_fungal, 1)

  nitrogen_graphs_microbe(C_microbe, N_microbe, C_SOM, NC_microbe_opt, NH4, NO3, Norg,
                          Temp, SoilWater, NC_Litter, imobilisation, assimilation,
                          N_limits_microbes, N_k_microbes, SWC_sat_microbes,
                          TRUE, Rm, Q10)
  nitrogen_graphs_microbe(C_microbe, N_microbe, C_FOM, NC_microbe_opt, NH4, NO3, Norg,
                          Temp, SoilWater, NC_Litter, imobilisation, assimilation,
                          N_limits_microbes, N_k_microbes, SWC_sat_microbes,
                          FALSE, Rm, Q10)

  mycofon_balence_graphs(C_roots, N_roots, optimal_root_fungal_biomass_ratio,
                         C_fungal, N_fungal, turnover_roots, turnover_roots_mycorrhized, turnover_mantle, turnover_ERM,
                         Rm, Q10, NH4, NO3, Norg, NC_in_fungai_opt, Temp, Tsb, SoilWater,
                         N_limits_plant, N_k_plant, SWC_sat_plant, N_limits_fungal, N_k_fungal, SWC_sat_fungal,
                         mantle_mass, ERM_mass, NH4_on_NO3, growth_C, growth_N,
                         max_C_allocation_CASSIA, to_CASSIA, TRUE)

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
  ###
  # Check that there is mass conservation
  ###

}
