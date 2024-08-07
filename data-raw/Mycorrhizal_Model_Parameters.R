###
# Parameters
###

### Root data
microbe_turnover = 0.016906 # Preeven 2013, Symphony model paper, rate of microbes going to SOM
turnover_mantle = 1/625 # Meyer, 2010 ERM turnover ranges between 10 days to 100 days, while mantle turnover time is 625 days
turnover_ERM = 1/50
turnover_roots = 1/365 # Meyer, 2010 Non-mycorrhizal roots have a turnover (dr_r) of 365 days
turnover_roots_mycorrhized = 1/625 # Meyer, 2010

### Optimal ratios assumed to be the average reported in Korhonen 2007
NC_in_root_opt =  200 ###(0.001*2875)/1881 # Helmisaari 1995: Nutroent cycling in Pinus sylvestris stands in eastern Finland, belowground component
NC_in_fungai_opt = 200 # TODO: find this!
NC_microbe_opt = 1 ## 28.73
# Average C/N ratio is 28.73, 1.28% N molecular content - Hyytälä average - Heinonsalo, 2015
# The stable SOM C:N ratio is about 1/6, if it is assumed that proteins are the main cause - Adamczyk, 2019
# Högberg 2017, the C/N ratio of microorganisms in the F + H layer is around 5
NC_fungal_opt = 120  ## 0.025 # Meyer 2010, CN ratio of 40 fungal
NC_Litter = 200

###
percentage_C_biomass = 0.5 # Used in CASSIA
optimal_root_fungal_biomass_ratio = 0.9 # TODO: Heinonsalo
growth_C = 0.3 # Franklin 2017
growth_N = 0.3 # Arbitrary - TODO
Cmic_Corg = 0.5 # TODO: find this value from hyytiälä

imobilisation = 0.1 # TODO:
assimilation = 0.1 # TODO:

### Parameters that link to CASSIA and will be replaced when the models are combined
max_C_allocation_CASSIA = 0.1 # TODO: link this to CASSIA and find a sensible value
to_CASSIA = 0.1 # TODO: link this to CASSIA and find a sensible value
Rm = 0.00958
Q10 = 2.5575


N_limits_plant = c(10, 10, 10) # Arbitrarty - Bayesian Callibration
N_limits_fungal = c(10, 10, 10) # Arbitrarty - Bayesian Callibration
N_limits_microbes = c(10, 10, 10, 10, 10) # Arbitrarty - Bayesian Callibration
N_k_plant = c(0.2, 0.2, 0.2) # Arbitrarty - Bayesian Callibration
N_k_fungal = c(0.2, 0.2, 0.2) # Arbitrarty - Bayesian Callibration
N_k_microbes = c(0.2, 0.2, 0.2, 0.2, 0.2) # Arbitrarty - Bayesian Callibration
SWC_sat_plant = c(0.4, 0.4, 0.4) # 0.4 value comes from the methods in Oyewole 2015, rest Arbitrarty - Bayesian Callibration
SWC_sat_fungal = c(0.5, 0.5, 0.5) # Arbitrarty - Bayesian Callibration
SWC_sat_microbes = c(0.5, 0.5, 0.5, 0.5, 0.5) # Arbitrarty - Bayesian Callibration
NH4_on_NO3 = c(10, 0.2, 0.5) # Arbitrarty - Bayesian Callibration

C_value_param_myco = 1 # TODO
N_value_param_myco = 1 # TODO
C_value_param_plant = 1 # TODO
N_value_param_plant = 1 # TODO

# TODO: check these!
parameters_R = c(microbe_turnover, NC_in_root_opt, NC_fungal_opt, NC_microbe_opt, percentage_C_biomass,
                 N_limits_plant, N_limits_fungal, N_limits_microbes,
                 N_k_plant, N_k_fungal, N_k_microbes,
                 SWC_sat_plant, SWC_sat_fungal, SWC_sat_microbes, NH4_on_NO3,
                 c(Rm, Q10, Rm, Q10, Rm, Q10),
                 optimal_root_fungal_biomass_ratio, turnover_mantle, turnover_ERM, turnover_roots, turnover_roots_mycorrhized,
                 growth_C, growth_N, C_value_param_myco, N_value_param_myco, C_value_param_plant, N_value_param_plant)

save(parameters_R, file = "./data/parameters_R.RData")

###
# State variables
###

Temp = 15
Tsb = 5
SoilWater = 0.153 # Sensible value from Hyytiälä tracker data

C_roots = 2.8/2 # CASSIA value from the starch and sugar 15_16 file
N_roots = C_roots * NC_in_root_opt # As initial condition assume that the values are at an optimum
C_fungal = C_roots * optimal_root_fungal_biomass_ratio
N_fungal = C_fungal * NC_fungal_opt
Norg_FOM = 1500 # TODO: this should be an input, make sure that the Norg and Norg_FOM are used consistently
micorization = N_roots / (C_roots * optimal_root_fungal_biomass_ratio)
NH4 = 0.31 # kg N ha−1 Korhonen, 2013
NO3 = 0.002 # kg N ha−1 Korhonen, 2013
Norg = 26.5 # kg N ha−1 Korhonen, 2013
C_SOM = 5.6 # kg m-2 C in 40cm hummus layer measured 2003 - Hyytiälä site description
C_microbe = C_SOM * Cmic_Corg
N_microbe = C_microbe * NC_microbe_opt
C_fungal_in = C_SOM # This is just for the graphs, in the model the value should be calculated by uptake
N_fungal_in = NH4 + NO3 + Norg # This is just for the graphs, in the model the value should be calculated by uptake
mantle_mass = 0.5*C_roots # TODO: what values should these have?
ERM_mass = 0.5*C_roots

C_FOM_needles = 1 # TODO: find a real value for this!
C_FOM_woody = 1 # TODO: find a real value for this!
C_FOM_roots = 1 # TODO: find a real value for this!
C_FOM_mantle = 1 # TODO: find a real value for this!
C_FOM_ERM = 1 # TODO: find a real value for this!
C_FOM = C_FOM_needles + C_FOM_woody + C_FOM_roots + C_FOM_mantle + C_FOM_ERM
C_SOM = 1 # TODO: find a real value for this!
N_SOM = 1 # TODO: find a real value for this!
C_decompose_FOM = 1 # TODO: find a real value for this!
C_decompose_SOM = 1 # TODO: find a real value for this!
N_decompose_FOM = 1 # TODO: find a real value for this!
N_decompose_SOM = 1 # TODO: find a real value for this!
Litter_needles = 1 # TODO: find a real value for this!
Litter_woody = 1 # TODO: find a real value for this!
Litter_roots = 1 # TODO: find a real value for this!
Litter_mycelium = 1 # TODO: find a real value for this!
N_FOM_needles = 1 # TODO: find a real value for this!
N_FOM_woody = 1 # TODO: find a real value for this!
N_FOM_roots = 1 # TODO: find a real value for this!
N_FOM_mycelium = 1 # TODO: find a real value for this!

NH4_used_Plant = 1 # TODO: This would be an output of the other models!
NH4_used_Fungal = 1 # TODO: This would be an output of the other models!
NO3_used_Plant = 1 # TODO: This would be an output of the other models!
NO3_used_Fungal = 1 # TODO: This would be an output of the other models!
FOM_Norg_used_Plant = 1 # TODO: This would be an output of the other models!
FOM_Norg_used_Fungal = 1 # TODO: This would be an output of the other models!
SOM_Norg_used = 1 # TODO: This would be an output of the other models!


