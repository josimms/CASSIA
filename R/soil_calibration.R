likelyhood_soil_models <- function(par, sum = T) {
  ### Test data
  data_path = "~/Downloads/CASSIA_Calibration/"
  # Data for the symphony model balance
  mineral_N = read.delim(paste0(data_path, "korhonen_mineral_N.csv"), sep = ",", dec = ",")
  data_mycofon

  ### Input data
  # TODO: add the weather data here

  ### Simulations
  simTable = Toy_Model(year, # TODO: this should also be within the cpp function in the end
                       C_roots, N_roots, C_fungal, N_fungal, Litter_mantle, Litter_ERM, # TODO: these values should be in the function by the end!
                       Hyde_weather, # TODO: check that this goes for the right dates
                       par)
  # TODO: this should include the liquid N pools, NH4, NO3 and DOC pools in the mineral soils for the korhonen data
  simVec_symphony = simTable # TODO: make this a subset
  simVec_mycofon = simTable # TODO: make this a subset

  ### Error
  # TODO: check the data error
  sdX = c(simVec_symphony, simVec_mycofon)*par

  ### Residuals
  res = c(data_symphony, data_mycofon) - c(simVec_symphony, simVec_mycofon)

  ###
  likeX <- sum(dhtn(res, mean = 0, sd = sdX, log = T))

  return(likeX)
}

soil_model_calibration <- function() {
  file_path = "~/Downloads/"
  almost_zero = 1e16
  bounds_lower # TODO: work out which parameters I am fitting and add the bounds here
  bounds_higher
  parameter_names

  prior = BayesianTools::createUniformPrior(bounds_higher, bounds_lower, best = NULL)

  BayesianSetiup = BayesianTools::createBayesianSetup(likelihood = likelyhood_soil_models,
                                                      prior = prior,
                                                      names = parameter_names,
                                                      parallel = F)

  settings = list(iterations = 1e6, thin = 100, nrChains = 1, message = T)

  out_uptake_n = BayesianTools::runMCMC(bayesianSetup = BayesianSetiup, sampler = "DREAMzs", settings = settings)
  # save(out_uptake_n, file = paste0(file_path, gsub(":", "_", Sys.time()), "_uptake_n"))

  summary(out_uptake_n)
  params_cali = BayesianTools::MAP(out_uptake_n)$parametersMAP

  # TODO: make a plot which has the values plotted as a time series
}
