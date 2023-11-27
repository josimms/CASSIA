likelyhood_nitrogen_uptake <- function(par, sum = T) {
  ### Import data
  data_path = "~/Downloads/"
  data = read.delim(paste0(data_path, "oyewole_2015_calibration_data.csv"), sep = ",", dec = ",")

  ### Simulations
  concentrations = rep(c(50, 500), length.out = 16)
  glycine <- nh4 <- no3 <- NULL
  for (i in 1:16) {
    glycine[i] = uptake_N(concentrations[i], 25, par[1], par[2], 0.52, 0.4)
    nh4[i] = uptake_N(concentrations[i], 25, par[3], par[4], 0.52, 0.4)
    no3[i] = uptake_N(concentrations[i], 25, par[5], par[6], 0.52, 0.4)
  }
  simVec = c(glycine, nh4, no3)
  # Temperature change with season
  # "There was no interaction effect of fertilisation and season n root uptake rates"
  # Temperature 25 as this means there is no temperature effect
  # 0.52 taken as the SWC concentration as in the paper it is said to be 0.47-0.57

  ### Error
  # Data error is linear with data value
  sdX = simVec*par[7]

  ### Residuals
  res = data$control[-c(1, 2, 9, 10)] - simVec

  ###
  likeX <- sum(dhtn(res, mean = 0, sd = sdX, log = T))

  return(likeX)
}

nitrogen_uptake_calibration <- function() {
  file_path = "~/Downloads/"
  almost_zero = 1e16
  bounds_lower = rep(almost_zero, 7)
  bounds_higher = c(rep(15, 6), 3)

  prior = BayesianTools::createUniformPrior(bounds_higher, bounds_lower, best = NULL)

  BayesianSetiup = BayesianTools::createBayesianSetup(likelihood = likelyhood,
                                                      prior = prior,
                                                      names = c("amino_acid_limit", "amino_acid_k", "NH4_limit", "NH4_k", "MO3_limit", "NO3_k", "variation"),
                                                      parallel = F)

  settings = list(iterations = 1e6, thin = 100, nrChains = 1, message = T)

  out_uptake_n = BayesianTools::runMCMC(bayesianSetup = BayesianSetiup, sampler = "DREAMzs", settings = settings)
  # save(out_uptake_n, file = paste0(file_path, gsub(":", "_", Sys.time()), "_uptake_n"))

  summary(out_uptake_n)
  params_cali = BayesianTools::MAP(out_uptake_n)$parametersMAP

  Temp = 25
  SoilWater = 0.52
  SWC_sat = 0.4
  nitrogen_graphs(50, Temp, params_cali[1], params_cali[2], SoilWater, SWC_sat)
  nitrogen_graphs(50, Temp, params_cali[3], params_cali[4], SoilWater, SWC_sat)
  nitrogen_graphs(50, Temp, params_cali[5], params_cali[6], SoilWater, SWC_sat)
}

