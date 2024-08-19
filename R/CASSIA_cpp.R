CASSIA_cpp <- function(
    #####
    ## Weather Inputs - input in a dataframe with date, temperature, Photosynthesis, soil temperature a and b horizon, soil moisture and precipitation
    #####

    weather,

    #####
    ## Site
    #####

    site,

    #####
    ## Parameters
    #####
    ratio_sugar = c(11, 1/3, 11, 11), # Determines the concentration difference between organs
    tau.myco = 3,
    tau.t.needles = 3,
    tau.t.phloem = 3,
    tau.t.roots = 3,
    tau.t.xylem.sh = 3,
    tau.t.xylem.st = 3,
    ratios = ratios_p, # "./data/ratio.csv",
    parameters = parameters_p,
    common = common_p,  # "./data/common.csv",
    sperling = sperling_p,
    repo = repo_p,
    pPREL = c(413.0, 0.450, 0.118, 3.0, 0.748464, 12.74915, -3.566967, 18.4513, -0.136732,
               0.033942, 0.448975, 0.500, -0.364, 0.33271, 0.857291, 0.041781,
               0.474173, 0.278332, 1.5, 0.33, 4.824704, 0.0, 0.0, 180.0, 0.0, 0.0, 10.0,
               -999.9, -999.9, -999.9),
    N_parameters = c(1/0.012, 0.0),

    #####
    ## Default values of the set up
    #####

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

    s.D0 = 79,					# DOY to start the calculation of temperature sum, 1=Jan 1; 69=March 1; 79=March 20 for diameter growth. Valid for Finland
    s.H0 = 1,					# and for shoot grwoth

    growth_photo_coef = 1,
    needle_mass_in = 4.467638,
    Throughfall = 1,
    etmodel = F,
    LOGFLAG = F) {
  ####
  # Testing the input and initial conditions
  ####

  # Is the site is the known sites?
  validate_site(site)

  # Are the model settings valid?
  model_settings <- list(
    myco_model = myco_model,
    sperling_model = sperling_model,
    PRELES_GPP = PRELES_GPP,
    mycorrhiza = mycorrhiza,
    phloem_trigger = phloem_trigger,
    xylogenesis = xylogenesis,
    LN_estim = LN_estim,
    trees_grow = trees_grow
  )
  # If not update them
  updated_settings <- update_model_settings(model_settings)
  list2env(updated_settings, envir = environment())

  # Is the weather data correct?
  validate_weather_data(weather, PRELES_GPP)

  #####
  ## Model conditions derived from model inputs
  #####
  # years from weather data
  start_year <- as.numeric(substring(weather$date[end(weather$date)], 1, 4))[2]
  end_year <- as.numeric(substring(weather$date[end(weather$date)], 1, 4))[1]

  out <- CASSIA_yearly(start_year, end_year, weather, GPP_ref,
                       c(pPREL, N_parameters), t(parameters), common, t(ratios), t(sperling), # site,
                       needle_mass_in,
                       Throughfall,
                       storage_reset, storage_grows,
                       LH_estim, LN_estim, mN_varies, LD_estim, sD_estim_T_count,
                       trees_grow, growth_decreases, needle_mass_grows,
                       mycorrhiza, root_as_Ding, sperling_model,
                       xylogenesis, environment_effect_xylogenesis,
                       temp_rise, drought, Rm_acclimation,
                       photosynthesis_as_input, TRUE,
                       etmodel, LOGFLAG)

  return(out)
}
