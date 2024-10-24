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
    parameters_R = parameters_R,
    common = common_p,  # "./data/common.csv",
    sperling = sperling_p,
    repo = repo_p,
    pPREL = c(413.0, 0.450, 0.118, 3.0, 0.748464, 12.74915, -3.566967, 18.4513, -0.136732,
               0.033942, 0.448975, 0.500, -0.364, 0.33271, 0.857291, 0.041781,
               0.474173, 0.278332, 1.5, 0.33, 4.824704, 0.0, 0.0, 180.0, 0.0, 0.0, 10.0,
               -999.9, -999.9, -999.9, 1/0.012, 0.0),

    #####
    ## Default values of the set up
    #####

    settings = settings_basic,

    s.D0 = 79,					# DOY to start the calculation of temperature sum, 1=Jan 1; 69=March 1; 79=March 20 for diameter growth. Valid for Finland
    s.H0 = 1,					# and for shoot grwoth

    growth_photo_coef = 1,
    needle_mass_in = 4.467638,
    Throughfall = 1,
    trenching_year = NA,
    soil = FALSE) {

  ####
  # Testing the input and initial conditions
  ####

  # Is the site is the known sites?
  validate_site(site)

  # Are the model settings valid?
  # If not update them
  updated_settings <- update_model_settings(settings_basic)

  # Is the weather data correct?
  validate_weather_data(weather, settings_basic$PRELES_GPP)

  #####
  ## Model conditions derived from model inputs
  #####
  # years from weather data
  start_year <- as.numeric(substring(weather$date[end(weather$date)], 1, 4))[2]
  end_year <- as.numeric(substring(weather$date[end(weather$date)], 1, 4))[1]

  if (soil) {
    out <- CASSIA_soil(start_year, end_year, weather, GPP_ref,
                         pPREL, t(parameters), common, t(ratios), t(sperling), parameters_R, # site,
                         needle_mass_in,
                         Throughfall, trenching_year,
                         updated_settings)

  } else {
    out <- CASSIA_yearly(start_year, end_year, weather, GPP_ref,
                         pPREL, t(parameters), common, t(ratios), t(sperling), # site,
                         needle_mass_in,
                         Throughfall,
                         updated_settings)
  }

  return(out)
}
