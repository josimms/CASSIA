CASSIA_cpp <- function(
    #####
    ## Weather Inputs - input in a dataframe with date, temperature, Photosynthesis, soil temperature a and b horizon, soil moisture and precipitation
    #####

    weather = Hyde_weather,

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
    ratios = ratios_p,
    parameters = parameters_p,
    common = common_p,
    sperling = sperling_p,
    repo = repo_p,

    #####
    ## Default values of the set up
    #####

    storage.reset = TRUE,			# storage.reset<-TRUE=Same initial storage each year, storage.reset<-False, The storage on the last day of year X is  postponded to the first day of the year X+1
    storage.grows = FALSE,			# TRUE if the critical storage level increases with tree size.

    LN.estim = TRUE,				# LN depends on the GPP during previous july-august
    mN.varies = TRUE,				# needle mass (in maintenance respiration) is 2/3 of the total during period 1.10 - 31.5.

    LD.estim = TRUE,				# LD depends on the GPP during March-August
    sD.estim.T.count = FALSE,			# sD depends on the number of days when g in growing window - analogue to needles

    LH.estim = TRUE,

    trees_grow = FALSE,				# can be false if mature trees are modelled and not for a very long period
    growth_decreases = FALSE,			# the height and diameter growth (alfa_S and alfaD) decrease during the simulation
    needle_mass_grows = FALSE,		# Is needle mass dynamic i.e. the modelled growth is also respiring etc and following for some years? If true, note that root mass is related to needle mass

    phloem.trigger = FALSE,    # Phloem controls bud burst rather than whole tree sugar

    mychorrhiza = TRUE, 			# If allocation to mychorrhiza is taken into account
    root_as_Ding = TRUE,

    sperling_model = FALSE,       # Dynamic sugar model using Sperling's enzyme dynamics
    myco_model = FALSE,           # Joanna's mycomodel development!
    xylogenesis = FALSE,

    PRELES_GPP = FALSE,
    environment_effect_xylogenesis = FALSE,

    photoparameters = 3,
    temp_rise = FALSE,
    drought = FALSE,
    Rm_acclimation = TRUE,

    s.D0 = 79,					# DOY to start the calculation of temperature sum, 1=Jan 1; 69=March 1; 79=March 20 for diameter growth. Valid for Finland
    s.H0 = 1,					# and for shoot grwoth

    growth_photo_coef = 1
    ) {

  #####
  # The input tests
  #####

  input_tests()

  #####
  # Initial conditions check
  #####

  if (myco_model) {
    if (!sperling_model) {
      sperling_model = TRUE
      warning("sperling_model has been changed to sperling_model = TRUE as the Sperling submodel should control the allocation")
    }
    if (PRELES_GPP) {
      PRELES_GPP == TRUE
      warning("PRELES_GPP has been changed to PRELES_GPP = TRUE as PRELES should control the nitrogen effect")
    }
  }

  if (sperling_model) {
    if (mychorrhiza) {
      mychorrhiza = FALSE
      warning("Mycorrhiza has been changed to mycorrhiza = FALSE as mycorrhiza is included explicitly in the Sperling submodel")
    }
  } else {
    if (phloem.trigger) {phloem.trigger == F}
    warning("phloem.trigger has be set to FALSE as this feature is only possible with the variability given by the sperling_model.")
  }

  if (xylogenesis) {
    # TODO: this is not set in stone, but fits the initial setup of Lettosuo
    warning("As xylogenesis is TRUE, LN.estim, trees_grow, myorrhiza and phloem.trigger set to FALSE, TRUE, FALSE and FALSE respectively")
    LN.estim = FALSE   # LN depends on the GPP during previous july-august
    trees_grow = TRUE  # can be false if mature trees are modelled and not for a very long period
    mycorrhiza = FALSE
    phloem.trigger = FALSE
  }

  #####
  ## Model conditions derived from model inputs
  #####
  # years from weather data
  years <- seq(as.numeric(substring(weather$date[end(weather$date)], 1, 4))[2],
               as.numeric(substring(weather$date[end(weather$date)], 1, 4))[1])
  if (sum(is.na(years)) != 0) {stop("NA in years vector: Check the date format in weather input")}

  #####
  ## Creating the output vectors
  #####

  total.days <- 366 * sum(years %in% leap_years) + 365 * (length(years) - sum(years %in% leap_years))

  ####
  # TODO: reevaluate whether I need this or not when I have coded the Cpp side
  ####

  if (sperling_model == T) {
    export_yearly <- data.frame(matrix(ncol=26, nrow=length(years)))
    export_daily <- data.frame(matrix(ncol=34, nrow=total.days))
    names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                              "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts",
                              "sugar.needles", "sugar.phloem", "sugar.xylem.sh", "sugar.xylem.st", "sugar.roots",
                              "starch.needles", "starch.phloem", "starch.xylem.sh", "starch.xylem.st", "starch.roots")
    names(export_daily) <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                             "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage_term", "to.mycorrhiza", "mycorrhiza.tot",
                             "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot",
                             "sugar.needles", "sugar.phloem", "sugar.xylem.sh", "sugar.xylem.st", "sugar.roots",
                             "starch.needles", "starch.phloem", "starch.xylem.sh", "starch.xylem.st", "starch.roots")
  } else if (xylogenesis == T) {
    export_yearly <- data.frame(matrix(ncol=20, nrow=length(years)))
    export_daily <- data.frame(matrix(ncol=27, nrow=total.days))
    names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                              "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "ring_width", "needle_mass", "sum.needle.cohorts",
                              "ew_width", "ring_density", "ew.cells_tot", "lw.cells_tot")
    names(export_daily) <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                             "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage", "to.mycorrhiza", "mycorrhiza.tot",
                             "P", "to_sugar", "to_starch", "daily.consumption", "ring_width", "GD.tot", "n.E.tot", "n.W.tot", "n.M.tot")
  } else {
    export_yearly <- data.frame(matrix(ncol=16, nrow=length(years)))
    export_daily <- data.frame(matrix(ncol=24, nrow=total.days))
    names(export_yearly) <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                              "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts")
    names(export_daily) <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                             "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage_term", "to.mycorrhiza", "mycorrhiza.tot",
                             "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot")
  }

  n.days.export <- 0
  n.year <- 1

  #####
  ## Non yearly dependent coefficients
  #####
  # height growth coefficient

  height_growth_coefficient <- diameter_growth_coefficient <- NULL
  height_growth_coefficient <- cbind(1997 : 2020, rep(ratios[c("height_growth_coefficient"),c(site)], length.out = length(1997 : 2020)))
  if (sum(is.na(height_growth_coefficient)) != 0) {warning("One height_growth_coefficient value is NA")} # TODO: is this necaisry as I will get a error later anyway
  diameter_growth_coefficient <- cbind(1997 : 2020, rep(ratios[c("diameter_growth_coefficient"),c(site)], length.out = length(1997 : 2020)))
  if (sum(is.na(diameter_growth_coefficient)) != 0) {warning("One diameter_growth_coefficient value is NA")} # TODO: is this necaisry as I will get a error later anyway

  if (growth_decreases == TRUE) {
    height_growth_coefficient <- cbind(1997 : 2020, seq(ratios[c("height_growth_coefficient_max"),c(site)], ratios[c("height_growth_coefficient_min"),c(site)], length.out = length(1997 : 2020)))
    diameter_growth_coefficient <- cbind(1997 : 2020, seq(ratios[c("diameter_growth_coefficient_max"),c(site)], ratios[c("diameter_growth_coefficient_min"),c(site)], length.out = length(1997 : 2020)))
    if (sum(is.na(height_growth_coefficient)) != 0) {warning("One height_growth_coefficient value is NA")} # TODO: is this necaisry as I will get a error later anyway
    if (sum(is.na(diameter_growth_coefficient)) != 0) {warning("One diameter_growth_coefficient value is NA")} # TODO: is this necaisry as I will get a error later anyway
  }

  if (xylogenesis == TRUE) {
    cell.l<-(parameters[c("cell.l.ew"),c(site)]+parameters[c("cell.l.lw"),c(site)])/2
    cell.d<-(parameters[c("cell.d.ew"),c(site)]+parameters[c("cell.d.lw"),c(site)])/2

    cell.volume.ew<-(parameters[c("cell.d.ew"),c(site)])^2*cell.l
    cell.volume.lw<-(parameters[c("cell.d.lw"),c(site)])^2*cell.l

    wall.volume.ew<-cell.volume.ew-(parameters[c("cell.d.ew"),c(site)]-2*parameters[c("wall.thickness.ew"),c(site)])^2*(cell.l-2*parameters[c("wall.thickness.ew"),c(site)])	# m3
    wall.volume.lw<-cell.volume.lw-(parameters[c("cell.d.lw"),c(site)]-2*parameters[c("wall.thickness.lw"),c(site)])^2*(cell.l-2*parameters[c("wall.thickness.lw"),c(site)])	# m3

    cell.wall.volume.growth.per.day.ew <- wall.volume.ew / parameters[c("tau.We"),c(site)]
    cell.wall.volume.growth.per.day.lw <- wall.volume.lw / parameters[c("tau.Wl"),c(site)]
  } else {
    cell.volume.ew<-(parameters[c("cell.d.ew"), c(site)])^2*parameters[c("cell.l.ew"),c(site)]
    cell.volume.lw<-(parameters[c("cell.d.lw"), c(site)])*(parameters[c("cell.d.lw"), c(site)])*parameters[c("cell.l.lw"),c(site)]		# the tangential width of the cell is the same for early and late wood

    wall.volume.ew<-cell.volume.ew-(parameters[c("cell.d.ew"), c(site)]-2*parameters[c("wall.thickness.ew"), c(site)])^2*(parameters[c("cell.l.ew"),c(site)]-2*parameters[c("wall.thickness.ew"), c(site)])	# m3
    wall.volume.lw<-cell.volume.lw-(parameters[c("cell.d.lw"), c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])*(parameters[c("cell.d.lw"), c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])*(parameters[c("cell.l.lw"),c(site)]-2*parameters[c("wall.thickness.lw"), c(site)])	# m3
  }

  if (environment_effect_xylogenesis == FALSE) {
    cell.density.ew<-parameters[c("cell.wall.density.ew"),c(site)]*wall.volume.ew/cell.volume.ew		# to calculate the wood density
    cell.density.lw<-parameters[c("cell.wall.density.lw"),c(site)]*wall.volume.lw/cell.volume.lw

    daily.rate.ew<-wall.volume.ew/parameters[c("tau.We"),c(site)]		#m3 cell wall day-1
    daily.rate.lw<-wall.volume.lw/parameters[c("tau.Wl"),c(site)]		#m3 cell wall day-1

    Carbon.daily.rate.ew<-daily.rate.ew*parameters[c("cell.wall.density.ew"),c(site)]	# kg C day-1	carbon to early wood wall formation
    Carbon.daily.rate.lw<-daily.rate.lw*parameters[c("cell.wall.density.lw"),c(site)]	# kg C day-1	carbon to late wood wall formation
  }

  if (site == "Hyde") {
    stem.no <- cbind(1997 : 2020, rep(1010, length.out = length(1997 : 2020)))	# Photosynthesis calculated with SPP-model for tree class 15-20 cm. Then compared with eddy GPP, determined with which stem no. the portion of eddy GPP is same as SPP estimate.
  } else if (site == "Lettosuo") {
    stem.no <- cbind(1997 : 2020, rep(500, length.out = length(1997 : 2020)))	# Photosynthesis calculated with PRELES or SPP-model.
  }

  # Carbon to height growth
  CH<-parameters[c("density_tree"),c(site)]*parameters[c("carbon_share"),c(site)]
  M.suc<-12*common[[c("M.C")]]+22*common[c("M.H")]+11*common[[c("M.O")]]


  #####
  ## Year loop
  #####

  LAI <- needle_mass <- NULL
  count <- 1

  out <- CASSIA_yearly()

  return(out)

}
