#' Run the CASSIA intra-annual tree growth model
#'
#' Simulates daily carbon allocation, organ growth, and optionally sugar
#' dynamics and mycorrhizal interactions for an individual Scots pine in boreal
#' conditions. Built-in parameter sets and weather data are included for
#' Hyytiälä (\code{"Hyde"}) and Lettosuo.
#'
#' @param weather Daily dataframe. Which columns are required depends on the
#'   photosynthesis sub-model selected; columns not needed by the active
#'   sub-model may be \code{NA} or absent. Always required:
#'   \describe{
#'     \item{dates}{Date (class \code{Date}, YYYY-MM-DD)}
#'     \item{T}{Air temperature (°C)}
#'     \item{TSA}{Soil temperature, A horizon (°C)}
#'     \item{TSB}{Soil temperature, B horizon (°C)}
#'     \item{MB}{Soil moisture, B horizon (m3 m-3)}
#'     \item{Rain}{Precipitation (mm day-1)}
#'   }
#'   Additional columns required by each sub-model:
#'   \describe{
#'     \item{P}{GPP per tree (g C m-2 day-1) —
#'       \code{photosynthesis_as_input = TRUE} (default)}
#'     \item{PAR, VPD, CO2, fAPAR}{Photosynthetically active radiation
#'       (mol m-2 day-1), vapour pressure deficit (kPa), CO2 concentration
#'       (ppm), fraction of absorbed PAR (0–1) —
#'       \code{preles = TRUE}}
#'     \item{PPFD, PPFD_max, VPD, CO2, Nitrogen, PA, SWP}{Photosynthetic
#'       photon flux density and daily maximum (µmol m-2 s-1), VPD (kPa),
#'       CO2 (ppm), leaf nitrogen (g m-2), atmospheric pressure (Pa),
#'       soil water potential (MPa) —
#'       \code{phydro = TRUE}}
#'   }
#' @param site Site name for built-in parameterisations: \code{"Hyde"}
#'   (Hyytiälä), \code{"Lettosuo"}, or \code{"HF_China"}. Add new sites by
#'   extending the \code{parameters_p} dataframe.
#' @param GPP_ref_in Reference GPP table used when \code{LN_estim = TRUE}.
#'   Defaults to the built-in \code{GPP_ref}.
#' @param ratio_sugar Numeric vector (length 4). Concentration difference ratios
#'   between organs used in the sugar transfer model. Default
#'   \code{c(11, 1/3, 11, 11)}.
#' @param tau.myco Numeric. Time constant for sugar transfer to mycorrhiza (days).
#'   Default \code{3}.
#' @param tau.t.needles Numeric. Time constant for sugar transfer to needles (days).
#'   Default \code{3}.
#' @param tau.t.phloem Numeric. Time constant for sugar transfer via phloem (days).
#'   Default \code{3}.
#' @param tau.t.roots Numeric. Time constant for sugar transfer to roots (days).
#'   Default \code{3}.
#' @param tau.t.xylem.sh Numeric. Time constant for sugar transfer to xylem
#'   shoot (days). Default \code{3}.
#' @param tau.t.xylem.st Numeric. Time constant for sugar transfer to xylem
#'   stem (days). Default \code{3}.
#' @param ratios Dataframe of growth ratio parameters. Defaults to \code{ratios_p}.
#' @param parameters Dataframe of site-specific parameters. Defaults to
#'   \code{parameters_p}.
#' @param parameters_R_in Dataframe of R-side parameters. Defaults to
#'   \code{parameters_R}.
#' @param common Dataframe of parameters shared across sites. Defaults to
#'   \code{common_p}.
#' @param sperling Dataframe of sugar model parameters. Defaults to
#'   \code{sperling_p}.
#' @param repo Dataframe of Repola allometric parameters. Defaults to
#'   \code{repo_p}.
#' @param pPREL Numeric vector of PRELES photosynthesis model parameters
#'   (length 32). Only used when \code{preles = TRUE}. Default is the
#'   published Hyytiälä calibration.
#' @param phydro_param Numeric vector. Parameters for the p-hydro photosynthesis
#'   model. Only used when \code{phydro = TRUE}.
#' @param storage_reset Logical. If \code{TRUE} (default), sugar and starch
#'   stores reset to initial values at the start of each year. If \code{FALSE},
#'   carry over from the previous simulated year.
#' @param storage_grows Logical. If \code{TRUE}, potential storage space grows
#'   with tree height and diameter. Default \code{FALSE}.
#' @param LN_estim Logical. If \code{TRUE} (default), needle elongation (LN) is
#'   scaled by the previous year's July-August GPP. If \code{FALSE}, a fixed
#'   coefficient is used.
#' @param mN_varies Logical. If \code{TRUE} (default), needle mass in
#'   maintenance respiration uses the Repola model. If \code{FALSE}, constant.
#' @param LD_estim Logical. If \code{TRUE} (default), LD is modified by GPP
#'   over the preceding March-August period.
#' @param sD_estim_T_count Logical. If \code{TRUE}, sD is the cumulative count
#'   of days with positive growth. If \code{FALSE} (default), sD is derived
#'   from the temperature-driven growth integral.
#' @param LH_estim Logical. If \code{TRUE} (default), LH is modified by GPP.
#' @param soil_moisture_effect_on_shoot Logical. Enable soil moisture limitation
#'   on shoot growth. Default \code{FALSE}.
#' @param soil_moisture_effect_on_needles Logical. Enable soil moisture
#'   limitation on needle growth. Default \code{FALSE}.
#' @param soil_moisture_effect_on_diameter Logical. Enable soil moisture
#'   limitation on diameter growth. Default \code{FALSE}.
#' @param driver_N Numeric vector of length 2. Scaling parameters controlling
#'   the effect of nitrogen availability on growth. Default \code{c(1, 1)}.
#' @param driver_H Numeric vector of length 2. Scaling parameters controlling
#'   height growth. Default \code{c(0, 2)}.
#' @param driver_D Numeric vector of length 2. Scaling parameters controlling
#'   diameter growth. Default \code{c(1, 1)}.
#' @param trees_grow Logical. If \code{TRUE}, D0 and h0 are updated
#'   year-to-year from modelled growth. Automatically \code{TRUE} when
#'   \code{xylogenesis = TRUE}. Default \code{FALSE}.
#' @param growth_decreases Logical. If \code{TRUE}, height and diameter growth
#'   coefficients decline linearly from 1997 to 2020. Default \code{FALSE}.
#' @param needle_mass_grows Logical. If \code{TRUE}, needle mass is updated
#'   dynamically each year via the Repola model. If \code{FALSE} (default),
#'   needle mass is fixed after the first iteration.
#' @param phloem_trigger Logical. If \code{TRUE}, bud burst is triggered when
#'   phloem sugar first drops below a threshold (requires
#'   \code{organ_level_sugar = TRUE}). Default \code{FALSE}.
#' @param mycorrhiza Logical. If \code{TRUE} (default), carbon allocation to
#'   mycorrhiza is included. Automatically \code{FALSE} when
#'   \code{organ_level_sugar = TRUE}.
#' @param root_as_Ding Logical. If \code{TRUE} (default), root growth follows
#'   Ding et al. (2020) with fibrous/non-fibrous root stages.
#' @param organ_level_sugar Logical. If \code{TRUE}, uses the organ-level sugar
#'   model with per-organ concentration pools (requires \code{sperling_p}
#'   parameters). If \code{FALSE} (default), uses the Schiestl-Aalto (2019)
#'   equilibrium approach with a single aggregate pool.
#' @param myco_model Logical. Enable experimental mycorrhizal model.
#'   Default \code{FALSE}.
#' @param xylogenesis Logical. If \code{TRUE}, xylogenesis controls cell growth
#'   with early/late-wood stages. Default \code{FALSE}.
#' @param environment_effect_xylogenesis Logical. If \code{TRUE}, temperature
#'   affects cell formation during xylogenesis. Default \code{FALSE}.
#' @param photosynthesis_as_input Logical. If \code{TRUE} (default), the
#'   \code{P} column in \code{weather} is used directly as GPP input.
#' @param photoparameters Integer (1-4). PRELES parameter set selector.
#'   Default \code{3}.
#' @param temp_rise Logical. If \code{TRUE} and \code{Rm_acclimation = TRUE},
#'   acclimation factor Rm_accl is set to 0.85. Default \code{FALSE}.
#' @param drought Logical. Currently unused. Default \code{FALSE}.
#' @param Rm_acclimation Logical. Apply maintenance respiration acclimation
#'   when \code{temp_rise = TRUE}. Default \code{TRUE}.
#' @param CASSIA_graphs Logical. Produce diagnostic plots during the run.
#'   Default \code{TRUE}.
#' @param etmodel Logical. If \code{TRUE}, use the Penman-Monteith
#'   evapotranspiration model inside PRELES. Only relevant when
#'   \code{preles = TRUE}. Default \code{FALSE}.
#' @param LOGFLAG Logical. If \code{TRUE}, write internal diagnostic output to
#'   the console. Intended for debugging. Default \code{FALSE}.
#' @param s.D0 Integer. Day-of-year to start the temperature sum for diameter
#'   growth (default 79, approx. March 20; valid for Finland).
#' @param s.H0 Integer. Day-of-year to start shoot growth. Default \code{1}.
#' @param growth_photo_coef Numeric. Direct photosynthesis growth coefficient
#'   override. Default \code{1} (no override).
#' @param needle_mass_in Numeric. Initial needle mass (kg C). Default
#'   \code{4.467638} (calibrated for Hyytiälä).
#' @param Throughfall Numeric. Throughfall fraction. Default \code{1}.
#' @param trenching_year Integer or \code{NA}. Year to apply trenching
#'   (mycorrhiza cut-off). Default \code{NA}.
#' @param nitrogen_balance Numeric. Baseline nitrogen balance value.
#'   Default \code{25}.
#' @param nitrogen_change Logical. Enable nitrogen addition/change treatment.
#'   Default \code{FALSE}.
#' @param nitrogen_contrast Logical. Enable nitrogen contrast treatment.
#'   Default \code{FALSE}.
#' @param surplus_c Logical. Enable surplus carbon scenario. Default
#'   \code{FALSE}.
#' @param soil Logical. If \code{TRUE}, runs the soil-coupled variant
#'   \code{CASSIA_soil}. Default \code{FALSE}.
#' @param ecoevolutionary Logical. If \code{TRUE}, runs the eco-evolutionary
#'   variant \code{CASSIA_eeo}. Default \code{FALSE}.
#' @param fAPAR_Tian Logical. Apply the fAPAR correction from Tian et al.
#'   (2021). Default \code{FALSE}.
#' @param phydro Logical. If \code{TRUE}, GPP is calculated via the p-hydro
#'   model. Default \code{FALSE}.
#' @param preles Logical. If \code{TRUE}, GPP is calculated via PRELES.
#'   Default \code{FALSE}.
#' @param tests Logical. Enable internal diagnostic tests. Default \code{FALSE}.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{Growth}{Dataframe of daily carbon fluxes and organ growth (bud,
#'       wall, needle, root, height growth; respiration; storage; mycorrhiza).}
#'     \item{Sugar}{Dataframe of daily sugar and starch pool values per organ.
#'       Populated when the sugar model is active.}
#'     \item{Preles}{Dataframe of daily photosynthesis outputs: GPP, ET,
#'       SoilWater.}
#'   }
#'
#' @references
#'   Schiestl-Aalto et al. (2015) CASSIA - a dynamic model for predicting
#'   intra-annual sink demand. \emph{New Phytologist} 206(2):647-659.
#'
#'   Schiestl-Aalto et al. (2019) Analysis of NSC storage dynamics in tree
#'   organs. \emph{Frontiers in Forests and Global Change} 2:17.
#'
#'   Ding et al. (2020) Temperature and moisture dependence of daily growth of
#'   Scots pine roots. \emph{Tree Physiology} 40(2):272-283.
#'
#'   Sperling et al. (2019) Predicting bloom dates by temperature mediated
#'   kinetics of carbohydrate metabolism. \emph{Agricultural and Forest
#'   Meteorology} 276:107643.
#'
#'   Tian et al. (2021) Disaggregating the effects of nitrogen addition on
#'   gross primary production. \emph{Agricultural and Forest Meteorology}
#'   301:108337.
#'
#' @examples
#' \dontrun{
#' weather_original$dates <- as.Date(
#'   strptime(paste(rep(2015:2017, times = c(365, 366, 365)), weather_original$X),
#'            format = "%Y %j"))
#' out <- CASSIA_cpp(weather = weather_original, site = "Hyde")
#' head(out$Growth)
#'
#' # Modify a parameter
#' p_new <- parameters_p
#' p_new["root.lifetime", "Hyde"] <- 2
#' out2 <- CASSIA_cpp(weather_original, site = "Hyde", parameters = p_new)
#' }
#' @export
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
    GPP_ref_in = GPP_ref,

    ratio_sugar = c(11, 1/3, 11, 11), # Determines the concentration difference between organs
    tau.myco = 3,
    tau.t.needles = 3,
    tau.t.phloem = 3,
    tau.t.roots = 3,
    tau.t.xylem.sh = 3,
    tau.t.xylem.st = 3,

    ratios = ratios_p,
    parameters = parameters_p,
    parameters_R_in = parameters_R,
    common = common_p,
    sperling = sperling_p,
    repo = repo_p,
    pPREL = c(413.0, 0.450, 0.118, 3.0,
              0.745700, 10.930000, -3.06300, 17.720000, -0.102700, 0.036730, 0.777900, 0.500, -0.364,
              0.271500, 0.835100, 0.073480, 0.999600, 0.442800,
              1.2, 0.33, 4.970496, 0.0, 0.0,
              160.0, 0.0, 0.0, 20.0,
              -999.9, -999.9, -999.9,
              1/0.012, 0.0),
    phydro_param = c(0.1008, 0.180496537959982, 5, 0.026263945805926, 0.011, 50,
               0.5, -0.857817410110663, 4.1311874912949e17, 1, 2.45e-2, 2.0, 1.1, 0.1,
               15, 10, 5, 0, 1, exp(-0.5 * 1.8), exp(-0.5 * 3.5), exp(-0.5 * 5.5)),

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

    soil_moisture_effect_on_shoot = FALSE,
    soil_moisture_effect_on_needles = FALSE,
    soil_moisture_effect_on_diameter = FALSE,
    driver_N = c(1, 1),
    driver_H = c(0, 2),
    driver_D = c(1, 1),

    trees_grow = FALSE,				# can be false if mature trees are modelled and not for a very long period
    growth_decreases = FALSE,			# the height and diameter growth (alfa_S and alfaD) decrease during the simulation
    needle_mass_grows = FALSE,		# Is needle mass dynamic i.e. the modelled growth is also respiring etc and following for some years? If true, note that root mass is related to needle mass

    phloem_trigger = FALSE,    # Phloem controls bud burst rather than whole tree sugar

    mycorrhiza = TRUE, 			# If allocation to mychorrhiza is taken into account
    root_as_Ding = TRUE,

    organ_level_sugar = FALSE,    # Organ-level sugar model with per-organ concentration pools
    myco_model = FALSE,           # Joanna's mycomodel development!
    xylogenesis = FALSE,

    environment_effect_xylogenesis = FALSE,

    photosynthesis_as_input = TRUE,

    photoparameters = 3,
    temp_rise = FALSE,
    drought = FALSE,
    Rm_acclimation = TRUE,

    CASSIA_graphs = TRUE,

    etmodel = F,
    LOGFLAG = F,

    s.D0 = 79,					# DOY to start the calculation of temperature sum, 1=Jan 1; 69=March 1; 79=March 20 for diameter growth. Valid for Finland
    s.H0 = 1,					# and for shoot grwoth

    growth_photo_coef = 1,
    needle_mass_in = 4.467638,
    Throughfall = 1,
    trenching_year = NA,

    nitrogen_balance = 25,
    nitrogen_change = FALSE,
    nitrogen_contrast = FALSE,
    surplus_c = FALSE,

    soil = FALSE,
    ecoevolutionary = FALSE,
    fAPAR_Tian = FALSE,

    phydro = FALSE,
    preles = FALSE,

    tests = FALSE) {

  ####
  # Testing the input and initial conditions
  ####

  # Is the site is the known sites?
  validate_site(site)

  no_trees = 1010
  if (site == "HF_China") {no_trees = 1044}

  # make the model settings into a list
  settings = list(
    "storage_reset" = storage_reset,
    "storage_grows" = storage_grows,
    "LN_estim" = LN_estim,
    "mN_varies" = mN_varies,
    "LD_estim" = LD_estim,
    "sD_estim_T_count" = sD_estim_T_count,
    "LH_estim" = LH_estim,
    "trees_grow" = trees_grow,
    "growth_decreases" = growth_decreases,
    "needle_mass_grows" = needle_mass_grows,
    "phloem_trigger" = phloem_trigger,
    "mycorrhiza" = mycorrhiza,
    "root_as_Ding" = root_as_Ding,
    "organ_level_sugar" = organ_level_sugar,
    "myco_model" = myco_model,
    "xylogenesis" = xylogenesis,
    "environment_effect_xylogenesis" = environment_effect_xylogenesis,
    "photosynthesis_as_input" = photosynthesis_as_input,
    "photoparameters" = photoparameters,
    "temp_rise" = temp_rise,
    "drought" = drought,
    "Rm_acclimation" = Rm_acclimation,
    "CASSIA_graphs" = CASSIA_graphs,
    "tests" = tests,
    "etmodel" = etmodel,
    "LOGFLAG" = LOGFLAG,
    "ecoevolutionary" = ecoevolutionary,
    "fAPAR_Tian" = fAPAR_Tian,
    "preles" = preles,
    "phydro" = phydro,
    "soil_moisture_effect_on_shoot" = soil_moisture_effect_on_shoot,
    "soil_moisture_effect_on_needles" = soil_moisture_effect_on_needles,
    "soil_moisture_effect_on_diameter" = soil_moisture_effect_on_diameter,
    "driver_N" = driver_N,
    "driver_H" = driver_H,
    "driver_D" = driver_D
  )

  # Are the model settings valid?
  # If not update them
  updated_settings <- update_model_settings(settings)

  # Is the weather data correct?
  if (sum(c(photosynthesis_as_input, preles, phydro)) > 1) {
    stop("Only one of the photosynthesis models can be use at a time.")
  }

  if (sum(c(photosynthesis_as_input, preles, phydro)) < 1) {
    stop("No photosynthesis model is selected.")
  }

  # TODO: add the weather tests again

  #####
  ## Model conditions derived from model inputs
  #####
  # years from weather data
  date_range = as.numeric(substring(weather$dates[c(1, nrow(weather))], 1, 4))
  if (sum(date_range %in% 0:2500) < 2) {
    stop("Dates are not between 0 and 2500. Is the column called dates?")
  }
  start_year <- date_range[1]
  end_year <- date_range[2]

  if (soil) {
    out <- CASSIA_soil(start_year, end_year, weather, GPP_ref_in,
                       pPREL, t(parameters), common, t(ratios), t(sperling), parameters_R_in, # site,
                       needle_mass_in,
                       Throughfall, surplus_c, nitrogen_change, nitrogen_contrast, nitrogen_balance, trenching_year,
                       updated_settings)
  } else if (ecoevolutionary) {
    out <- CASSIA_eeo(start_year, end_year, weather, GPP_ref_in,
                      pPREL, t(parameters), common, t(ratios), t(sperling), parameters_R_in, phydro_param, # site,
                      needle_mass_in,
                      Throughfall, surplus_c, nitrogen_change, nitrogen_contrast, nitrogen_balance, trenching_year,
                      updated_settings)
  } else {
    out <- CASSIA_yearly(start_year, end_year, weather, GPP_ref_in,
                         pPREL, t(parameters), common, t(ratios), t(sperling), no_trees,
                         needle_mass_in,
                         Throughfall, surplus_c, nitrogen_change, nitrogen_contrast, nitrogen_balance,
                         updated_settings)
  }

  return(out)
}

