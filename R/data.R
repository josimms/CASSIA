#' Daily weather data for Hyytiälä SMEAR II, 2015–2017
#'
#' Three years of daily weather observations from the Hyytiälä SMEAR II
#' measurement station in southern Finland. Used as the default weather input
#' for running CASSIA with built-in Hyytiälä parameters.
#'
#' Before passing to \code{\link{CASSIA_cpp}} you must add a \code{dates}
#' column:
#' \preformatted{
#' weather_original$dates <- as.Date(
#'   strptime(paste(rep(2015:2017, times = c(365, 366, 365)),
#'                  weather_original$X),
#'            format = "\%Y \%j"))
#' }
#'
#' @format A data frame with 1096 rows (365 + 366 + 365 days) and 12 columns:
#' \describe{
#'   \item{X}{Day of year (1–365 or 1–366).}
#'   \item{T}{Daily mean air temperature (°C).}
#'   \item{P}{Daily gross primary production used as photosynthesis input when
#'     \code{photosynthesis_as_input = TRUE} (g C m\eqn{^{-2}} d\eqn{^{-1}}).}
#'   \item{TSA}{Daily mean soil temperature, A horizon (°C).}
#'   \item{TSB}{Daily mean soil temperature, B horizon (°C).}
#'   \item{MB}{Daily mean volumetric soil moisture, B horizon (m\eqn{^3}
#'     m\eqn{^{-3}}).}
#'   \item{Rain}{Daily precipitation (mm d\eqn{^{-1}}).}
#'   \item{Nitrogen}{Atmospheric nitrogen input (g N m\eqn{^{-2}} d\eqn{^{-1}}).}
#'   \item{PAR}{Daily mean photosynthetically active radiation
#'     (\eqn{\mu}mol m\eqn{^{-2}} s\eqn{^{-1}}).}
#'   \item{VPD}{Daily mean vapour pressure deficit (kPa).}
#'   \item{CO2}{Atmospheric CO\eqn{_2} concentration (ppm).}
#'   \item{fAPAR}{Fraction of absorbed PAR (dimensionless, 0–1).}
#' }
#' @source SMEAR II station, University of Helsinki. Processed from raw
#'   half-hourly measurements by the scripts in \code{data-raw/Hyytiala.R}.
#'   See \code{vignette("Weather")} for details on the processing workflow.
"weather_original"


#' Main CASSIA model parameters
#'
#' A matrix of site-specific parameters used by \code{\link{CASSIA_cpp}} to
#' control carbon allocation, phenology, and growth of each organ. One column
#' per site; rows are named parameters.
#'
#' To modify a parameter, copy the matrix and change the value:
#' \preformatted{
#' p <- parameters_p
#' p["root.lifetime", "Hyde"] <- 2
#' out <- CASSIA_cpp(weather = weather_original, site = "Hyde",
#'                   parameters = p)
#' }
#'
#' @format A numeric matrix with 75 rows (parameters) and 3 columns (sites:
#'   \code{Hyde}, \code{Lettosuo}, \code{Väriö}). Row names are parameter
#'   names. All parameter meanings are described in the CASSIA instruction
#'   booklet (\file{CASSIA_Instruction_Booklet.pdf}) and in
#'   Schiestl-Aalto et al. (2015).
#' @source Built by \code{data-raw/building_data.R}.
#' @seealso \code{\link{ratios_p}}, \code{\link{common_p}},
#'   \code{\link{sperling_p}}, \code{\link{repo_p}}
"parameters_p"


#' Allocation ratio parameters
#'
#' Site-specific allocation ratio parameters used by \code{\link{CASSIA_cpp}}.
#' These control ratios between organ growth rates (e.g. needle to fine-root
#' allocation, sapwood share, growth coefficients for height and diameter).
#'
#' @format A numeric matrix with 9 rows (parameters) and 2 columns (sites:
#'   \code{Hyde}, \code{Lettosuo}). Row names are:
#' \describe{
#'   \item{form_factor}{Tree form factor for volume calculations.}
#'   \item{needle_fineroot_ratio}{Ratio of needle mass to fine-root mass.}
#'   \item{sapwood.share}{Fraction of stemwood that is sapwood.}
#'   \item{height_growth_coefficient}{Base coefficient for height growth.}
#'   \item{diameter_growth_coefficient}{Base coefficient for diameter growth.}
#'   \item{height_growth_coefficient_max}{Upper bound on height growth coefficient.}
#'   \item{height_growth_coefficient_min}{Lower bound on height growth coefficient.}
#'   \item{diameter_growth_coefficient_max}{Upper bound on diameter growth coefficient.}
#'   \item{diameter_growth_coefficient_min}{Lower bound on diameter growth coefficient.}
#' }
#' @source Built by \code{data-raw/building_data.R}.
#' @seealso \code{\link{parameters_p}}
"ratios_p"


#' Common (site-independent) CASSIA parameters
#'
#' Parameters that are shared across all sites and passed to
#' \code{\link{CASSIA_cpp}} via the \code{common_parameters} argument. These
#' include soil physical constants, gas constants, and model-wide constants.
#'
#' @format A numeric matrix with 1 row and 23 columns. Column names are the
#'   parameter names:
#'   \code{a}, \code{b}, \code{TR0}, \code{abs_zero}, \code{b.s},
#'   \code{theetta.FC}, \code{phi.e}, \code{K.sat}, \code{R.length},
#'   \code{M.H20}, \code{r.cyl}, \code{r.root}, \code{ypsilon},
#'   \code{Rg.N}, \code{Rg.S}, \code{Rg.R}, \code{gas.const},
#'   \code{M.C}, \code{M.H}, \code{M.O}, \code{osmotic.sugar.conc},
#'   \code{m_N}, \code{Uggla}.
#' @source Built by \code{data-raw/building_data.R}.
#' @seealso \code{\link{parameters_p}}
"common_p"


#' Sugar and xylogenesis model parameters
#'
#' Parameters for the organ-level sugar sub-model and the xylogenesis
#' sub-model. Passed to \code{\link{CASSIA_cpp}} via the
#' \code{sperling_parameters} argument.
#'
#' @format A numeric matrix with 55 rows (parameters) and 3 columns (sites:
#'   \code{Hyde}, \code{Lettosuo}, \code{Flakaliden_c}). Row names are
#'   parameter names. For a step-by-step description of how the xylogenesis
#'   parameters are used, see \code{vignette("ring_width_process")}.
#' @source Built by \code{data-raw/building_data.R}.
#' @seealso \code{\link{parameters_p}}, \code{\link{CASSIA_cpp}}
"sperling_p"


#' Repola allometry parameters
#'
#' Parameters for the Repola (2009) allometric equations that convert CASSIA
#' growth outputs to biomass components. Passed internally by
#' \code{\link{repola}}.
#'
#' @format A numeric matrix with 1 row and 5 columns:
#' \describe{
#'   \item{b0.repo}{Intercept for biomass equation.}
#'   \item{b1.repo}{Diameter-at-breast-height coefficient.}
#'   \item{b2.repo}{Height coefficient.}
#'   \item{uk.repo}{Uncertainty scaling factor.}
#'   \item{eki.repo}{Exponent correction term.}
#' }
#' @source Repola, J. (2009). Biomass equations for Scots pine and Norway
#'   spruce in Finland. *Silva Fennica* 43(4): 625–647.
#' @seealso \code{\link{repola}}
"repo_p"


#' Reference daily GPP for Hyytiälä
#'
#' A vector of 365 daily gross primary production values representing a
#' climatological reference year at Hyytiälä. Used internally by
#' \code{\link{CASSIA_cpp}} to scale inter-annual GPP variation when
#' \code{LN_estim = TRUE}.
#'
#' @format A named numeric vector of length 365. Each element is the mean GPP
#'   for that day of year (g C m\eqn{^{-2}} d\eqn{^{-1}}).
#' @source Derived from multi-year Hyytiälä SMEAR II eddy-covariance
#'   measurements. Built by \code{data-raw/Hyytiala.R}.
"GPP_ref"


#' Soil and nitrogen sub-model parameters
#'
#' A numeric vector of parameters used by the soil and mycorrhizal sub-models
#' in \code{\link{CASSIA_cpp}}. These control nitrogen uptake kinetics,
#' microbial dynamics, and soil carbon pools.
#'
#' @format A numeric vector of length 58. Individual elements are indexed
#'   positionally; parameter meanings are described in the CASSIA instruction
#'   booklet (\file{man/CASSIA_Instruction_Booklet.pdf}) and in
#'   \code{data-raw/Mycorrhizal_Model_Parameters.R}.
#' @source Built by \code{data-raw/Mycorrhizal_Model_Parameters.R}.
#' @seealso \code{\link{parameters_p}}
"parameters_R"
