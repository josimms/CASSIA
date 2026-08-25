#' Low-level C++ model runners
#'
#' \code{CASSIA_yearly} and \code{CASSIA_soil} are the underlying compiled
#' model runners called internally by \code{\link{CASSIA_cpp}}. Under normal
#' use you should call \code{CASSIA_cpp} instead, which validates inputs,
#' selects the correct runner based on the \code{soil} and
#' \code{ecoevolutionary} toggles, and returns a consistent output format.
#'
#' \code{CASSIA_yearly} runs the standard (no-soil) annual model.
#' \code{CASSIA_soil} runs the soil-coupled variant (equivalent to calling
#' \code{CASSIA_cpp(..., soil = TRUE)}).
#'
#' @param start_year Integer. First simulation year.
#' @param end_year Integer. Last simulation year.
#' @param weather Dataframe of daily weather inputs (same format as
#'   \code{\link{CASSIA_cpp}}).
#' @param GPP_ref Numeric vector of 365 reference GPP values
#'   (see \code{\link{GPP_ref}}).
#' @param pPREL Numeric vector of PRELES parameters (length 32).
#' @param pCASSIA_parameters Transposed parameter matrix from
#'   \code{\link{parameters_p}}.
#' @param pCASSIA_common Common parameter matrix from \code{\link{common_p}}.
#' @param pCASSIA_ratios Transposed ratio matrix from \code{\link{ratios_p}}.
#' @param pCASSIA_sperling Transposed xylogenesis matrix from
#'   \code{\link{sperling_p}}.
#' @param no_trees Integer. Stand density (trees per hectare). Default 1010
#'   for Hyytiälä.
#' @param needle_mass_in Numeric. Initial needle mass (kg C).
#' @param Throughfall Numeric. Throughfall fraction.
#' @param settings Named list of model toggle flags. Built by
#'   \code{CASSIA_cpp} from its individual logical arguments.
#' @param parameters_R Numeric vector of soil sub-model parameters (only in
#'   \code{CASSIA_soil}; see \code{\link{parameters_R}}).
#' @param trenching_year Integer or \code{NA}. Year to apply trenching (only
#'   in \code{CASSIA_soil}).
#'
#' @return A named list with elements \code{Growth}, \code{Sugar}, and
#'   \code{Preles} (same structure as \code{\link{CASSIA_cpp}}).
#'
#' @seealso \code{\link{CASSIA_cpp}} for the recommended high-level interface.
#' @name CASSIA_yearly
NULL
