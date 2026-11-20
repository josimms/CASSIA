#' Step-wise CASSIA run
#'
#' Run CASSIA one day (or a few days) at a time instead of a whole period in one call, e.g. to
#' couple it day by day with another model. \code{cassia_init()} takes the same arguments as
#' \code{\link{CASSIA_cpp}} and returns a run handle; \code{cassia_step()} advances it;
#' \code{cassia_outputs()} returns the same data frames as \code{CASSIA_cpp()} for the days simulated
#' so far. Stepping a run to its end gives exactly the result of \code{CASSIA_cpp()} with the same
#' arguments. Each run has its own state (also PRELES's soil water, snow and canopy water), so several
#' runs can be stepped side by side in one R session.
#'
#' @param weather,site,... As for \code{\link{CASSIA_cpp}} (standard path only:
#'   \code{soil = FALSE}, \code{ecoevolutionary = FALSE}).
#' @param run A run handle from \code{cassia_init()}.
#' @param n Number of days to advance (stops at the end of the weather).
#' @param values Named list or vector of drivers for the next day to simulate, using the weather
#'   column names of \code{CASSIA_cpp} (e.g. \code{TSA}, \code{TSB}, \code{MB}, \code{T},
#'   \code{Rain}, \code{P}, \code{PAR}, \code{VPD}, \code{CO2}, \code{fAPAR}); they replace that day's
#'   values in the weather given to \code{cassia_init()}.
#' @return \code{cassia_init()}: a run handle (class \code{cassia_run}). \code{cassia_step()}: the
#'   number of days simulated so far (invisibly). \code{cassia_day()}: list with \code{date} of the
#'   next day to simulate (\code{NA} when finished), \code{days_done}, \code{days_total},
#'   \code{done}. \code{cassia_outputs()}: list of data frames as \code{CASSIA_cpp()}, one row per day
#'   simulated so far.
#' @examples
#' \dontrun{
#' run <- cassia_init(weather = weather_original, site = "Hyde")
#' cassia_step(run, 30)
#' cassia_set_forcing(run, list(TSA = 4.2))
#' cassia_step(run)
#' head(cassia_outputs(run)$Growth)
#' }
#' @name cassia_step_api
NULL

#' @rdname cassia_step_api
#' @export
cassia_init <- function(weather, site, ...) {
  CASSIA_cpp(weather = weather, site = site, ..., .return_run = TRUE)
}

#' @rdname cassia_step_api
#' @export
cassia_step <- function(run, n = 1) {
  stopifnot(inherits(run, "cassia_run"), n >= 0)
  invisible(cassia_step_cpp(run$ptr, as.integer(n)))
}

#' @rdname cassia_step_api
#' @export
cassia_day <- function(run) {
  stopifnot(inherits(run, "cassia_run"))
  cursor <- cassia_cursor_cpp(run$ptr)
  list(date = if (cursor$done) as.Date(NA) else run$dates[cursor$days_done + 1],
       days_done = cursor$days_done, days_total = cursor$days_total, done = cursor$done)
}

#' @rdname cassia_step_api
#' @export
cassia_set_forcing <- function(run, values) {
  stopifnot(inherits(run, "cassia_run"), !is.null(names(values)))
  cassia_set_forcing_cpp(run$ptr, as.list(values))
  invisible(run)
}

#' @rdname cassia_step_api
#' @export
cassia_outputs <- function(run) {
  stopifnot(inherits(run, "cassia_run"))
  n <- cassia_cursor_cpp(run$ptr)$days_done
  lapply(cassia_outputs_cpp(run$ptr), function(df) {
    if (n >= nrow(df)) return(df)
    df <- df[seq_len(n), , drop = FALSE]
    rownames(df) <- NULL
    df
  })
}
