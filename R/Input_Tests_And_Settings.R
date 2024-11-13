validate_site <- function(site) {
  valid_sites <- c("Hyde", "Lettosuo", "Flakaliden_c", "HF_China")
  if (!(site %in% valid_sites)) {
    stop(paste("Unknown site. Please choose from:", paste(valid_sites, collapse = ", ")))
  }
}

update_model_settings <- function(settings) {
  if (settings$myco_model) {
    if (!settings$sperling_model) {
      settings$sperling_model <- TRUE
      warning("sperling_model set to TRUE as it should control allocation in myco_model.")
    }
    if (!settings$PRELES_GPP) {
      settings$PRELES_GPP <- TRUE
      warning("PRELES_GPP set to TRUE as PRELES should control the nitrogen effect in myco_model.")
    }
  }

  if (settings$sperling_model && settings$mycorrhiza) {
    settings$mycorrhiza <- FALSE
    warning("mycorrhiza set to FALSE as it's included explicitly in the Sperling submodel.")
  }

  if (!settings$sperling_model && settings$phloem_trigger) {
    settings$phloem.trigger <- FALSE
    warning("phloem.trigger set to FALSE as this feature requires the Sperling model.")
  }

  if (settings$xylogenesis) {
    settings$LN_estim <- FALSE
    settings$trees_grow <- TRUE
    settings$mycorrhiza <- FALSE
    settings$phloem_trigger <- FALSE
    warning("xylogenesis is TRUE: LN.estim, trees_grow, mycorrhiza, and phloem.trigger adjusted accordingly.")
  }

  return(settings)
}

validate_weather_data <- function(weather, PRELES_GPP, ecoevolutionary) {
  # Function to generate warning messages
  check_values <- function(column, min_val, max_val) {
    if (any(column < min_val | column > max_val, na.rm = TRUE)) {
      warning(paste("Values out of range for", deparse(substitute(column)), "Check input."))
    }
  }

  # Check temperature values
  check_values(weather$T, -30, 40)

  # Check PF values
  check_values(weather$PF, -10, 20)

  # Check Tsa values
  check_values(weather$Tsa, -10, 20)

  # Check Tsb values
  check_values(weather$Tsb, -10, 20)

  # Check M.soil values
  check_values(weather$M.soil, 0, 1)

  # Check Rain values
  check_values(weather$Rain, 0, 50)

  # Function to check for NA values and issue warnings
  check_na <- function(column) {
    if (any(is.na(column))) {
      stop(paste(deparse(substitute(column)), "has NA values. Check input."))
    }
  }

  # Check if weather data has the required columns
  if (PRELES_GPP) {
    required_columns <- c("date", "Date", "T", "P", "TSA", "TSB", "MB", "Rain", "PAR", "VPD", "fAPAR")
    if (!all(required_columns %in% names(weather))) {
      stop("Incomplete weather data - incorrect variables, or named incorrectly")
    }

    check_na(weather$T)
    check_na(weather$P)
    check_na(weather$Tsa)
    check_na(weather$Tsb)
    check_na(weather$M.soil) # MB?
    check_na(weather$Rain)
    check_na(weather$PAR)
    check_na(weather$VPD)
    check_na(weather$fAPAR)
  } else if (ecoevolutionary) {
    ### TODO: required columns for the weather data
  } else {
    required_columns <- c("date", "T", "P", "TSA", "TSB", "MB", "Rain")
    if (!all(required_columns %in% names(weather))) {
      stop("Incomplete weather data - incorrect variables, or named incorrectly")
    }

    check_na(weather$T)
    check_na(weather$P)
    check_na(weather$Tsa)
    check_na(weather$Tsb)
    check_na(weather$M.soil) # MB?
    check_na(weather$Rain)
  }
}

# Function to set up export data frames based on model type
setup_export_data <- function(model_type, years, total_days) {
  if (model_type == "sperling") {
    yearly_columns <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                        "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts",
                        "sugar.needles", "sugar.phloem", "sugar.xylem.sh", "sugar.xylem.st", "sugar.roots",
                        "starch.needles", "starch.phloem", "starch.xylem.sh", "starch.xylem.st", "starch.roots")
    daily_columns <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                       "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage_term", "to.mycorrhiza", "mycorrhiza.tot",
                       "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot",
                       "sugar.needles", "sugar.phloem", "sugar.xylem.sh", "sugar.xylem.st", "sugar.roots",
                       "starch.needles", "starch.phloem", "starch.xylem.sh", "starch.xylem.st", "starch.roots")
  } else if (model_type == "xylogenesis") {
    yearly_columns <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                        "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "ring_width", "needle_mass", "sum.needle.cohorts",
                        "ew_width", "ring_density", "ew.cells_tot", "lw.cells_tot")
    daily_columns <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                       "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage", "to.mycorrhiza", "mycorrhiza.tot",
                       "P", "to_sugar", "to_starch", "daily.consumption", "ring_width", "GD.tot", "n.E.tot", "n.W.tot", "n.M.tot")
  } else {
    yearly_columns <- c("year", "starch", "sugar", "wall.tot", "height.tot", "needle.tot", "root.tot", "tot.Rm", "tot.Rg",
                        "tot.P", "cumsum.PF", "cum.Daily.H.tot", "cum.Daily.N.tot", "tot.mm", "needle_mass", "sum.needle.cohorts")
    daily_columns <- c("date", "year", "day", "bud.tot.growth", "wall.tot.growth", "needle.tot.growth", "root.tot.growth", "height.tot.growth",
                       "Rg.tot", "Rm.tot", "height.tot", "wall.tot", "storage", "sugar", "starch", "storage_term", "to.mycorrhiza", "mycorrhiza.tot",
                       "P", "to_sugar", "to_starch", "Daily.H.tot", "Daily.N.tot", "GD.tot")
  }

  export_yearly <- data.table(matrix(ncol = length(yearly_columns), nrow = length(years)))
  export_daily <- data.table(matrix(ncol = length(daily_columns), nrow = total_days))
  names(export_yearly) <- yearly_columns
  names(export_daily) <- daily_columns

  return(list(export_yearly = export_yearly, export_daily = export_daily))
}
