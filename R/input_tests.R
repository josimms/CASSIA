function(weather) {

  #####
  ## Input tests!
  #####
  # Check that the sites are within the sites allowed
  if ((site %in% c("Hyde", "Lettosuo", "Flakaliden_c")) == F) {stop("Unknown site: Please pick between Hyde and Lettosuo")}

  ### TODO: Change these tests so they are for the dataframe
  if (sum(Temp < -30) + sum(Temp > 30) != 0) {warning(paste("Temp < -30 or Tempp > 30 in year", year,  ". Values not impossible, but unlikely check input."))}
  if (sum(PF < -10) + sum(PF > 20) != 0) {warning(paste("PF < -10 or PF > 20 in year", year, ". Values not impossible, but unlikely check input."))}
  if (sum(Tsa < -10) + sum(Tsa > 20) != 0) {warning(paste("Tsa < -10 or Tsa > 20 in year", year, ". Values not impossible, but unlikely check input."))}
  if (sum(Tsb < -10) + sum(Tsb > 20) != 0) {warning(paste("Tsb < -10 or Tsb > 20 in year", year, ". Values not impossible, but unlikely check input."))}
  if (sum(M.soil < 0) + sum(M.soil > 1) != 0) {warning(paste("M.soil < 0 or M.soil > 1 in year", year, ". Values not impossible, but unlikely check input."))}
  if (sum(Rain < 0) + sum(Rain > 50) != 0) {warning(paste("Rain < 0 or Rain > 50 in year", year, ". Values not impossible, but unlikely check input."))}

  if (sum(is.na(Temp))) {warning(paste("Temp has NA values check input."))}
  if (sum(is.na(PF))) {warning(paste("PF has NA values check input."))}
  if (sum(is.na(Tsa))) {warning(paste("Tsa has NA values check input."))}
  if (sum(is.na(Tsb))) {warning(paste("Tasb has NA values check input."))}
  if (sum(is.na(M.soil))) {warning(paste("M.soil has NA values check input."))}
  if (sum(is.na(Rain))) {warning(paste("Rain has NA values check input."))}

}
