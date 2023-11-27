PRELES_plot <- function(data_format, N_parameters) {

  # install.packages("~/Documents/Rprebasso-master/",repos=NULL, type="source")

  # devtools::install_github("ForModLabUHel/Rprebasso")
  # library(Rprebasso)

  original <- PRELES(TAir = weather_original$T,
                     Precip = weather_original$Rain,
                     CO2 = weather_original$CO2,
                     PAR = weather_original$PAR,
                     VPD = weather_original$VPD,
                     Nitrogen = weather_original$Nitrogen,
                     fAPAR = rep(0.7, length = nrow(weather_original)),
                     p = c(pPREL, N_parameters))

  plot(original$GPP)

  data_format_nit_const = cbind(data_format, Nitrogen = rep(0.012, nrow(data_format)))

  # TODO: check the inputs!
  GPP <- ET <- SWC <- S <- NULL
  for (day in 1:nrow(data_format_nit_const)) {
    GPP[day] <- preles_test_cpp(nrow(data_format_nit_const), day, data_format_nit_const[day,], pPREL, N_parameters[1], N_parameters[2])$GPP
    ET[day] <- preles_test_cpp(nrow(data_format_nit_const), day, data_format_nit_const[day,], pPREL, N_parameters[1], N_parameters[2])$ET
    SWC[day] <- preles_test_cpp(nrow(data_format_nit_const), day, data_format_nit_const[day,], pPREL, N_parameters[1], N_parameters[2])$SWC
    S[day] <- preles_test_cpp(nrow(data_format_nit_const), day, data_format_nit_const[day,], pPREL, N_parameters[1], N_parameters[2])$S
  }

  Nitrgen = 0:1500
  data_format_nit = data.frame(cbind(T = rep(15, length = length(Nitrgen)),
                                     Rain = rep(2.6, length = length(Nitrgen)),
                                     PAR = rep(14, length = length(Nitrgen)),
                                     CO2 = rep(366, length = length(Nitrgen)),
                                     VPD = rep(0.2571383, length = length(Nitrgen)),
                                     fAPAR = rep(0.7, length = length(Nitrgen)),
                                     Nitrogen = Nitrgen))

  count = 1
  GPP_N <- ET_N <- SWC_N <- S_N <- NULL
  for (nit in 1:1501) {
    GPP_N[count] <- preles_test_cpp(nrow(data_format_nit), count, data_format_nit[count,], pPREL, N_parameters[1], N_parameters[2])$GPP
    ET_N[count] <- preles_test_cpp(nrow(data_format_nit), count, data_format_nit[count,], pPREL, N_parameters[1], N_parameters[2])$ET
    SWC_N[count] <- preles_test_cpp(nrow(data_format_nit), count, data_format_nit[count,], pPREL, N_parameters[1], N_parameters[2])$SWC
    S_N[count] <- preles_test_cpp(nrow(data_format_nit), count, data_format_nit[count,], pPREL, N_parameters[1], N_parameters[2])$S
    count = count + 1
  }

  par(mfrow = c(2, 2))
  plot(data_format$Date, GPP, xlab = "Date", ylab = "GPP", main = "Normal: GPP", type = "l")
  plot(data_format$Date, ET, xlab = "Date", ylab = "ET", main = "Normal: ET", type = "l")
  plot(data_format$Date, SWC, xlab = "Date", ylab = "SWC", main = "Normal: SWC", type = "l")
  plot(data_format$Date, S, xlab = "Date", ylab = "S", main = "Normal: S", type = "l")

  par(mfrow = c(2, 2))
  plot(Nitrgen, GPP_N, xlab = "Date", ylab = "GPP", main = "Nitrogen Effect: GPP", type = "l")
  plot(Nitrgen, ET_N, xlab = "Date", ylab = "ET", main = "Nitrogen Effect: ET", type = "l")
  plot(Nitrgen, SWC_N, xlab = "Date", ylab = "SWC", main = "Nitrogen Effect: SWC", type = "l")
  plot(Nitrgen, S_N, xlab = "Date", ylab = "S", main = "Nitrogen Effect: S", type = "l")
}
