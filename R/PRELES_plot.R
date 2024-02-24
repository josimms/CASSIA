PRELES_plot <- function(data_format, N_parameters) {

  # install.packages("~/Documents/Rprebasso-master/",repos=NULL, type="source")

  # devtools::install_github("ForModLabUHel/Rprebasso")
  # library(Rprebasso)

  data_format_nit_const = cbind(data_format, Nitrogen = rep(0.012, nrow(data_format)))

  original <- Rprebasso::PRELES(TAir = data_format_nit_const$T,
                                Precip = data_format_nit_const$Rain,
                                CO2 = data_format_nit_const$CO2,
                                PAR = data_format_nit_const$PAR,
                                VPD = data_format_nit_const$VPD,
                                fAPAR = data_format_nit_const$fAPAR,
                                p = pPREL, returncols = c("GPP", "ET", "SW", "fS", "fW", "fE"))


  cpp_model <- NULL
  N_parameters = c(1/0.012, 0)
  cpp_model <- preles_test_cpp(2005, 2023,
                               data_format_nit_const,
                               c(pPREL, N_parameters[1], N_parameters[2]),
                               0)

  # TODO: why is the new model deleting one data point?

  par(mfrow = c(2, 2))
  plot(data_format$Date, original$GPP, xlab = "Day of the Year", ylab = "GPP", main = "GPP", type = "l")
  lines(data_format$Date, cpp_model$GPP, col = "blue")
  plot(data_format$Date, original$ET, xlab = "Day of the Year", ylab = "ET", main = "ET", type = "l")
  lines(data_format$Date, cpp_model$ET, col = "blue")
  plot(data_format$Date, original$SW, xlab = "Day of the Year", ylab = "SWC", main = "SWC", type = "l")
  lines(data_format$Date, cpp_model$SWC, col = "blue")
  plot(data_format$Date, original$fS, xlab = "Day of the Year", ylab = "fS", main = "S", type = "l")
  lines(data_format$Date, cpp_model$fS, col = "blue")

  plot(data_format$Date, original$fW, xlab = "Day of the Year", ylab = "fW", main = "fW", type = "l")
  lines(data_format$Date, cpp_model$fW, col = "blue") # TODO: why doesn't this exist?
  plot(data_format$Date, original$fE, xlab = "Day of the Year", ylab = "fE", main = "fE", type = "l")
  lines(data_format$Date, cpp_model$fE, col = "blue")
  plot(data_format$Date, cpp_model$fN, xlab = "Day of the Year", ylab = "fN", main = "fN (only CPP)",
       type = "l", col = "blue")
  plot(data_format$Date, original$fS, xlab = "Day of the Year", ylab = "fS", main = "fS", type = "l")
  lines(data_format$Date, cpp_model$fS, col = "blue")

  par(mfrow = c(4, 2))
  plot(cpp_model$GPP, original$GPP, xlab = "GPP Nitrogen", ylab = "GPP", main = "Normal: GPP")
  abline(0, 1, col = "red")
  plot(data_format$Date, cpp_model$GPP - original$GPP, xlab = "Dates", ylab = "Difference", main = "Residuals")
  plot(cpp_model$ET, original$ET, xlab = "ET Nitrogen", ylab = "ET", main = "Normal: ET")
  abline(0, 1, col = "red")
  plot(data_format$Date, cpp_model$ET - original$ET, xlab = "Dates", ylab = "Difference", main = "Residuals")
  plot(cpp_model$SWC, original$SWC, xlab = "SWC Nitrogen", ylab = "SWC", main = "Normal: SWC")
  abline(0, 1, col = "red")
  plot(data_format$Date, cpp_model$SWC - original$SW, xlab = "Dates", ylab = "Difference", main = "Residuals")
  plot(cpp_model$fS, original$S, xlab = "fS Nitrogen", ylab = "fS", main = "Normal: S")
  abline(0, 1, col = "red")
  plot(data_format$Date, cpp_model$fS - original$S, xlab = "Dates", ylab = "Difference", main = "Residuals")

  par(mfrow = c(2, 2))
  plot(data_format_nit_const$Nitrogen, cpp_model$GPP, xlab = "Nitrogen", ylab = "GPP", main = "Nitrogen Effect: GPP")
  plot(data_format_nit_const$Nitrogen, cpp_model$ET, xlab = "Nitrogen", ylab = "ET", main = "Nitrogen Effect: ET")
  plot(data_format_nit_const$Nitrogen, cpp_model$SWC, xlab = "Nitrogen", ylab = "SWC", main = "Nitrogen Effect: SWC")
  plot(data_format_nit_const$Nitrogen, cpp_model$fN, xlab = "Nitrogen", ylab = "fN", main = "Nitrogen Effect: N")
}

