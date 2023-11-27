respiration_plot <- function() {

  B0 = 20

  RmN_daily <- RmS_daily <- RmR_daily <- Rm_daily <- NULL
  for (day in 1:nrow(data_format)) {
    RmN_daily[day] = respiration_test_cpp(t(parameters_p), common_p, t(ratios_p), t(sperling_p), sperling_extras,
                                           day,
                                           data_format$T[day], data_format$TSA[day],
                                           temp_rise, Rm_acclimation, mN_varies,
                                           B0)$RmN
    RmS_daily[day] = respiration_test_cpp(t(parameters_p), common_p, t(ratios_p), t(sperling_p), sperling_extras,
                                          day,
                                          data_format$T[day], data_format$TSA[day],
                                          temp_rise, Rm_acclimation, mN_varies,
                                          B0)$RmS
    RmR_daily[day] = respiration_test_cpp(t(parameters_p), common_p, t(ratios_p), t(sperling_p), sperling_extras,
                                          day,
                                          data_format$T[day], data_format$TSA[day],
                                          temp_rise, Rm_acclimation, mN_varies,
                                          B0)$RmR
    Rm_daily[day] = respiration_test_cpp(t(parameters_p), common_p, t(ratios_p), t(sperling_p), sperling_extras,
                                          day,
                                          data_format$T[day], data_format$TSA[day],
                                          temp_rise, Rm_acclimation, mN_varies,
                                          B0)$Rm_a
  }

  ###
  # Pure results
  ###

  par(mfrow = c(2, 1), xpd=TRUE)
  plot(data_format$Date, Rm_daily, type = "l",
       main = "Respiration Output", ylab = "", xlab = "Date", col = "blue")
  lines(data_format$Date, RmS_daily, col = "brown")
  legend(data_format$Date[1], y = -max(Rm_daily)/2, c("Wood", "All"),
         col = c("brown", "blue"), lty = c(1, 1), bty = "n")

  plot(data_format$Date, RmN_daily, col = "green", type = "l",
       main = "Respiration Output", ylab = "", xlab = "Date")
  lines(data_format$Date, RmR_daily, col = "black")
  legend(data_format$Date[1], y = -max(RmN_daily)/2, c("Neeldes", "Roots"),
         col = c("green", "black"), lty = c(1, 1), bty = "n")
}
