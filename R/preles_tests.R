test_PRELES_N_function <- function() {
  
  output <- list()
  nitrogen_range <- seq(0, 10, by = 1)
  count = 1
  for (n in nitrogen_range) {
    output[[count]] <- Rprebasso::PRELES(data_format$PAR, 
                                        data_format$T, 
                                        data_format$VPD, 
                                        data_format$Rain,
                                        data_format$CO2,
                                        data_format$fAPAR, 
                                        rep(n, length = nrow(data_format)))$GPP
    count = count + 1 
  }
  
  par(mfrow = c(1, 1))
  plot(nitrogen_range, unlist(lapply(output, sum, na.rm = T)), 
       xlab = "Nitrogen Input, constant",
       ylab = "Sum of photosynthesis per year")
  
}