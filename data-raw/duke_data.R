duke_directory = "~/Documents/CASSIA_Calibration/Raw_Data/duke_data/data/"

###
# Functions
###

make_date <- function(x) {
  if (sum(names(x) %in% "DOY")) {
    x$Date = as.Date(paste(x$Year, x$DOY, sep = "-"), format = "%Y-%j")
  }
  return(x)
}

###
# Read in data
###

duke_data <- lapply(paste0(duke_directory, list.files(paste0(duke_directory), pattern = "*.csv")), read.csv)
names(duke_data) <- gsub(".csv", "", gsub("DukeFACE_", "", list.files(paste0(duke_directory), pattern = "*.csv")))
duke_data <- lapply(duke_data, make_date)

###
# Plot data
###

par(mfrow = c(2, 1))
plot(duke_data$DukaFACE_LAI$Date[duke_data$DukaFACE_LAI$N == "CONT"], 
     duke_data$DukaFACE_LAI$LAI_total[duke_data$DukaFACE_LAI$N == "CONT"], 
     xlab = "Date", ylab = "LAI", main = "LAI", 
     type = "l", ylim = c(1, 9.5)) # TODO: colours
lines(duke_data$DukaFACE_LAI$Date[duke_data$DukaFACE_LAI$N != "CONT"], 
      duke_data$DukaFACE_LAI$LAI_total[duke_data$DukaFACE_LAI$N != "CONT"],
      col = 2)
legend("topleft", c("Control", "Fertilised"), col = c(1, 2), pch = 15, bty = "n")

plot(duke_data$SoilCO2Efflux$Date, duke_data$SoilCO2Efflux$R1C, 
     xlab = "Date", ylab = "CO2 Efflux", main = "CO2 Efflux Control Plots", 
     type = "l")
for (i in 3:ncol(duke_data$SoilCO2Efflux)) {
  lines(duke_data$SoilCO2Efflux$Date, duke_data$SoilCO2Efflux[,i], col = if (i <= 10) {1} else {2})
}

par(mfrow = c(2, 2))
plot(duke_data$NPP$Year[duke_data$NPP$N == "CONT"], 
     duke_data$NPP$NPP_branch_pine[duke_data$NPP$N == "CONT"], 
     xlab = "Date", ylab = "NPP", main = "NPP: Branch", 
     type = "l", 
     ylim = c(min(duke_data$NPP$NPP_branch_pine, na.rm = T),
              max(duke_data$NPP$NPP_branch_pine, na.rm = T)))
legend("topleft", c("Control", "Fertilised"), col = c(1, 2), pch = 15, bty = "n")
lines(duke_data$NPP$Year[duke_data$NPP$N != "CONT"], 
      duke_data$NPP$NPP_branch_pine[duke_data$NPP$N != "CONT"], col = 2)

plot(duke_data$NPP$Year[duke_data$NPP$N == "CONT"], 
     duke_data$NPP$NPP_stem_pine[duke_data$NPP$N == "CONT"],
     xlab = "Year", ylab = "Biomass", main = "NPP: Stem", type = "l",
     ylim = c(min(duke_data$NPP$NPP_stem_pine, na.rm = T), 
              max(duke_data$NPP$NPP_stem_pine, na.rm = T)))
lines(duke_data$NPP$Year[duke_data$NPP$N != "CONT"], 
      duke_data$NPP$NPP_stem_pine[duke_data$NPP$N != "CONT"], col = 2)

plot(duke_data$NPP$Year[duke_data$NPP$N == "CONT"], 
     duke_data$NPP$NPP_foliage_pine[duke_data$NPP$N == "CONT"],
     xlab = "Year", ylab = "Biomass", main = "NPP: Foliage", type = "l",
     ylim = c(min(duke_data$NPP$NPP_foliage_pine, na.rm = T), 
              max(duke_data$NPP$NPP_foliage_pine, na.rm = T)))
lines(duke_data$NPP$Year[duke_data$NPP$N != "CONT"], 
      duke_data$NPP$NPP_foliage_pine[duke_data$NPP$N != "CONT"], col = 2)

plot(duke_data$NPP$Year[duke_data$NPP$N == "CONT"], 
     duke_data$NPP$NPP_coarseroot_pine[duke_data$NPP$N == "CONT"],
     xlab = "Year", ylab = "Biomass", main = "NPP: Coarse Root", type = "l",
     ylim = c(min(duke_data$NPP$NPP_coarseroot_pine, na.rm = T), 
              max(duke_data$NPP$NPP_coarseroot_pine, na.rm = T)))
lines(duke_data$NPP$Year[duke_data$NPP$N != "CONT"], 
      duke_data$NPP$NPP_coarseroot_pine[duke_data$NPP$N != "CONT"], col = 2)


par(mfrow = c(2, 2))
plot(duke_data$Biomass$Year[duke_data$Biomass$N == "CONT"],
     duke_data$Biomass$Biomass_branch_pine[duke_data$Biomass$N == "CONT"], 
     xlab = "Year", ylab = "Biomass", main = "Biomass: Branch", type = "l",
     ylim = c(min(duke_data$Biomass$Biomass_branch_pine, na.rm = T),
              max(duke_data$Biomass$Biomass_branch_pine, na.rm = T)))
legend("topleft", c("Control", "Fertilised"), col = c(1, 2), pch = 15, bty = "n")
lines(duke_data$Biomass$Year[duke_data$Biomass$N != "CONT"],
      duke_data$Biomass$Biomass_branch_pine[duke_data$Biomass$N != "CONT"], col = 2)

plot(duke_data$Biomass$Year[duke_data$Biomass$N == "CONT"], 
     duke_data$Biomass$Biomass_stem_pine[duke_data$Biomass$N == "CONT"],
     xlab = "Year", ylab = "Biomass", main = "Biomass: Stem", type = "l",
     ylim = c(min(duke_data$Biomass$Biomass_stem_pine, na.rm = T),
              max(duke_data$Biomass$Biomass_stem_pine, na.rm = T)))
lines(duke_data$Biomass$Year[duke_data$Biomass$N != "CONT"],
      duke_data$Biomass$Biomass_stem_pine[duke_data$Biomass$N != "CONT"], col = 2)

plot(duke_data$Biomass$Year[duke_data$Biomass$N == "CONT"], 
     duke_data$Biomass$Biomass_foliage_pine[duke_data$Biomass$N == "CONT"],
     xlab = "Year", ylab = "Biomass", main = "Biomass: Foliage", type = "l",
     ylim = c(min(duke_data$Biomass$Biomass_foliage_pine, na.rm = T),
              max(duke_data$Biomass$Biomass_foliage_pine, na.rm = T)))
lines(duke_data$Biomass$Year[duke_data$Biomass$N != "CONT"],
      duke_data$Biomass$Biomass_foliage_pine[duke_data$Biomass$N != "CONT"], col = 2)

plot(duke_data$Biomass$Year[duke_data$Biomass$N == "CONT"],
     duke_data$Biomass$Biomass_root_pine[duke_data$Biomass$N == "CONT"],
     xlab = "Year", ylab = "Biomass", main = "Biomass: Root", type = "l",
     ylim = c(min(duke_data$Biomass$Biomass_root_pine, na.rm = T),
              max(duke_data$Biomass$Biomass_root_pine, na.rm = T)))
lines(duke_data$Biomass$Year[duke_data$Biomass$N != "CONT"],
      duke_data$Biomass$Biomass_root_pine[duke_data$Biomass$N != "CONT"], col = 2)


