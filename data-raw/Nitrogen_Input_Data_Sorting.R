direct = "~/Documents/CASSIA_Calibration/Processed_Data/"

###
# Read data
###
nitrogen = read.csv(paste0(direct, "korhonen_mineral_N.csv"), dec = ",")
nitrogen$date = as.Date(nitrogen$date, format = "%d.%m.%Y")
save(nitrogen, file = paste0(direct, "nitrogen.RData"))

###
# Unit conversions. g N m-2 will be the input
###
# Assuming the first 10cm of the soil
# 1200 kg / m3 seems to be a rule of thumb for soil
# For the m3 to m2 conversion assume that there is a constant a

# mg/kg = 0.001 g / kg = 0.001 g / (m3 / 1200) = 1200 * 0.001 g / m3 = 1200 * 0.001 * 0.01 g / m2

