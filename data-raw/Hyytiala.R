library(Tandem)
Hyde_weather_CASSIA <- Tandem::loadRData(paste0("/home/joanna/Asiakirjat/Hyytiälä/", "Hyde_weather_CASSIA.RData"))

# Gapfill with FMI
summary(Hyde_weather_CASSIA)
