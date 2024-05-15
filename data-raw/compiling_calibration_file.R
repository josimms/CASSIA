####
# Directory and function set up
####

raw_data_directory <- "~/Documents/CASSIA_Calibration/Raw_Data/"
processes_data_directory <- "~/Documents/CASSIA_Calibration/Processed_Data/"

## Reading the data and printing the file for debugging and sorting out the units
excel_Reading_function <- function(name, sheet) {
  print(name)
  out_data = readxl::read_excel(name, sheet = sheet, na = "-")
  return(out_data)
}

## Reference adder function
reference_adder <- function(dataframe, string) {
  dataframe$Reference = string
  return(dataframe)
}

###
# Nitrogen data and litter data
###

# Inport
# TODO: This is reading in wrong, but sort out when plotting the data
nitrogen_balence_data_raw_lys_mökki = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-c(1, 2)]), excel_Reading_function, sheet = "lys mökki")
# lys mökki has the NH4, NO3 etc.
nitrogen_balence_data_raw_Lumi = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "Lumi")
nitrogen_balence_data_raw_oksa = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "oksa")
nitrogen_balence_data_raw_karike = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "karike")
nitrogen_balence_data_raw_pato = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "Pato")
nitrogen_balence_data_raw_keskisade = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "keskisade")
nitrogen_balence_data_raw_mökki_runko = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "mökki runko")

## TODO: work out exactly which ones of these are actually useful! Import extra ones if needed!
# nitrogen_balence_data_raw_tensiometrit = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "Tensiometrit")
# nitrogen_balence_data_raw_sade_mökki = lapply(paste0(raw_data_directory, "ancillary_ecosystem/", list.files(paste0(raw_data_directory, "ancillary_ecosystem/"))[-1]), excel_Reading_function, sheet = "sade mökki")

# Add reference!
nitrogen_balence_data_raw_lys_mökki_referenced = lapply(nitrogen_balence_data_raw_lys_mökki, reference_adder, "ancillary_ecosystem")
nitrogen_balence_data_raw_Lumi_referenced = lapply(nitrogen_balence_data_raw_Lumi, reference_adder, "ancillary_ecosystem")
nitrogen_balence_data_raw_oksa_referenced = lapply(nitrogen_balence_data_raw_oksa, reference_adder, "ancillary_ecosystem")
nitrogen_balence_data_raw_karike_referenced = lapply(nitrogen_balence_data_raw_karike, reference_adder, "ancillary_ecosystem")
nitrogen_balence_data_raw_pato_referenced = lapply(nitrogen_balence_data_raw_pato, reference_adder, "ancillary_ecosystem")
nitrogen_balence_data_raw_keskisade_referenced = lapply(nitrogen_balence_data_raw_keskisade, reference_adder, "ancillary_ecosystem")
nitrogen_balence_data_raw_mökki_runko_referenced = lapply(nitrogen_balence_data_raw_mökki_runko, reference_adder, "ancillary_ecosystem")

# Change the names of the columns so they are consistent
nitrogen_balence_data_raw_karike_dataframe = data.table::rbindlist(nitrogen_balence_data_raw_karike_referenced, fill = T)
# TODO: get the data to imports
nitrogen_balence_data_raw_lys_mökki_dataframe = data.table::rbindlist(nitrogen_balence_data_raw_lys_mökki_referenced, fill = T)

### Make averages of the litter data
# Taking out averages
nitrogen_balence_data_raw_karike_dataframe_karikekeräin <- nitrogen_balence_data_raw_karike_dataframe[!is.na(nitrogen_balence_data_raw_karike_dataframe$karikekeräin),]
# Summarising by day
nitrogen_balence_data_raw_karike_dataframe_daily <- aggregate(nitrogen_balence_data_raw_karike_dataframe_karikekeräin[,c(3:8, 10, 14:17)],
                                                              nitrogen_balence_data_raw_karike_dataframe_karikekeräin[,c("pvm.")],
                                                              mean, na.rm = T)

###
# Plot the data!
###

plot(nitrogen_balence_data_raw_karike_dataframe_daily$pvm.,
     nitrogen_balence_data_raw_karike_dataframe_daily$oksa,
     main = "Biomass Hyytiälä",
     xlab = "Date",
     ylab = "Biomass, g m-2 per irregular time interval",
     pch = 1)
for (i in 2:ncol(nitrogen_balence_data_raw_karike_dataframe_daily)) {
  points(nitrogen_balence_data_raw_karike_dataframe_daily$pvm.,
         nitrogen_balence_data_raw_karike_dataframe_daily[,c(i)],
         pch = i,
         col = i)
}
legend("topright",
       names(nitrogen_balence_data_raw_karike_dataframe_daily)[-1],
       bty = "n",
       pch = 2:nrow(nitrogen_balence_data_raw_karike_dataframe_daily),
       col = 2:nrow(nitrogen_balence_data_raw_karike_dataframe_daily),
       ncol = 2,
       cex = 0.75)

###
# Check biomass with Jussi's data!
###

processed_data_directory <- "~/Documents/CASSIA_Calibration/Processed_Data/"
jussi_biomass_data <- readxl::read_excel(paste0(processed_data_directory, "Jussi_data.xlsx"))

plot(jussi_biomass_data$date[17:nrow(jussi_biomass_data)],
     jussi_biomass_data$amount[17:nrow(jussi_biomass_data)],
     col = as.factor(jussi_biomass_data$variable[17:nrow(jussi_biomass_data)]),
     pch = as.numeric(as.factor(jussi_biomass_data$variable[17:nrow(jussi_biomass_data)])),
     main = "Jussi Biomass Data - used as a calculation check",
     xlab = "Date",
     ylab = "Biomass")

# Numerically check - it's the same!!! :D
summary(jussi_biomass_data$amount[jussi_biomass_data$variable == "Average of oksat"])
summary(nitrogen_balence_data_raw_karike_dataframe_daily$oksa[jussi_biomass_data$date %in% nitrogen_balence_data_raw_karike_dataframe_daily$pvm.])

###
# Growth data!
###

#
processed_data_directory <- "~/Documents/CASSIA_Calibration/Processed_Data/"
SMEARII_data <- readxl::read_excel(paste0(processed_data_directory, "smearII_data.xlsx"))

plot(SMEARII_data$date[grepl("LAI", SMEARII_data$variable, fixed=TRUE)],
     SMEARII_data$amount[grepl("LAI", SMEARII_data$variable, fixed=TRUE)],
     main = "LAI",
     col = as.numeric(as.factor(gsub("_ICOS", "", gsub("_allsided", "", SMEARII_data$variable[grepl("LAI", SMEARII_data$variable, fixed=TRUE)])))),
     pch = as.numeric(grepl("ICOS", SMEARII_data$variable[grepl("LAI", SMEARII_data$variable, fixed=TRUE)], fixed=TRUE))+2,
     ylab = "LAI",
     xlab = "Dates")
legend("topleft",
       c("Allsided", "ICOS"),
       pch = 2:3,
       bty = "n")
legend("bottomleft",
       levels(as.factor(gsub("_ICOS", "", gsub("_allsided", "", SMEARII_data$variable[grepl("LAI", SMEARII_data$variable, fixed=TRUE)])))),
       col = 1:3, # TODO: sort out the colours!
       bty = "n",
       pch = "+")

# Height and DBH

plot(SMEARII_data$date[grepl("Arithmetic mean height", SMEARII_data$variable, fixed=TRUE)],
     SMEARII_data$amount[grepl("Arithmetic mean height", SMEARII_data$variable, fixed=TRUE)],
     main = "Arithmetic mean height",
     col = as.numeric(as.factor(SMEARII_data$variable[grepl("Arithmetic mean height", SMEARII_data$variable, fixed=TRUE)])),
     pch = 2,
     ylab = "Arithmetic mean height, m",
     xlab = "Dates")
legend("topleft",
       levels(as.factor(gsub("Arithmetic mean height", "", SMEARII_data$variable[grepl("Arithmetic mean height", SMEARII_data$variable, fixed=TRUE)]))),
       col = 1:3,
       bty = "n",
       pch = 2)

plot(SMEARII_data$date[grepl("Arithmetic mean DBH", SMEARII_data$variable, fixed=TRUE)],
     SMEARII_data$amount[grepl("Arithmetic mean DBH", SMEARII_data$variable, fixed=TRUE)],
     main = "Arithmetic mean DBH",
     col = as.numeric(as.factor(SMEARII_data$variable[grepl("Arithmetic mean DBH", SMEARII_data$variable, fixed=TRUE)])),
     pch = 2,
     ylab = "Arithmetic mean DBH, cm",
     xlab = "Dates")
legend("topleft",
       levels(as.factor(gsub("Arithmetic mean DBH", "", SMEARII_data$variable[grepl("Arithmetic mean DBH", SMEARII_data$variable, fixed=TRUE)]))),
       col = 1:3, # TODO: sort out the colours!
       bty = "n",
       pch = 2)

####
# Reading in a csv file
####

smearII_data = read.csv(paste0(processes_data_directory, "/smearII.csv"))
smearII_data$pmv_kaikki = as.Date(paste0(smearII_data$Unnamed..1, smearII_data$pvm))
smearII_data$TOT.N.mg.l = paste0(smearII_data$TOT.N..mg.l., smearII_data$TOT.N..mg.l..1)
smearII_data$NO3.N_0.01[is.na(smearII_data$NO3.N_0.01)] = FALSE
smearII_data$NO3.N_0.001[is.na(smearII_data$NO3.N_0.001)] = FALSE
smearII_data$In.meaurement..NO3.N = smearII_data$NO3.N_0.01 + smearII_data$NO3.N_0.001
smearII_data$NO3.N.mg.l = paste0(smearII_data$NO3.N..mg.l., smearII_data$NO3.N..mg.l..1)
smearII_data <- smearII_data[,-which(names(smearII_data) %in% c("X", "Unnamed..1", "Unnamed..5", "Unnamed..6", "pvm", "TOT.N..mg.l.", "TOT.N..mg.l..1", "NO3.N..mg.l.", "NO3.N..mg.l..1"))]
names(smearII_data)[1] = c("PIT")
smearII_data$pmv_kaikki_monthly = as.Date(paste(2024, substring(smearII_data$pmv_kaikki, 6, 10), sep = "-"))

### Yearly graph

par(mfrow = c(3, 1))
plot(smearII_data$pmv_kaikki, as.numeric(smearII_data$NH4..µg.l.),
     xlab = "Date", ylab = "NH4, ug/l", ylim = c(0, 50), main = "Nitrogen Balance")
points(smearII_data$pmv_kaikki, as.numeric(smearII_data$NH4.N..mg.l.)*1.28786/0.001,
       pch = as.logical(smearII_data$In.measurement..NH4.N..mg.l.)+1)
legend("topleft", c("Measurements", "Measurement Sensitivity"), pch = c(1, 2), bty = "n")

plot(smearII_data$pmv_kaikki, as.numeric(smearII_data$NO3.N.mg.l)*4.42664/0.001,
     xlab = "Date", ylab = "NO3, ug/l")
points(smearII_data$pmv_kaikki, smearII_data$nitraatti/0.001, col = "blue")
points(smearII_data$pmv_kaikki, as.numeric(smearII_data$NO3.N)*4.42664/0.001,
       pch = smearII_data$In.meaurement..NO3.N+1, col = "red")
legend("topleft", c("Measurements", "Measurement Sensitivity"), pch = c(1, 2), bty = "n")

plot(smearII_data$pmv_kaikki, smearII_data$TOT.N.mg.l,
     xlab = "Date", ylab = "Norg")
points(smearII_data$pmv_kaikki, as.numeric(smearII_data$TN..mg.l.))

### Monthly graph

par(mfrow = c(2, 1))
plot(smearII_data$pmv_kaikki_monthly, as.numeric(smearII_data$NH4..µg.l.),
     xlab = "Date", ylab = "NH4, ug/l", ylim = c(0, 50), main = "Nitrogen Balance")
points(smearII_data$pmv_kaikki_monthly, as.numeric(smearII_data$NH4.N..mg.l.)*1.28786/0.001,
       pch = as.logical(smearII_data$In.measurement..NH4.N..mg.l.)+1)
legend("topleft", c("Measurements", "Measurement Sensitivity"), pch = c(1, 2), bty = "n")

plot(smearII_data$pmv_kaikki_monthly, as.numeric(smearII_data$NO3.N.mg.l)*4.42664/0.001,
     xlab = "Date", ylab = "NO3, ug/l")
points(smearII_data$pmv_kaikki_monthly, smearII_data$nitraatti/0.001, col = "blue")
points(smearII_data$pmv_kaikki_monthly, as.numeric(smearII_data$NO3.N)*4.42664/0.001,
       pch = smearII_data$In.meaurement..NO3.N+1, col = "red")
legend("topleft", c("Measurements", "Measurement Sensitivity"), pch = c(1, 2), bty = "n")

