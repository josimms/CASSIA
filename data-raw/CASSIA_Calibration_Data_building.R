####
# Read the data
#
# The data was generated from the original model that Pauliina sent to me, and should be used for calibration purposes only!
####
Hyde_daily_original = read.table("./data-raw/daily_Hyde.txt", sep = " ", header = T)
Hyde_yearly_original = read.table("./data-raw/yearly_Hyde.txt", sep = " ", header = T)

####
# Save as a Rdata file so can be used easily in the package
####

Hyde_daily_original$day[Hyde_daily_original$year == 2016 & Hyde_daily_original$day > 59] <- 61:366

Hyde_daily_original$date <- as.Date(as.Date(Hyde_daily_original$day,
                                            origin = as.Date(paste0(Hyde_daily_original$year,
                                                                    "-01-01"))))

save(Hyde_daily_original, file = "./data/Hyde_daily_original.RData")
save(Hyde_yearly_original, file = "./data/Hyde_yearly_original.RData")

####
# Reading the CASSIA preles data
####


