# ####################################################################################################
# # Convert CETRAD precip data to textfiles
# # Author: Drew Gower
# # Date: January 23, 2016
# ####################################################################################################
## Clear variables and load required packages
rm(list=ls(all=TRUE))
library(gdata)
library(lubridate)
library(data.table)


## Import and manipulate precip data
# Set Working Directory
setwd("~/Dropbox/Kenya Project/Data/Precipitation/CETRAD")

# Define seasons
sref <- data.table(mnth = seq(1, 12), seas = c(rep("S1", 2), rep("S2", 3), 
                                               rep("S3", 4), rep("S4", 3)))

# Read in excel sheet and convert dates and station names
pre.dt <- data.table(read.xls("170116_cetradpre.xlsx", sheet = 1))
pre.dt$RDate <- as.Date(as.character(pre.dt$RDate), format = "%Y-%m-%d", 
                        tz = "Africa/Nairobi")
nms <- colnames(pre.dt)[2:ncol(pre.dt)]
nms <- gsub("\\.\\.", "\\.(", nms)
nms <- gsub("\\.$", ")", nms)
nms <- gsub("\\.", " ", nms)
colnames(pre.dt) <- c("Date", nms)
cor.dt <- copy(pre.dt)

# Remove outlier values
mns <- as.numeric(pre.dt[, lapply(.SD, function(x) {mean(x[x > 0], na.rm = T)}), 
                        .SDcols = nms])
tvl <- -log(1 / (100 * 365)) * mns
for(i in 1:length(nms)) {
  cor.dt[, (nms[i]) := lapply(.SD, function(x) {replace(x, x > tvl[i], NA)}), 
                   .SDcols = nms[i]]
}


## Create summary table of yearly and seasonal values
# Create table structure
sea.dt <- melt(cor.dt, id = 1, measure.vars = 2:ncol(pre.dt), 
               variable.name = "Station", value.name = "Precip")
sea.dt[, `:=` (Year = as.numeric(format(sea.dt$Date, '%Y')),
               Seas = sref$seas[as.numeric(format(sea.dt$Date, '%m'))])]
sea.dt[, Date := NULL]
setcolorder(sea.dt, c("Year", "Seas", "Station", "Precip"))
ann.dt <- copy(sea.dt)

# Calculate annual values
ann.dt[, Seas := "All"]
ann.dt[, `:=` (Ntot = sum(!is.na(Precip)), 
               Npos = sum(Precip > 0, na.rm = T),
               Sum = sum(Precip, na.rm = T),
               Alpha = mean(Precip[which(Precip > 0)]),
               Lambda = sum(Precip > 0, na.rm = T) / sum(!is.na(Precip))), 
       by = c("Station", "Year")] 
ann.dt[, `:=` (Npos = replace(Npos, Ntot == 0, NA),
               Sum = replace(Sum, Ntot == 0, NA),
               Alpha = replace(Alpha, Ntot == 0, NA),
               Lambda = replace(Lambda, Ntot == 0, NA))]
ann.dt[, Alpha := replace(Alpha, is.nan(Alpha), 0)]
ann.dt <- unique(ann.dt[, Precip := NULL])

# Calculate seasonal values
sea.dt[, `:=` (Ntot = sum(!is.na(Precip)), 
               Npos = sum(Precip > 0, na.rm = T),
               Sum = sum(Precip, na.rm = T),
               Alpha = mean(Precip[which(Precip > 0)]),
               Lambda = sum(Precip > 0, na.rm = T) / sum(!is.na(Precip))), 
       by = c("Station", "Seas", "Year")] 
sea.dt[, `:=` (Npos = replace(Npos, Ntot == 0, NA),
               Sum = replace(Sum, Ntot == 0, NA),
               Alpha = replace(Alpha, Ntot == 0, NA),
               Lambda = replace(Lambda, Ntot == 0, NA))]
sea.dt[, Alpha := replace(Alpha, is.nan(Alpha), 0)]
sea.dt <- unique(sea.dt[, Precip := NULL])

# Combine annual and seasonal values
sea.dt <- setorder(rbind(sea.dt, ann.dt), Station, Year, Seas)


## Write data to text files
setwd("~/Dropbox/Kenya Project/Data/Precipitation/CETRAD/Output")

# Retrieve individual data sets and print
for(i in 1:length(nms)) {
  rec.dt <- cor.dt[, .SD, .SDcols = c("Date", nms[i]) ]
  colnames(rec.dt) <- c("Date", "Precip")
  write.table(rec.dt, paste(nms[i], "txt", sep="."), sep="\t", 
              row.names = F, col.names = F)
}

write.table(sea.dt, "summary.txt", sep="\t", row.names = F, col.names = T)
