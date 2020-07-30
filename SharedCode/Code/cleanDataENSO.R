# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')


ensoMonthlyFile  <- 'cpc_nino34.txt'
ensoSeasonalFile <- 'cpc_nino34_seasonal.txt'
enoOniFile       <- 'cpc_oni.txt'

ninoMonthly  <- read.table(sprintf('Data/Raw/%s', ensoMonthlyFile), 
					header=TRUE)
ninoSeasonal <- read.table(sprintf('Data/Raw/%s', ensoSeasonalFile),
					header=FALSE, stringsAsFactors=FALSE)
oniSeasonal  <- read.table(sprintf('Data/Raw/%s', enoOniFile),
					header=TRUE, stringsAsFactors=FALSE)
#
# ONI from https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt

# first clean the monthly series
timeMap <- as.matrix(ninoMonthly[,1:2])
colnames(timeMap) <- NULL
timeMap <- cbind(timeMap, timeMap[,1] + (timeMap[,2]-1)/12)

# remove 2020 (incomplete year)
inds <- which(timeMap[,1]<2020)
anomSeries <- ninoMonthly$ANOM[inds]
timeMap <- timeMap[inds,]

mySeasonal <- seasonalAverage12(ninoMonthly$ANOM,unique(timeMap[,1]))


# Now clean the seasonal series (already cleaned through 2019, remove 1950)
cpcSeasonal <- ninoSeasonal[-(1:12),3]


# Finally, clean the ONI (also already cleaned through 2019, remove 1950)
cpcOni <- oniSeasonal[-(1:12),4]

# final clean of the time map
timeMapFinal <- timeMap[-c(1:12),]

# plot up comparisons
yMax <- max(abs(c(mySeasonal, cpcSeasonal, cpcOni)))

plot(timeMapFinal[,3],c(t(mySeasonal)),ylim=c(-yMax,yMax),type='l',col='black')
points(timeMapFinal[,3],cpcSeasonal,type='l',col='blue')
points(timeMapFinal[,3],cpcOni,type='l',col='red')
grid()
abline(h=c(-0.5,0.5),lty=3)
abline(h=0)

# difference series
plot(timeMapFinal[,3], c(t(mySeasonal)) - cpcSeasonal, col='blue', type='l')
points(timeMapFinal[,3], c(t(mySeasonal)) - cpcOni, col='red', type='l')
grid()
abline(h=0)

# Instinct is to use the ONI: same as the monthly Nino3.4 averaged to seasonal
# using a climatology that updates every 5 years. The only difference from the
# WMO definition will be (potentially) the climatology update rate


ninaSeasonal <- list(nina34 = cpcOni, timeMap = timeMapFinal)

save(ninaSeasonal, file='Data/RawProcessed/ninaSeasonal.Rda')