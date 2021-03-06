# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')


# Link the raw netcdf file
handle <- nc_open(sprintf('Data/Raw/%s',gpcpFile))

# Lon in the grid and data and flip to the format I'm working with
lonRaw <- ncvar_get(handle,'X')
lon <- lonRaw - 180

latRaw <- ncvar_get(handle,'Y')
lat <- rev(latRaw)

pptRaw <- ncvar_get(handle,'precip')

midLon <- length(lon)/2
pptMoveMed <- pptRaw[c((midLon+1):length(lon),1:midLon),,]

ppt <- array(NA, dim=dim(pptRaw))
for(i in 1:nrow(pptRaw)){
	for(k in 1:dim(pptRaw)[3]){
		ppt[i,,k] <- rev(pptMoveMed[i,,k])
	}
}



# deal with the time variable (code is for full years of data)
timeRaw <- ncvar_get(handle,'T')

startYear <- 1979
startMonth <- 1

nYears <- length(timeRaw)/12
endYear <- startYear+nYears -1

timeMap <- cbind(rep(startYear:endYear,each=12),
				 rep(1:12, nYears))
timeMap <- cbind(timeMap,timeMap[,1] + (timeMap[,2]-1)/12)


# convert ppt to monthly values from (mm/day)

startDate <- as.Date("1979-01-01", "%Y-%m-%d")
endDate   <- as.Date("2017-12-01", "%Y-%m-%d")

numDays <- days_in_month(seq(startDate,endDate,'month'))

monthlyPrecip <- array(NA, dim=dim(ppt))

for(i in 1:dim(monthlyPrecip)[3]){
	monthlyPrecip[,,i] <- ppt[,,i]*numDays[i]
}

# convert ppt to seasonal sums

seasonalPrecip <- seasonalField12(monthlyPrecip,startYear:endYear,sum)


# mask to land locations
load('Data/RawProcessed/cruSeasonal_2.5.Rda')
landMask <- ifelse(is.na(observationList$ppt[,,1]),NA,1)
rm(observationList)

seasonalPrecipLand <- array(NA, dim=dim(seasonalPrecip))

for(i in 1:dim(seasonalPrecip)[3]){
	seasonalPrecipLand[,,i] <- seasonalPrecip[,,i]*landMask
}

timeMapFinal <- timeMap[-(1:12),]

observationList <- list(lon=lon,
						lat=lat,
						observations=seasonalPrecipLand,
						timeMap = timeMapFinal)

save(observationList,file=sprintf("Data/RawProcessed/gpcpSeasonal_2.5.Rda"))

