# Read in the netcdf file for ppt
handle <- nc_open(sprintf('Data/Raw/%s',gpccPptFile))

lon <- ncvar_get(handle,'X')
latRaw <- ncvar_get(handle,'Y')

lat <- rev(latRaw)

# set up all of the necessary time objects to properly subset the data
time <- ncvar_get(handle,'T')

tFull <- seq(1950,2013+1,length=length(time)+1)[1:length(time)]
year <- floor(tFull)
month <- rep(1:12,length(time)/12)

fullTimeMat <- cbind(year,month)

# subset based on starting year and make all of the necessary time objects
timeInds <- which(year >= gpccStartYear & year <= gpccEndYear)

timeMat <- fullTimeMat[timeInds,]

# grab the ppt monthly array
pptRaw <- ncvar_get(handle, 'prcp', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )

# flip around to have same arrangement as other data
yInds <- ncol(pptRaw):1

ppt <- array(NA, dim=dim(pptRaw))

for(i in 1:nrow(pptRaw)){
	for(k in 1:length(time)){
		ppt[i,,k] <- rev(pptRaw[i,,k])
	}
}

nc_close(handle)

# build a land mask to mask the counts data
landMask <- ifelse(is.na(ppt[,,1]),NA,1)

# Do the same for the station counts data
handle2 <- nc_open(sprintf('Data/Raw/%s',gpccStnFile))

countsRaw <- ncvar_get(handle2, 'cnts', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )

# flip around to have same arrangement as other data
yInds <- ncol(pptRaw):1

countsFull <- array(NA, dim=dim(pptRaw))

for(i in 1:nrow(countsRaw)){
	for(k in 1:length(time)){
		countsFull[i,,k] <- rev(countsRaw[i,,k])
	}
}

counts <- array(NA, dim=dim(pptRaw))

for(k in 1:length(time)){
	counts[,,k] <- countsFull[,,k]*landMask
}

nc_close(handle2)

# Get the needed year vector for the seasonal function
tYear <- gpccStartYear:gpccEndYear

# Take seasonal statistics
pptSeasonal <- seasonalField12(ppt,tYear,sum)
countsSeasonal <- seasonalField12(counts,tYear,min)

# get the correct time objects
timeMapWorking <- timeMat[13:(nrow(timeMat)),]
timeMap <- cbind(timeMapWorking,timeMapWorking[,1] + (timeMapWorking[,2]-1)/12)

# build the lists with all information (save as the generic 'observation list' to
# allow for code to easily run with multiple ppt datasets)
observationList <- list(lon=lon,lat=lat,ppt=pptSeasonal,
					counts=countsSeasonal,timeMap=timeMap)

save(observationList, file=sprintf('Data/RawProcessed/%s',gpccOutFileName))