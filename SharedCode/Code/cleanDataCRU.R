# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# Read in the netcdf file
handle <- nc_open(sprintf('Data/Raw/%s',cruFile))

lon <- ncvar_get(handle,'lon')
lat <- ncvar_get(handle,'lat')

# set up all of the necessary time objects to properly subset the data
time <- ncvar_get(handle,'time')

tFull <- seq(1901,2016+1,length=length(time)+1)[1:length(time)]
year <- floor(tFull)
month <- rep(1:12,length(time)/12)

fullTimeMat <- cbind(year,month)

# subset based on starting year and make all of the necessary time objects
timeInds <- which(year >= cruStartYear & year <= cruEndYear)

timeMat <- fullTimeMat[timeInds,]

# grab the ppt monthly time series
ppt <- ncvar_get(handle, 'pre', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )

counts <- ncvar_get(handle, 'stn', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )
nc_close(handle)

# deal with NAs
counts[counts==-999] <- NA


# Get the needed year vector for the seasonal function
tYear <- cruStartYear:cruEndYear

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

save(observationList, file='Data/RawProcessed/cruSeasonal_0.5.Rda')