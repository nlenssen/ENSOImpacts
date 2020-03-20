# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# load the data and do the usual iridl cleaning process

handle <- nc_open(sprintf('Data/Raw/%s',newFile))

# load in all the data (aleady in proper format!!)
lon <- ncvar_get(handle,'X')
lat <- ncvar_get(handle,'Y')

ppt <- ncvar_get(handle,'prcp')

# deal with the time variable (code is hardcoded since the data will
# not get longer or shorter since it is a deprecated dataset)

tYear <- newStartYear:newEndYear

timeMap <- cbind(rep(tYear,each=12),rep(1:12,length(tYear)),NA)
timeMap[,3] <- timeMap[,1] + (timeMap[,2]-1)/12


# make the data into seaasonal sums!
seasonalPrecip <- seasonalField12(ppt,tYear,sum)


# we make a fake counts array as the New dataset did not provide this quantity
counts <- seasonalPrecip
counts[!is.na(seasonalPrecip)] <- 8

# clean up the time map and save everything
timeMapFinal <- timeMap[-(1:12),]

observationList <- list(lon=lon,lat=lat,ppt=seasonalPrecip,
					counts=counts,timeMap=timeMap)

save(observationList,file=sprintf("Data/RawProcessed/new_0.5.Rda"))
