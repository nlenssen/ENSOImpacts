# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# get forecast obj for grid info
load(file='Data/RawProcessed/forecastObjIRI.Rda')

# get seasonally processed CRU on original grid
load('Data/RawProcessed/cruSeasonal_0.5.Rda')

#unlist the necessary parts and remove rest
pptSeasonal    <- observationList$ppt
countsSeasonal <- observationList$counts
timeMap        <- observationList$timeMap
rm(observationList)

# build the empty arrays
obsArray25 <- array(NA, dim=c(length(forecastList$lon),
							  length(forecastList$lat),
							  nrow(timeMap)))

countArrayMin25 <- array(NA, dim=c(length(forecastList$lon),
							  length(forecastList$lat),
							  nrow(timeMap)))

countArrayMed25 <- array(NA, dim=c(length(forecastList$lon),
							  length(forecastList$lat),
							  nrow(timeMap)))

for(i in 1:length(forecastList$lon)){
	for(j in 1:length(forecastList$lat)){

		tempPpt <- pptSeasonal[((i-1)*5 + 1):(i*5),((j-1)*5 + 1):(j*5),]
		tempCounts <- countsSeasonal[((i-1)*5 + 1):(i*5),((j-1)*5 + 1):(j*5),]
		
		if(all(is.na(tempPpt))) next

		obsArray25[i,j,] <- apply(tempPpt,3,sum,na.rm=TRUE)

		countArrayMin25[i,j,] <- apply(tempCounts,3,min,na.rm=TRUE)

		countArrayMed25[i,j,] <- apply(tempCounts,3,median,na.rm=TRUE)	

	}
}

countArrayMin25[countArrayMin25==Inf] <- NA
countArrayMed25[countArrayMed25==Inf] <- NA

observationList <- list(lon=forecastList$lon,
						lat=forecastList$lat,
						ppt=obsArray25,
						counts=countArrayMed25,
						countsMin=countArrayMin25,
						timeMap = timeMap)

save(observationList,file='Data/RawProcessed/cruSeasonal_2.5.Rda')