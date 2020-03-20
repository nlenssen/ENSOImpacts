# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# get forecast obj for grid info
load(file='Data/RawProcessed/forecastObjIRI.Rda')

# get seasonally processed CRU on original grid
load('Data/RawProcessed/landProp_0.5.Rda')

#unlist the necessary parts
propMat  <- landPropList$prop

#
propMat25 <- matrix(NA, nrow=length(forecastList$lon),
						ncol=length(forecastList$lat))


for(i in 1:length(forecastList$lon)){
	for(j in 1:length(forecastList$lat)){

		tempProp <- propMat[((i-1)*5 + 1):(i*5),((j-1)*5 + 1):(j*5)]
		
		if(all(is.na(tempProp))) next

		propMat25[i,j] <- mean(tempProp,na.rm=TRUE)
	}
}

landPropList <- list(lon=forecastList$lon,
					 lat=forecastList$lat,
					 prop=propMat25)

save(landPropList,file='Data/RawProcessed/landProp_2.5.Rda')