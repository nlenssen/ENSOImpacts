# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# Read the data
rawData <- fromJSON(file = sprintf('Data/Raw/%s',iriEnsoFile))[[1]]

# Season code names
seasonCode <- c('DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ',
				'JJA', 'JAS', 'ASO', 'SON', 'OND', 'NDJ')

# Build the time map
timeMapIssed <- matrix(nrow=0,ncol=2)
timeMapForecast <- matrix(NA, nrow=0,ncol=2)
forecastFull <- matrix(nrow=0,ncol=3)
colnames(forecastFull) <- c('EN', 'NU', 'LN')

nYears <- length(rawData)
tYear <- rep(NA, nYears)

for(yr in 1:nYears){
	tempList <- rawData[[yr]]
	tYear[yr] <- tempList$year

	nMonths <- length(tempList$months)
	tempMonths <- rep(NA, nMonths)

	tempForecast <- matrix(NA,nrow=nMonths,ncol=3)

	for(mn in 1:nMonths){
	tempMonths[mn] <- tempList$months[[mn]]$month

	fcastNext <- rawData[[yr]]$months[[mn]]$probabilities[[lead]]
	tempTargetSeason <- which(seasonCode==fcastNext$season)
	tempForecast[mn,]     <- c(fcastNext$elnino,fcastNext$neutral,fcastNext$lanina)/100

	tempTargetYear <- tYear[yr]
	if(tempTargetSeason==1) tempTargetYear <- tempTargetYear + 1

	timeMapForecast <- rbind(timeMapForecast,c(tempTargetYear,tempTargetSeason))
	}

	timeMapIssed <- rbind(timeMapIssed,cbind(rep(tempList$year,nMonths),tempMonths))
	forecastFull <- rbind(forecastFull,tempForecast)
}

# handle the issued forecat data
rowInds <- which(timeMapForecast[,1]<=endYear)
timeMapClip <- timeMapForecast[rowInds,]

# Formatting work
timeMapIssuedFinal<- timeMapIssed[1:nrow(timeMapClip),]
timeMapIssuedFinal[,2] <- timeMapIssuedFinal[,2]+1
timeMapIssuedFinal <- cbind(timeMapIssuedFinal,timeMapIssuedFinal[,1]+((timeMapIssuedFinal[,2]-1)/12))

# Make the final target time map
timeMap <- timeMapClip
timeMap <- cbind(timeMap,timeMap[,1] + (timeMap[,2]-1)/12)

forecast <- forecastFull[rowInds,]

ensoForecastList <- list(fcast = forecast, timeMap=timeMap, timeMapIssued = timeMapIssuedFinal)

save(ensoForecastList,file='Data/RawProcessed/ensoForecastIRI.Rda')