source('Forecast/Namelists/namelist_Test.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')


test <- rpss(iriForecastList,observedTercile)

fcast1 <- iriForecastList
obs    <- observedTercile
fcastBaseline <- NULL
landPropMat   <- NULL



lon  <- fcast1$lon
nlon <- length(lon)
lat  <- fcast1$lat
nlat <- length(lat)

if(is.null(landPropMat)) landPropMat <- landPropList$prop

nForecasts <- nrow(fcast1$timeMap)

# build the objects that are filled by the loop
rpssVec <- rpssTropics <- rep(NA, length=nForecasts)

# store each of the  the spatial patterns of the skill
brierForecastLow    <- array(NA, dim=c(nlon,nlat,nForecasts))
brierForecastHigh   <- array(NA, dim=c(nlon,nlat,nForecasts))
rpsForecast         <- array(NA, dim=c(nlon,nlat,nForecasts))

brierBaselineLow  <- array(NA, dim=c(nlon,nlat,nForecasts))
brierBaselineHigh <- array(NA, dim=c(nlon,nlat,nForecasts))
rpsBaseline       <- array(NA, dim=c(nlon,nlat,nForecasts))

#  lat-cosine scaling matrix to multiply before summing over the globe
cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))

# climatology forecast
climatology <- rep(1/3,3)
cumClimatology <- cumsum(climatology)


for(tInd in 1:nForecasts){
	outcomeInd <- which(obs$timeMap[,3] == fcast1$timeMap[tInd,3])
	tempForecast <- fcast1$forecast[,,tInd,]
	tempObs  <- obs$tercile[,,outcomeInd]

	if(!is.null(fcastBaseline)){
		altInd  <- which(fcastBaseline$timeMap[,3] == fcastBaseline$timeMap[tInd,3])
		tempBaseline <- fcastBaseline$forecast[,,altInd,]
	}

	for(i in 1:nlon){
	for(j in 1:nlat){

		pointForecast <- tempForecast[i,j,1:3]
		cumForecast   <- cumsum(pointForecast)

		# make sure that we have a valid forecast
		if(!is.na(pointForecast[1]) & !is.na(tempObs[i,j])){
			pointOutcome  <- rep(0,3)
			pointOutcome[tempObs[i,j] + 2 ]  <- 1

			cumOutcome <- cumsum(pointOutcome)

			rpsForecast[i,j,tInd] <- sum((cumForecast - cumOutcome)^2)

			# switch to allow for different baselines from climatology
			if(!is.null(fcastBaseline)){
				cumBaseline <- cumsum(tempBaseline[i,j,1:3])
				rpsBaseline[i,j,tInd] <- sum((cumBaseline - cumOutcome)^2)
			} else{
				rpsBaseline[i,j,tInd] <- sum((cumClimatology - cumOutcome)^2)
			}

		}
	}
	}

}


# calculate the field RPS and RPSS
rpsForecastField <- apply(rpsForecast,c(1,2),mean,na.rm=T)
rpsBaselineField <- apply(rpsBaseline,c(1,2),mean,na.rm=T)

rpssField <- 1 - rpsForecastField/rpsBaselineField

# calculate the series RPS and RPSS
tropicsLats <- which(lat < 30 & lat > -30)

weightedForecastField <- arrMult(rpsForecast,cosMat)
weightedBaselineField <- arrMult(rpsBaseline,cosMat)

rpsForecastGlobalSeries <- apply(weightedForecastField,c(3),mean,na.rm=T)
rpsBaselineGlobalSeries <- apply(weightedBaselineField,c(3),mean,na.rm=T)

rpsForecastTropicsSeries <- apply(weightedForecastField[,tropicsLats,],c(3),mean,na.rm=T)
rpsBaselineTropicsSeries <- apply(weightedBaselineField[,tropicsLats,],c(3),mean,na.rm=T)

rpssGlobalSeries  <- 1 - rpsForecastGlobalSeries/rpsBaselineGlobalSeries
rpssTropicsSeries <- 1 - rpsForecastTropicsSeries/rpsBaselineTropicsSeries

# calculate best season fields

tSeason <- fcast1$timeMap[,2]
seasonArr <- array(NA, dim=c(nrow(rpsForecast),ncol(rpsForecast),12))

for(i in 1:12){
	workingForecast <- rpsForecast[,,which(tSeason==i)]
	workingBaseline <- rpsBaseline[,,which(tSeason==i)]
 
	seasonArr[,,i] <- 1 - apply(workingForecast,c(1,2),mean,na.rm=T)/
							apply(workingBaseline,c(1,2),mean,na.rm=T)
}

suppressWarnings(bestMat <- apply(seasonArr,c(1,2),max,na.rm=T))
bestMat[bestMat==-Inf] <- NA

bestSeason <- apply(seasonArr,c(1,2),function(x) if(!all(is.na(x))) which.max(x) else NA)

# build the outlist
global  		<- list(field=rpssField,series=rpssGlobalSeries)
tropics 		<- list(series=rpssTropicsSeries)

outList <- list(lon=lon,lat=lat,timeMap=fcast1$timeMap,
				global=global,tropics=tropics,
				bestMat = bestMat, 
				bestSeason = bestSeason,
				seasonalArr = seasonArr)

return(outList)










