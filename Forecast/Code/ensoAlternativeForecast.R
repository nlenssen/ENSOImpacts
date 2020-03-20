# load in the ENSO forecasts
load('Data/RawProcessed/ensoForecastIRI.Rda')

# load the IRI forecasts
load('Data/RawProcessed/forecastObjIRI.Rda')

# load the IRI forecasts
load('Data/RawProcessed/ninaSeasonal.Rda')

# load the count array to get spatial signal
# reminder that countArray has dimensions
# lon x lat x season x el/la/neutral x high/low anom
load(sprintf('Data/Impacts/%s/maskedAnalysisResults.Rda',impactsAnalysisName))

# load in the ENSO events as detected by the ONI method used in the Impacts
load('Data/Impacts/Cru_0.5_1951_2016_Data/ninaSeasonalInd.Rda')


# unpack the results
sampleSize <- resultsList$sampleSize
sigCounts  <- resultsList$sigCounts
ensoYears  <- resultsList$ensoYears
countsArray <- resultsList$countsArray

# time objs
forecastInds <- which(forecastList$timeMap[,4]==lead &
				      forecastList$timeMap[,1] >= startYear &
				      forecastList$timeMap[,1] <= endYear)

subTimeMap <- forecastList$timeMap[forecastInds,]


# rearrange the nino ind series to be in the format needed for this code
yearRange <- as.numeric(rownames(ninaIndicator))
timeMapNina <- cbind(rep(yearRange,each=12),rep(1:12,times=length(yearRange)),NA)
timeMapNina[,3] <- timeMapNina[,1] + (timeMapNina[,2]-1)/12

nina34 <- ninaSeasonal$nina34[which(ninaSeasonal$timeMap[,1] %in% yearRange)]

isNina <- c(t(ninaIndicator))

isNina[isNina== -1] <- 2

# optional cutoff for only strong el nino/la nina events
absoluteNinoCutoff <- 0.5
isStrongVec <- abs(nina34) >= absoluteNinoCutoff

isNina[!isStrongVec] <- 0

# extract the relevant info for the ENSO forecasts
eForecastInds <- which(subTimeMap[,3] %in% ensoForecastList$timeMap[,3])
badInds <- which(!(subTimeMap[,3] %in% ensoForecastList$timeMap[,3]))

eForecastTimeMap <- subTimeMap[-badInds,]

# using the count array to set deterministic forecasts
ensoForecastArray <- array(NA, dim=c(length(forecastList$lon),
									 length(forecastList$lat),
									 length(forecastInds),
									 FCAST_CATS))

ensoForecastArrayProb <- ensoForecastArray
ensoForecastArrayReal <- ensoForecastArray
ensoForecastArrayRealDet <- ensoForecastArray


climatology <- rep(1/3,3)
landMask <- ifelse(is.na(sampleSize[,,1]),NA,1)


# build a forecast for each time point

for(tInd in 1:length(forecastInds)){
	# pull the correct Nino3.4 index
	ninaInd <- which(subTimeMap[tInd,1] == timeMapNina[,1] & 
					 subTimeMap[tInd,2] == timeMapNina[,2])

	# determine the oracle ENSO state for the season of the forecast
	ninaState <- isNina[ninaInd]
	
	# build objects for the loop
	tempArray <- array(NA, dim=c(length(forecastList$lon),length(forecastList$lat),3))
	tempArrayProb <- tempArray

	fcastMask <- matrix(ifelse(!(forecastList$forecast[,,forecastInds[tInd],] < 0),
						 1, NA),length(forecastList$lon),length(forecastList$lat))

	##########################################
	# handle the two 'oracle-sst' forecasts
	##########################################

	if(ninaState==0){
		for(i in 1:3){
			tempArray[,,i] <- landMask * fcastMask * climatology[i]
			tempArrayProb[,,i] <- landMask * fcastMask * climatology[i]
		}
	} else{
		
		tempSeason <- forecastList$timeMap[forecastInds[tInd],2]

		highLow <- sigCounts[,,tempSeason,ninaState,]
		highLow[highLow %in% c(-1,0)] <- NA

		countsTemp <- countsArray[,,tempSeason,ninaState,]
		numYearsTemp <- ensoYears[,,tempSeason,ninaState]

		for(i in 1:length(forecastList$lon)){
			for(j in 1:length(forecastList$lat)){
				fcastVec <- highLow[i,j,]

				countsVec  <- countsTemp[i,j,]
				sampleSizeVec <- numYearsTemp[i,j]

				if(all(is.na(fcastVec))){
					tempArray[i,j,] <- climatology
					tempArrayProb[i,j,] <- climatology
				} else{
					fcast  <- rep(0,3)

					predInd <- ifelse(which(!is.na(fcastVec))==1,3,1)
					
					# handle ties in the counts
					if(length(predInd)>1){
						tempArray[i,j,] <- climatology
					} else{
						fcast[predInd] <- 1
						# check to make sure we caught all the ties
						if(sum(fcast)>1){
							print(paste('Forecast with proab greater than one!',tInd,i,j))
							break
						}
						
						tempArray[i,j,] <- fcast
					}
					
 					# deal with prop forecast by assigning the empirical
 					# counts plus 1/3 to handle zeros
 					pfcast <- rep(1/3,3)

 					midCounts <- sampleSizeVec - sum(countsVec)

 					# put in L/M/H form to aggree with the IRI forecast!!!
 					fullCounts <- c(countsVec[2], midCounts,countsVec[1])+1/3

 					tempArrayProb[i,j,] <- fullCounts/sum(fullCounts)

				}	
			}
		}

		for(k in 1:3){
			tempArray[,,k] <- tempArray[,,k] * landMask * fcastMask
			tempArrayProb[,,k] <- tempArrayProb[,,k] * landMask * fcastMask
		}
	}

	ensoForecastArray[,,tInd,] <- tempArray
	ensoForecastArrayProb[,,tInd,] <- tempArrayProb

	##########################################
	# handle the 'realistic' sst forecast case
	##########################################

	if(tInd %in% eForecastInds){
		# pull in the correct forecast	
		tempEForecastInd <- which(eForecastTimeMap[,3] == subTimeMap[tInd,3])
		tempEnsoForecast <- ensoForecastList$fcast[tempEForecastInd,]

		# working is lon x lat x El/Nu/La x Wet/Nu/Dry
		tempArrayWorking <- array(NA, dim=c(length(forecastList$lon),length(forecastList$lat),3,3))
		
		# 
		tempArrayReal    <- array(0, dim=c(length(forecastList$lon),length(forecastList$lat),3))
		tempArrayRealDet <- tempArrayReal

		# multiply the probabilistic forecast by the enso forecast probability
		for(state in 1:3){
			tempSeason <- forecastList$timeMap[forecastInds[tInd],2]

			# Handle the E/L and Neutral cases separately
			if(state != 2){
				ninaState <- ifelse(state==1,1,2)

				highLow <- sigCounts[,,tempSeason,ninaState,]
				highLow[highLow %in% c(-1,0)] <- NA

				countsTemp <- countsArray[,,tempSeason,ninaState,]
				numYearsTemp <- ensoYears[,,tempSeason,ninaState]


				for(i in 1:length(forecastList$lon)){
					for(j in 1:length(forecastList$lat)){
						fcastVec <- highLow[i,j,]

						countsVec  <- countsTemp[i,j,]
						sampleSizeVec <- numYearsTemp[i,j]

						if(all(is.na(fcastVec))){
							tempArrayWorking[i,j,state,] <- climatology
						} else{
							fcast  <- rep(0,3)

							predInd <- ifelse(which(!is.na(fcastVec))==1,3,1)
														
		 					# deal with prop forecast by assigning the empirical
		 					# counts plus 1/3 to handle zeros
		 					pfcast <- rep(1/3,3)

		 					midCounts <- sampleSizeVec - sum(countsVec)

		 					# put in L/M/H form to aggree with the IRI forecast!!!
		 					fullCounts <- c(countsVec[2], midCounts,countsVec[1])+1/3

		 					tempArrayWorking[i,j,state,] <- fullCounts/sum(fullCounts) 

						}	
					}
				}

			} else{ # If the state is neutral, we just issue a climatology forecast!
				tempArrayWorking[,,state,] <- 1/3
			}
			
		}

		# finally, take the weighted sum (by enso forecast) of the three state forecasts at
		# each tercile

		for(state in 1:3){
			for(terc in 1:3){
				tempArrayReal[,,terc] <- tempArrayReal[,,terc] + tempArrayWorking[,,state,terc] * tempEnsoForecast[state] *landMask * fcastMask
			}
		}


		ensoForecastArrayReal[,,tInd,] <- tempArrayReal

		# deal with the ensoForecastArrayRealDet

		for(i in 1:length(forecastList$lon)){
			for(j in 1:length(forecastList$lat)){
				tempReal <- tempArrayReal[i,j,]

				if(all(is.na(tempReal))) next

				tempFcast <- rep(0,3)
				fcastInd <- which.max(tempReal)
				
				# issue a det forecast if the prob is greater than detCutoff
				if(tempReal[fcastInd] > detCutoff){
					tempFcast[fcastInd] <- 1
				} else{
					tempFcast <- climatology
				}

				tempArrayRealDet[i,j,] <- tempFcast
			}
		}

		ensoForecastArrayRealDet[,,tInd,] <- tempArrayRealDet

	}


}

ensoForecastList <- list(lon=forecastList$lon,lat=forecastList$lat,
						 forecast=ensoForecastArray,
						 timeMap=subTimeMap)

ensoProbForecastList <- list(lon=forecastList$lon,lat=forecastList$lat,
						 forecast=ensoForecastArrayProb,
						 timeMap=subTimeMap)

ensoRealForecastList <- list(lon=forecastList$lon,lat=forecastList$lat,
						 forecast=ensoForecastArrayReal,
						 timeMap=subTimeMap,timeInds = eForecastInds)

ensoRealDetForecastList <- list(lon=forecastList$lon,lat=forecastList$lat,
						 forecast=ensoForecastArrayRealDet,
						 timeMap=subTimeMap,timeInds = eForecastInds)

save(ensoForecastList,ensoProbForecastList,ensoRealForecastList,ensoRealDetForecastList,
		file=sprintf('%s/ensoForecast.Rda',ddir))

