lead <- leads[i]


makeRealisticForecast <- function(lead, ensoForecastList,
	ninaSeasonal, resultsList, forecastList){
	# unpack the results from the Impacts study
	sampleSize <- resultsList$sampleSize
	sigCounts  <- resultsList$sigCounts
	ensoYears  <- resultsList$ensoYears
	countsArray <- resultsList$countsArray

	# time objs
	forecastInds <- which(forecastList$timeMap[,4]==lead &
		forecastList$timeMap[,1] >= startYear &
		forecastList$timeMap[,1] <= endYear &
		forecastList$timeMap[,3] %in% ensoForecastList$timeMap[,3])

	subTimeMap <- forecastList$timeMap[forecastInds,]

	eForecastInds <- which(ensoForecastList$timeMap[,3] %in% subTimeMap[,3])

	# make the nina indicator series (1 for El Nino, 2 for La Nina)
	nina34 <- ninaSeasonal$nina34
	timeMapNina <- ninaSeasonal$timeMap

	isNina <- rep(0, length(nina34))
	isNina[nina34 > ninoCutoff] <- 1
	isNina[nina34 < ninaCutoff] <- 2

	# extract the relevant info for the ENSO forecasts


	# using the count array to set deterministic forecasts
	ensoForecastArrayReal <- array(NA, dim=c(length(forecastList$lon),
		length(forecastList$lat),
		length(forecastInds),
		FCAST_CATS))

	climatology <- rep(1/3,3)
	landMask <- ifelse(is.na(sampleSize[,,1]),NA,1)


	# build a forecast for each time point

	for(tInd in 1:length(forecastInds)){
		# pull the correct Nino3.4 index
		ninaInd <- which(subTimeMap[tInd,1] == timeMapNina[,1] & 
			subTimeMap[tInd,2] == timeMapNina[,2])

		eInd <- which(subTimeMap[tInd,3] == ensoForecastList$timeMap[,3])

		
		# build objects for the loop
		tempArray <- array(NA, dim=c(length(forecastList$lon),length(forecastList$lat),3))
		tempArrayProb <- tempArray

		fcastMask <- matrix(ifelse(!(forecastList$forecast[,,forecastInds[tInd],] < 0),
			1, NA),length(forecastList$lon),length(forecastList$lat))


		##########################################
		# handle the 'realistic' sst forecast case
		##########################################

		
		# pull in the correct ENSO forecast	
		tempEnsoForecast <- ensoForecastList$fcast[eInd,]

		# get the season we are working on
		tempSeason <- forecastList$timeMap[forecastInds[tInd],2]

		# working is lon x lat x El/Nu/La x Wet/Nu/Dry
		tempArrayWorking <- array(NA, dim=c(length(forecastList$lon),length(forecastList$lat),3,3))
		
		# loop array
		tempArrayReal    <- array(0, dim=c(length(forecastList$lon),length(forecastList$lat),3))

		# multiply the probabilistic forecast by the enso forecast probability
		for(state in 1:3){
			
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





		}


		ensoRealForecastList <- list(lon=forecastList$lon,lat=forecastList$lat,
			forecast=ensoForecastArrayReal,
			timeMap=subTimeMap,timeInds = eForecastInds)

		return(ensoRealForecastList)
	}