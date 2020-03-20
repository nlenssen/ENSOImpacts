# get the forecasts and outcomes into y and z notation from [Tödter and Ahrens, 2012]
yzForecast <- function(fcast,obs){
	lon  <- fcast$lon
	nlon <- length(lon)
	lat  <- fcast$lat
	nlat <- length(lat)

	nForecasts <- nrow(fcast$timeMap)

	# match the two time series
	tIndObs <- which(obs$timeMap[,3] %in% fcast$timeMap[,3])
	subObs <- obs$tercile[,,tIndObs]

	# 
	yArr <- array(NA,dim=c(nlon,nlat,nForecasts,2))
	zArr <- array(NA,dim=c(nlon,nlat,nForecasts,2))
	climArr <- array(NA,dim=c(nlon,nlat,nForecasts,2))

	clim <- rep(1/3,3)

	# make a grid list
	grid <- expand.grid(list(lon=lon,lat=lat))

	for(cutoff in 1:2){
		tempClim <- c(sum(clim[1:cutoff]),sum(clim[(cutoff+1):3]))

		for(tInd in 1:nForecasts){
			for(i in 1:nlon){
				for(j in 1:nlat){
					tempFcastTerc <- fcast$forecast[i,j,tInd,1:3]

					tempFcast <- c(sum(tempFcastTerc[1:cutoff]),sum(tempFcastTerc[(cutoff+1):3]))

					# skip if no forecast or obs at a location
					if(all(is.na(tempFcast)) | is.na(subObs[i,j,tInd])){
						next
					}

					yArr[i,j,tInd,cutoff] <- tempFcast[2]
					climArr[i,j,tInd,cutoff] <- tempClim[2]

					tempObsTerc <- rep(0,3)
					tempObsTerc[subObs[i,j,tInd]+2] <- 1
					tempObs <- c(sum(tempObsTerc[1:cutoff]),sum(tempObsTerc[(cutoff+1):3]))
					zArr[i,j,tInd,cutoff] <- tempObs[2]
				}
			}
		}
	}

	return(list(yArr = yArr, zArr = zArr, climArr = climArr, 
				lon = lon, lat = lat, nForecasts = nForecasts,
				grid = grid))
}


# Function to calculate the RIGN (no decomp) for a field
rignField <- function(yzList){
	# unpack the yzList
	yArr <- yzList$yArr
	zArr <- yzList$zArr
	climArr <- yzList$climArr
	lon <- yzList$lon
	lat <- yzList$lat
	nForecasts <- yzList$nForecasts
	grid <- yzList$grid
	
	nlon <- length(lon)
	nlat <- length(lat)

	# build the empty objects
	igFieldArr         <- array(NA, dim=c(nlon,nlat,2))
	climFieldArr        <- array(NA, dim=c(nlon,nlat,2))

	igFieldArr <- array(NA,dim=c(nlon,nlat,2))

	# calculate the IGN scores for each location/cutoff/forecast
	for(cutoff in 1:2){
		for(i in 1:nlon){
			for(j in 1:nlat){
				igFieldArr[i,j,cutoff] <- - mean(log2(abs(yArr[i,j,,cutoff]-(1-zArr[i,j,,cutoff]))),na.rm=T)
				climFieldArr[i,j,cutoff] <- - mean(log2(abs(climArr[i,j,,cutoff]-(1-zArr[i,j,,cutoff]))),na.rm=T)
			}
		}
	}

	# calculate the RIGN fields for the two forecasts
	field <- apply(igFieldArr, c(1,2), sum)
	climField <- apply(climFieldArr, c(1,2), sum)

	# qualitatively the same, by quantitatively different than the eir
	eirField <- 2^(climField - field ) - 1

	outList <- list(rign = field, eir = eirField, lon=lon, lat=lat)

	return(outList)
}


# calculate a RIGN series
rignSeries <- function(yzList, landPropMat, tropics = FALSE){

	# unpack the yzList
	yArr <- yzList$yArr
	zArr <- yzList$zArr
	climArr <- yzList$climArr
	lon <- yzList$lon
	lat <- yzList$lat
	nForecasts <- yzList$nForecasts
	grid <- yzList$grid

	nlon <- length(lon)
	nlat <- length(lat)

	# build a vector of weights
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
	weightMat <- cosMat*landPropMat
	weightVec <- c(weightMat)

	# build the empty arrays to store results
	igSeriesSplit <- matrix(NA, nForecasts, 2)
	climSeriesSplit <- matrix(NA, nForecasts, 2)

	for(cutoff in 1:2){
		for(tInd in 1:nForecasts){
			yTempRaw <- c(yArr[,,tInd,cutoff])
			zTempRaw <- c(zArr[,,tInd,cutoff])
			climTempRaw <- c(climArr[,,tInd,cutoff])

			domainBadInds <- NULL
			if(tropics){
				domainBadInds <- which(abs(grid[,2])>30)
			}

			yBadInds <- which(is.na(yTempRaw))
			zBadInds <- which(is.na(zTempRaw))

			badInds <- union(yBadInds,union(zBadInds,domainBadInds))

			if(length(badInds)>0){
				yTemp <- yTempRaw[-badInds]
				zTemp <- zTempRaw[-badInds]
				wTemp <- weightVec[-badInds]
				climTemp <- climTempRaw[-badInds]
			} else{
				yTemp <- yTempRaw
				zTemp <- zTempRaw
				wTemp <- weightVec
				climTemp <- climTempRaw
			}

			igSeriesSplit[tInd,cutoff] <- -weighted.mean(log2(abs(yTemp-(1-zTemp))), weights=wTemp)
			climSeriesSplit[tInd,cutoff] <- -weighted.mean(log2(abs(climTemp - (1-zTemp))), weight=wTemp)
		}
	}

	rigSeries <- rowSums(igSeriesSplit)
	climSeries <- rowSums(climSeriesSplit)

	# calculate the eir
	eirSeries <- 2^(climSeries - rigSeries) - 1

	outList <- list(rign = rigSeries, clim = climSeries, eir = eirSeries)

	return(outList)
}

# calculate the total rign
rignTotal <- function(yzList, landPropMat, tropics = FALSE){

	# unpack the yzList
	yArr <- yzList$yArr
	zArr <- yzList$zArr
	climArr <- yzList$climArr
	lon <- yzList$lon
	lat <- yzList$lat
	nForecasts <- yzList$nForecasts
	grid <- yzList$grid

	nlon <- length(lon)
	nlat <- length(lat)

	# build a vector of weights
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
	weightMat <- cosMat*landPropMat


	# build a larger vector of weights to line up
	weightVec <- rep(c(weightMat),nForecasts)

	# build the empty arrays to store results
	igPart   <- rep(NA, 2)
	climPart <- rep(NA, 2)

	for(cutoff in 1:2){
		yTempRaw <- c(yArr[,,,cutoff])
		zTempRaw <- c(zArr[,,,cutoff])
		climTempRaw <- c(climArr[,,,cutoff])
		
		domainBadInds <- NULL
		if(tropics){
			domainBadInds <- which(rep(abs(grid[,2]),nForecasts)>30)
		}

		yBadInds <- which(is.na(yTempRaw))
		zBadInds <- which(is.na(zTempRaw))

		badInds <- union(yBadInds,union(zBadInds,domainBadInds))

		if(length(badInds)>0){
			yTemp <- yTempRaw[-badInds]
			zTemp <- zTempRaw[-badInds]
			wTemp <- weightVec[-badInds]
			climTemp <- climTempRaw[-badInds]

		} else{
			yTemp <- yTempRaw
			zTemp <- zTempRaw
			wTemp <- weightVec
			climTemp <- climTempRaw

		}


		igPart[cutoff] <- -weighted.mean(log2(abs(yTemp-(1-zTemp))), weights=wTemp)
		climPart[cutoff] <- -weighted.mean(log2(abs(climTemp - (1-zTemp))), weight=wTemp)
	}

	rig <- sum(igPart)
	clim <- sum(climPart)

	eir <- 2^(clim - rig) - 1

	outList <- list(rign = rig, clim = clim, eir = eir)

	return(outList)
}


# field decomp according to [Tödter and Ahrens, 2012]
rignFieldDecomp <- function(yzList, binSize=0.05){

	# unpack the yzList
	yArr <- yzList$yArr
	zArr <- yzList$zArr
	climArr <- yzList$climArr
	lon <- yzList$lon
	lat <- yzList$lat
	nForecasts <- yzList$nForecasts
	grid <- yzList$grid


	nlon <- length(lon)
	nlat <- length(lat)


	# create the forecast bins
	bins <- seq(0-(binSize)/2, 1+(binSize)/2, by=binSize)
	centers <- bins[-1]-(binSize)/2

	# create the empty arrays
	igFieldArr <- array(NA, dim=c(nlon,nlat,2))

	relField <- array(NA, dim=c(nlon,nlat,2))
	resField <- array(NA, dim=c(nlon,nlat,2))
	uncField <- array(NA, dim=c(nlon,nlat,2))

	for(i in 1:nlon){
		for(j in 1:nlat){
			for(cutoff in 1:2){
				# make sure we have no NA values
				yTempRaw <- yArr[i,j,,cutoff]
				yBadInds <- which(is.na(yTempRaw))

				zTempRaw <- zArr[i,j,,cutoff]
				zBadInds <- which(is.na(zTempRaw))

				badInds <- union(yBadInds,zBadInds)

				if(length(badInds)>0){
					yTemp <- yTempRaw[-badInds]
					zTemp <- zTempRaw[-badInds]
				} else{
					yTemp <- yTempRaw
					zTemp <- zTempRaw
				}
				
				# calculate the pieces needed for the decomposition
				pyi <- rep(NA, length(centers))
				zBari <- rep(NA, length(centers))

				for(b in 1:length(centers)){
					# get the inds where the forecast falls into the bin
					subInds <- which(yTemp > bins[b] & yTemp <= bins[b+1])

					pyi[b] <- length(subInds)/length(yTemp)

					if(length(subInds)>0){
						zBari[b] <- mean(zTemp[subInds])
					} else{
						zBari[b] <- 0
					}
				}

				zBar <- mean(zTemp)

				# now, calculate the rel, res, and unc for the grid box
				relVec <- pyi * (zBari*log2(zBari/centers) + (1-zBari)*log2((1-zBari)/(1-centers)))
				rel <- sum(relVec[relVec<Inf],na.rm=T)
				
				resVec <- pyi * (zBari*log2(zBari/zBar) + (1-zBari)*log2((1-zBari)/(1-zBar)))
				res <- sum(resVec[resVec<Inf],na.rm=T)

				unc <- -zBar * log2(zBar) - (1-zBar)*log2(1-zBar)

				relField[i,j,cutoff] <- rel
				resField[i,j,cutoff] <- res
				uncField[i,j,cutoff] <- unc

				igFieldArr[i,j,cutoff] <- rel - res + unc
			}
		}
	}

	# sum the values of the cutoffs
	rigFieldMat <- apply(igFieldArr,c(1,2),sum)

	relTot <- apply(relField, c(1,2), sum)
	resTot <- apply(resField, c(1,2), sum)
	uncTot <- apply(uncField, c(1,2), sum)

	# skill calc from eqn (20) of T+A 2012
	rignSS <- (resTot - relTot)/uncTot

	# build the output list
	outList <- list(rign = rigFieldMat, rel=relTot, res = resTot, unc = uncTot, skill = rignSS)

	return(outList)
}


rignSeriesDecomp <- function(yzList, landPropMat, binSize=0.05, tropics=FALSE){

	# unpack the yzList
	yArr <- yzList$yArr
	zArr <- yzList$zArr
	climArr <- yzList$climArr
	lon <- yzList$lon
	lat <- yzList$lat
	nForecasts <- yzList$nForecasts
	grid <- yzList$grid


	nlon <- length(lon)
	nlat <- length(lat)

	# build a vector of weights
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
	weightMat <- cosMat*landPropMat
	weightVec <- c(weightMat)

	# create the forecast bins
	bins <- seq(0-(binSize)/2, 1+(binSize)/2, by=binSize)
	centers <- bins[-1]-(binSize)/2

	# create the empty objects
	resSeries <- matrix(NA, nForecasts, 2)
	relSeries <- matrix(NA, nForecasts, 2)
	uncSeries <- matrix(NA, nForecasts, 2)

	# loop through cutoffs and time points
	for(cutoff in 1:2){
		for(tInd in 1:nForecasts){

			yTempRaw <- c(yArr[,,tInd,cutoff])
			zTempRaw <- c(zArr[,,tInd,cutoff])

			domainBadInds <- NULL
			if(tropics){
				domainBadInds <- which(abs(grid[,2])>30)
			}

			yBadInds <- which(is.na(yTempRaw))
			zBadInds <- which(is.na(zTempRaw))

			badInds <- union(yBadInds,union(zBadInds,domainBadInds))

			if(length(badInds)>0){
				yTemp <- yTempRaw[-badInds]
				zTemp <- zTempRaw[-badInds]
				wTemp <- weightVec[-badInds]
			} else{
				yTemp <- yTempRaw
				zTemp <- zTempRaw
				wTemp <- weightVec
			}

			# calculate the zBar

			zBar <- weighted.mean(zTemp,w=wTemp)

			pyi <- rep(NA, length(centers))
			zBari <- rep(NA, length(centers))

			for(b in 1:length(centers)){
				# get the inds where the forecast falls into the bin
				subInds <- which(yTemp > bins[b] & yTemp <= bins[b+1])

				pyi[b] <- sum(wTemp[subInds])/sum(wTemp)

				if(length(subInds)>0){
					zBari[b] <- weighted.mean(zTemp[subInds],w=wTemp[subInds])
				} else{
					zBari[b] <- 0
				}
			}


			relVec <- pyi * (zBari*log2(zBari/centers) + (1-zBari)*log2((1-zBari)/(1-centers)))
			relSeries[tInd,cutoff] <- sum(relVec[relVec<Inf],na.rm=T)

			resVec <- pyi * (zBari*log2(zBari/zBar) + (1-zBari)*log2((1-zBari)/(1-zBar)))
			resSeries[tInd,cutoff] <- sum(resVec[resVec<Inf],na.rm=T)

			uncSeries[tInd,cutoff] <- -zBar * log2(zBar) - (1-zBar)*log2(1-zBar)
		}
	}

	resSeriesTot <- rowSums(resSeries)
	relSeriesTot <- rowSums(relSeries)
	uncSeriesTot <- rowSums(uncSeries)

	igSeries <- relSeriesTot - resSeriesTot + uncSeriesTot

	outList <- list(rign=igSeries, rel=relSeriesTot, res=resSeriesTot, unc=uncSeriesTot)
	return(outList)
}

rignTotalDecomp <- function(yzList, landPropMat, binSize=0.05, tropics=FALSE){
	# unpack the yzList
	yArr <- yzList$yArr
	zArr <- yzList$zArr
	climArr <- yzList$climArr
	lon <- yzList$lon
	lat <- yzList$lat
	nForecasts <- yzList$nForecasts
	grid <- yzList$grid


	nlon <- length(lon)
	nlat <- length(lat)

	# build a matrix of weights
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
	weightMat <- cosMat*landPropMat

	# create the forecast bins
	bins <- seq(0-(binSize)/2, 1+(binSize)/2, by=binSize)
	centers <- bins[-1]-(binSize)/2

	# build a larger vector of weights to line up
	weightVec <- rep(c(weightMat),nForecasts)

	# store output here
	resPart <- rep(NA, 2)
	relPart <- rep(NA, 2)
	uncPart <- rep(NA, 2)

	for(cutoff in 1:2){
		yTempRaw <- c(yArr[,,,cutoff])
		zTempRaw <- c(zArr[,,,cutoff])

		domainBadInds <- NULL
		if(tropics){
			domainBadInds <- which(rep(abs(grid[,2]),nForecasts)>30)
		}

		yBadInds <- which(is.na(yTempRaw))
		zBadInds <- which(is.na(zTempRaw))

		badInds <- union(yBadInds,union(zBadInds,domainBadInds))

		if(length(badInds)>0){
			yTemp <- yTempRaw[-badInds]
			zTemp <- zTempRaw[-badInds]
			wTemp <- weightVec[-badInds]
		} else{
			yTemp <- yTempRaw
			zTemp <- zTempRaw
			wTemp <- weightVec
		}

		# calculate the zBar
		zBar <- weighted.mean(zTemp,w=wTemp)

		# calculate the other pieces
		pyi <- rep(NA, length(centers))
		zBari <- rep(NA, length(centers))

		for(b in 1:length(centers)){
			# get the inds where the forecast falls into the bin
			subInds <- which(yTemp > bins[b] & yTemp <= bins[b+1])

			pyi[b] <- sum(wTemp[subInds])/sum(wTemp)

			if(length(subInds)>0){
				zBari[b] <- weighted.mean(zTemp[subInds],w=wTemp[subInds])
			} else{
				zBari[b] <- 0
			}
		}


		relVec <- pyi * (zBari*log2(zBari/centers) + (1-zBari)*log2((1-zBari)/(1-centers)))
		relPart[cutoff] <- sum(relVec[relVec<Inf],na.rm=T)

		resVec <- pyi * (zBari*log2(zBari/zBar) + (1-zBari)*log2((1-zBari)/(1-zBar)))
		resPart[cutoff] <- sum(resVec[resVec<Inf],na.rm=T)

		uncPart[cutoff] <- -zBar * log2(zBar) - (1-zBar)*log2(1-zBar)
	}

	resTot <- sum(resPart)
	relTot <- sum(relPart)
	uncTot <- sum(uncPart)

	rig <- relTot - resTot + uncTot

	outList <- list(rig = rig, res = resTot, rel = relTot, unc = uncTot)

	return(outList)
}


# updated/corrected all-in-one EIR calculation
eirCalcFull <- function(fcast,obs,landPropMat=NULL){
	# eventual function
	lon  <- fcast$lon
	nlon <- length(lon)
	lat  <- fcast$lat
	nlat <- length(lat)

	tropicsLats <- which(lat < 30 & lat > -30)

	# stuff for area weighting
	if(is.null(landPropMat)) landPropMat <- landPropList$prop
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
	weightMat <- cosMat*landPropMat


	nForecasts <- nrow(fcast$timeMap)


	# match the two time series
	tIndObs <- which(obs$timeMap[,3] %in% fcast$timeMap[,3])
	subObs <- obs$tercile[,,tIndObs]


	# build the array of event probabilities
	outcomeProbArr <- array(NA, dim=c(nlon,nlat,nForecasts))

	for(tInd in 1:nForecasts){
		for(i in 1:nlon){
			for(j in 1:nlat){
				tempFcast <- fcast$forecast[i,j,tInd,1:3]
				tempObs <- subObs[i,j,tInd]
			
				# skip if no forecast or obs at a location
				if(all(is.na(tempFcast)) | all(is.na(tempObs))){
					next
				}

				# store the probability of the event that occured
				outcomeProbArr[i,j,tInd] <- tempFcast[tempObs+2]
			}
		}
	}

	# calculate the field
	igField <- apply(outcomeProbArr,c(1,2),igCalc)

	# calculate the climatology ignorance field
	climArr <- ifelse(is.na(outcomeProbArr),NA, 1/3)
	climField <- apply(climArr,c(1,2),igCalc)

	# calculate the eir field
	eirField <- 2^(climField - igField) - 1

	######################
	# calculate the series
	######################

	eirSeriesGlobal  <- rep(NA, nForecasts)
	eirSeriesTropics <- rep(NA, nForecasts)

	for(tInd in 1:nForecasts){

		igMat <- - log2(outcomeProbArr[,,tInd])

		climMat <- ifelse(is.na(outcomeProbArr[,,tInd]),NA,1/3)
		igClimMat <- - log2(climMat)

		eirSeriesGlobal[tInd] <- weighted.mean(2^(c(igClimMat) - c(igMat)) - 1,
			weights=c(weightMat),na.rm=T)

		eirSeriesTropics[tInd] <- weighted.mean(2^(c(igClimMat[,tropicsLats]) - c(igMat[,tropicsLats])) - 1,
			weights=c(weightMat[,tropicsLats]),na.rm=T)
	}


	#####################
	# calculate the total
	#####################

	# first, we need to get a count of the number of forecasts at each location
	countPropMat <- apply(outcomeProbArr,c(1,2), function(x) sum(!is.na(x)))/dim(outcomeProbArr)[3]

	#create a new weight matrix based on forecast quantity and area
	weightPropMat <- countPropMat * weightMat

	# now, take the weighted mean of the EIR field for the total values
	eirTotalGlobal <- weighted.mean(eirField,weights=weightPropMat,na.rm=T)
	eirTotalTropics <- weighted.mean(eirField[,tropicsLats],weights=weightPropMat[,tropicsLats],na.rm=T)

	outList <- list(field = eirField, globalSeries = eirSeriesGlobal, tropicsSeries = eirSeriesTropics,
					globalTotal = eirTotalGlobal, tropicsTotal = eirTotalTropics)

	return(outList)
}


# helper funciton to calculate weighted ignorance
igCalc <- function(x,weights=NULL){
	if(is.null(weights)) weights <- rep(1,length(x))
	goodInds <- !is.na(x)
	-weighted.mean(log(x[goodInds],base=2),weights[goodInds])
}

