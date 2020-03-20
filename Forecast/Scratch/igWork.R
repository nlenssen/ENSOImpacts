# source('Forecast/Namelists/namelist_CMAP.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')


# helper functions
eirCalc <- function(x){
	2^(log2(x)-log2(1/3)) - 1
}

igCalc <- function(x,weights=NULL){
	if(is.null(weights)) weights <- rep(1,length(x))
	goodInds <- !is.na(x)
	-weighted.mean(log(x[goodInds],base=2),weights[goodInds])
}

# first create the ignorance field calculation (unweighted, no decomp)


# inputs
fcast       <- iriForecastList
obs         <- observedTercile
landPropMat <- landPropList$prop

# eventual function
lon  <- fcast$lon
nlon <- length(lon)
lat  <- fcast$lat
nlat <- length(lat)

tropicsLats <- which(lat < 30 & lat > -30)

if(is.null(landPropMat)) landPropMat <- landPropList$prop

nForecasts <- nrow(fcast$timeMap)


# match the two time series
tIndObs <- which(obs$timeMap[,3] %in% fcast$timeMap[,3])
subObs <- obs$tercile[,,tIndObs]


# stuff for area weighting
cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
weightMat <- cosMat*landPropMat


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

# calculate the EIR field
eirField <- apply(outcomeProbArr,c(1,2),function(x) mean(eirCalc(x),na.rm=T))

# calculate the series
igSeriesGlobal   <- rep(NA, nForecasts)
igSeriesTropics  <- rep(NA, nForecasts)

eirSeriesGlobal  <- rep(NA, nForecasts)
eirSeriesTropics <- rep(NA, nForecasts)

for(tInd in 1:nForecasts){
	# calculate the ignorance series
	igSeriesGlobal[tInd]  <- igCalc(c(outcomeProbArr[,,tInd]),c(weightMat))
	igSeriesTropics[tInd] <- igCalc(c(outcomeProbArr[,tropicsLats,tInd]),
								c(weightMat[,tropicsLats]))
	# calculate the EIR series
	eirSeriesGlobal[tInd]  <- weighted.mean(c(eirCalc(outcomeProbArr[,,tInd])),
		weights=c(weightMat),na.rm=T)

	eirSeriesTropics[tInd] <- weighted.mean(c(eirCalc(outcomeProbArr[,tropicsLats,tInd])),
		weights=c(weightMat[,tropicsLats]),na.rm=T)
}

# calculate the ig and EIR for the entire forecast
igTotalGlobal  <- igCalc(c(outcomeProbArr),
						weights=rep(c(weightMat),nForecasts))
igTotalTropics <- igCalc(c(outcomeProbArr[,tropicsLats,]),
						weights=rep(c(weightMat[,tropicsLats]),nForecasts))

eirTotalGlobal  <- weighted.mean(eirField,weights=weightMat,na.rm=T)
eirTotalTropics <- weighted.mean(eirField[,tropicsLats],weights=weightMat[,tropicsLats],na.rm=T)


global <- list(field  = igField, series = igSeriesGlobal)
tropics <- list(series = igSeriesTropics)
ignorance <- list(global=global, tropics = tropics)

global <- list(field  = eirField, series = eirSeriesGlobal)
tropics <- list(series = eirSeriesTropics)
eir <- list(global=global, tropics = tropics)

outList <- list(ignorance=ignorance,eir=eir)







###
# Now calculate the ignorance as the sum of two binary ignorance scores
outcomeProbArrBin <- array(NA, dim=c(nlon,nlat,nForecasts,2))
obsArrBin         <- array(NA, dim=c(nlon,nlat,nForecasts,2))
climArrBin        <- array(NA, dim=c(nlon,nlat,nForecasts,2))

clim <- rep(1/3,3)

for(cutoff in 1:2){
	tempClim <- c(sum(clim[1:cutoff]),sum(clim[(cutoff+1):3]))

	for(tInd in 1:nForecasts){
		for(i in 1:nlon){
			for(j in 1:nlat){
				tempFcastTerc <- fcast$forecast[i,j,tInd,1:3]

				tempFcast <- c(sum(tempFcastTerc[1:cutoff]),sum(tempFcastTerc[(cutoff+1):3]))

				tempObsTerc <- rep(0,3)
				tempObsTerc[subObs[i,j,tInd]+2] <- 1
				tempObs <- c(sum(tempObsTerc[1:cutoff]),sum(tempObsTerc[(cutoff+1):3]))
				# skip if no forecast or obs at a location
				if(all(is.na(tempFcast)) | is.na(subObs[i,j,tInd])){
					next
				}

				# store the probability of the event that occured
				outcomeProbArrBin[i,j,tInd,cutoff] <- tempFcast[which(tempObs==1)]
				obsArrBin[i,j,tInd,cutoff] <- tempObs[2]==1
				climArrBin[i,j,tInd,cutoff] <- tempClim[which(tempObs==1)]

			}
		}
	}
}

igFieldLower <- apply(outcomeProbArrBin[,,,1],c(1,2),igCalc)
igFieldUpper <- apply(outcomeProbArrBin[,,,2],c(1,2),igCalc)

rigField1 <- igFieldLower + igFieldUpper

rigClimField <- apply(climArrBin[,,,1],c(1,2),igCalc) + apply(climArrBin[,,,2],c(1,2),igCalc)

# qualitatively the same, by quantitatively different than the eir
eirField2 <- 2^(rigClimField - rigField1 ) - 1




# get the forecasts and outcomes into y and z notation from [Tödter and Ahrens, 2012]
yArr <- array(NA,dim=c(nlon,nlat,nForecasts,2))
zArr <- array(NA,dim=c(nlon,nlat,nForecasts,2))

for(cutoff in 1:2){
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


				tempObsTerc <- rep(0,3)
				tempObsTerc[subObs[i,j,tInd]+2] <- 1
				tempObs <- c(sum(tempObsTerc[1:cutoff]),sum(tempObsTerc[(cutoff+1):3]))
				zArr[i,j,tInd,cutoff] <- tempObs[2]
			}
		}
	}
}

# now, calculate the ignorance according to [Tödter and Ahrens, 2012]
# it validates to the prior calculation to machine precision
igFieldLower2 <- matrix(NA,nlon,nlat)
igFieldUpper2 <- matrix(NA,nlon,nlat)

for(i in 1:nlon){
	for(j in 1:nlat){
		igFieldLower2[i,j] <- - mean(log2(abs(yArr[i,j,,1]-(1-zArr[i,j,,1]))),na.rm=T)
		igFieldUpper2[i,j] <- - mean(log2(abs(yArr[i,j,,2]-(1-zArr[i,j,,2]))),na.rm=T)
	}
}

# now, we need to bin the forecasts
binSize <- 0.05

bins <- seq(0-(binSize)/2,1+(binSize)/2,by=binSize)
centers <- bins[-1]-(binSize)/2

# bins <- seq(0-binSize,1,by=binSize)
# centers <- bins[-1]-(binSize)/2

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

# bunch o test plots
image.plot(igFieldArr[,,1])	
image.plot(igFieldLower2)

image.plot(igFieldArr[,,1] - igFieldLower2)



set.panel(1,2)
zr <- range(igFieldArr[,,1],igFieldLower2,na.rm=T)
image.plot(igFieldArr[,,1],zlim=zr)
image.plot(igFieldLower2,zlim=zr)


set.panel(1,2)
zr <- range(igFieldArr[,,2],igFieldUpper2,na.rm=T)
image.plot(igFieldArr[,,2],zlim=zr)
image.plot(igFieldUpper2,zlim=zr)


# make the rigField and relTot, resTot, and uncTot fields

rigFieldMat <- apply(igFieldArr,c(1,2),sum)

relTot <- apply(relField, c(1,2), sum)
resTot <- apply(resField, c(1,2), sum)
uncTot <- apply(uncField, c(1,2), sum)

# use the ign skill score formula given as eqn (20) in Todter and Ahrens


#####
# Weighted (Series) calculation
#####
weightVec <- c(weightMat)

# first check the work by calculating the RIGN with the two cats
tropics <- FALSE

igSeriesSplit <- matrix(NA, nForecasts, 2)

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

		igSeriesSplit[tInd,cutoff] <- -weighted.mean(log2(abs(yTemp-(1-zTemp))), weights=wTemp)
	}
}

igSeriesSplitTot <- rowSums(igSeriesSplit)

# check the work
plot(fcast$timeMap[,3],igSeriesGlobal,type='l')

points(fcast$timeMap[,3],igSeriesSplitTot,type='l',col='red')


# Now try the weighted decomp
tropics <- FALSE

weightVec <- c(weightMat)

resSeries <- matrix(NA, nForecasts, 2)
relSeries <- matrix(NA, nForecasts, 2)
uncSeries <- matrix(NA, nForecasts, 2)


grid <- expand.grid(list(lon=lon,lat=lat))

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



####
# Now do the total calculation
####

tropics <- FALSE



# Comparison of the different series calculations. I would expect that all three should
# be the same. The IGN is calculated using the 
plot(fcast$timeMap[,3],igSeriesGlobal,type='l',ylim=c(1.4,2))

points(fcast$timeMap[,3],igSeries,type='l',col='red')

points(fcast$timeMap[,3],igSeriesSplitTot,type='l',col='blue')
legend('bottomright', c('RIGN Decomp', 'RIGN', 'IGN'), col=c('red','blue','black'),lty=1)


#####
# Test functions
#####
source('Forecast/Scratch/igFunctions.R')

# make the Y/Z decomp
yzTest <- yzForecast(iriForecastList,observedTercile)

# Test the field function
fieldTest <- rignField(yzTest)
image.plot(fieldTest$rign - rigField1)
image.plot(fieldTest$eir - eirField2)

# Test the series function
seriesTest <- rignSeries(yzTest,landPropMat, tropics=FALSE)

# Test the field decomp
fdTest <- rignFieldDecomp(yzTest)

# doesn't quite equal...
image.plot(fdTest$rign - fieldTest$rign)

# look at the skill patterns accord
sTemp <- fdTest$skill
sTemp[sTemp<0] <- NA
image.plot(sTemp)

# Test the series decomp
sdTest <- rignSeriesDecomp(yzTest, landPropMat, binSize=0.05, tropics=FALSE)

plot(sdTest$rign, type='l',col='red')
points(seriesTest$rign, type='l',col='blue')
## TODO

# (DONE) Get the weighted ignorance score decomposition figured out
# (DONE) Make nice functions to calculate all the things I need for the paper
# (3) Create a script to make the final figures below
# (4) Use the ignorance score for the lead analysis and make figures
# (5) Write the paper up yo

## Final Figures to make (All RANKED Ignorance Scores)

# Figure 6

# series resolution (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob

# series reliability (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob


# Figure 7

# resolution fields for
# 	- IRI
# 	- ENSO Prop

# Figure 9

# EIR Global Series (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob

# EIR Tropics Series (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob


# Figure 10

# EIR global fields for
# 	- IRI
# 	- ENSO Prop









