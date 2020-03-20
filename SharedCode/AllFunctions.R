#  add the nino indicator rectangles (need ninaIndicator loaded)
addNino <- function(y){
	yearRange <- as.numeric(rownames(ninaIndicator))
	timeMapNina <- cbind(rep(yearRange,each=12),rep(1:12,times=length(yearRange)),NA)
	timeMapNina[,3] <- timeMapNina[,1] + (timeMapNina[,2]-1)/12
	isNina <- c(t(ninaIndicator))

	xNino34 <- ninaSeasonal$timeMap[,3]


	activeNinoInds <- which(isNina == 1)
	activeNinaInds <- which(isNina == -1)

	ninoRec <- matrix(NA, nrow=length(activeNinoInds),ncol=4)
	ninaRec <- matrix(NA, nrow=length(activeNinaInds),ncol=4)

	for(i in 1:length(activeNinoInds)){
		ind <- activeNinoInds[i]
		ninoRec[i,] <- c(rep(xNino34[ind],2),rep(xNino34[ind+1],2))
	}

	for(i in 1:length(activeNinaInds)){
		ind <- activeNinaInds[i]
		ninaRec[i,] <- c(rep(xNino34[ind],2),rep(xNino34[ind+1],2))
	}

	colElNino <- adjustcolor('red', alpha=0.2)
  	colLaNina <- adjustcolor('blue', alpha=0.2)

	# add polygons
	for(i in 1:nrow(ninoRec)){
		polygon(x=ninoRec[i,],y=y,col=colElNino,border=NA)
	}

	for(i in 1:nrow(ninaRec)){
		polygon(x=ninaRec[i,],y=y,col=colLaNina,border=NA)
	}
}


# element-multiply each field in an array by a matrix (INEFFICIENT)
arrMult <- function(arr,mat){
	outArr <- array(NA, dim=dim(arr))
	for(i in 1:dim(arr)[3]){
		outArr[,,i] <- arr[,,i] * mat
	}
	return(outArr)
}

# get the best season of a (non-GROC) skill statistic
bestSeason <- function(arr,tSeason){
	seasonArr <- array(NA, dim=c(nrow(arr),ncol(arr),12))

	for(i in 1:12){
		tempSlice <- arr[,,which(tSeason==i)]
		seasonArr[,,i] <- apply(tempSlice,c(1,2),mean,na.rm=T)
	}

	suppressWarnings(bestMat <- apply(seasonArr,c(1,2),max,na.rm=T))
	bestMat[bestMat==-Inf] <- NA

	bestSeason <- apply(seasonArr,c(1,2),function(x) if(!all(is.na(x))) which.max(x) else NA)

	outList <- list(bestMat = bestMat, bestSeason = bestSeason, seasonalArr = seasonArr)
}

# Brier decomposition according to Murphy 1973 following the notation of
# Stephenson et al. 2008.
brierDecomp73 <- function(forecast,outcome){
	n <- nrow(forecast)

	# to make sure unique doesn't give us something funny
	forecast <- round(forecast,4)

	nK <- table(forecast[,2])
	names(nK) <- NULL
	forecastProb <- sort(unique(forecast[,2]))

	oBarK <- rep(NA, length(nK))

	for(k in 1:length(oBarK)){
		oBarK[k] <- mean(outcome[forecast[,2]==forecastProb[k]])
	}

	oBar <- 1/n*sum(nK*oBarK)

	rel <- 1/n * sum(nK*(forecastProb - oBarK)^2)
	res <- 1/n * sum(nK*(oBarK-oBar)^2)
	unc <- oBar * (1-oBar)

	brier <- rel - res + unc

	outVec <- c(brier,rel,res,unc)
	names(outVec) <- c('Brier', 'REL', 'RES', 'UNC')

	return(outVec)
}

# Weighted brier decomp following the notation of [Young 2010]
brierDecomp73Weighted <- function(forecast,outcome,weights){
	W <- sum(weights)

	# to make sure unique doesn't give us something funny
	forecast <- round(forecast,4)

	forecastProb <- sort(unique(forecast[,2]))

	wt <- rep(NA, length(forecastProb))
	dBart <- rep(NA, length(forecastProb))

	# get weighted relative observed occurence frequecies
	# get weighted occurences of probability cats
	for(k in 1:length(dBart)){
		subInds <- which(forecast[,2]==forecastProb[k])
		wt[k] <- sum(weights[subInds])
		dBart[k] <- 1/wt[k]*sum(outcome[subInds]*weights[subInds])
	}

	dBar <- 1/W*sum(weights*outcome)

	rel <- 1/W * sum(wt*(forecastProb - dBart)^2)
	res <- 1/W * sum(wt*(dBart-dBar)^2)
	unc <- dBar * (1-dBar)

	brier <- rel - res + unc

	outVec <- c(brier,rel,res,unc)
	names(outVec) <- c('Brier', 'REL', 'RES', 'UNC')

	return(outVec)
}

brierDecomp10WeightedBinned <- function(forecast,outcome,weights,binSize=0.05){
	W <- sum(weights)

	# make bins
	bins <- seq(0,1,by=binSize)
	centers <- bins[-1]-(binSize)/2

	# determine the bin counts
	nk <- rep(NA, length(centers))
	fkBar <- rep(NA,length(centers))
	wt <- rep(NA, length(centers))
	dBart <- rep(NA, length(centers))

	# calculate dBar (or oBar in 2010 notation)
	dBar <- 1/W*sum(weights*outcome)

	# Do first bin outside loop to include zero
	tempInds <- which(forecast[,2]>=bins[1] & forecast[,2]<=bins[2])
	nk[1] <- length(tempInds)
	fkBar[1] <- mean(forecast[tempInds,2])
	wt[1] <- sum(weights[tempInds])
	dBart[1] <- 1/wt[1]*sum(outcome[tempInds]*weights[tempInds])

	# loop through the rest of the bins
	for(i in 2:length(nk)){
		tempInds <- which(forecast[,2]>bins[i] & forecast[,2]<=bins[i+1])
		nk[i] <- length(tempInds)
		fkBar[i] <- mean(forecast[tempInds,2])
		wt[i] <- sum(weights[tempInds])
		dBart[i] <- 1/wt[i]*sum(outcome[tempInds]*weights[tempInds])
	}


	# calculate the traditional three components
	rel <- 1/W * sum(wt*(fkBar - dBart)^2,na.rm=T)
	res <- 1/W * sum(wt*(dBart-dBar)^2,na.rm=T)
	unc <- dBar * (1-dBar)

	# calculate the two bin error components accourding to Stephenson 2010
	brier <- rel - res + unc

	outVec <- c(brier,rel,res,unc)
	names(outVec) <- c('Brier', 'REL', 'RES', 'UNC')

	return(outVec)
}

# Brier decomposition according to Murphy 1986
brierDecomp86 <- function(forecast,outcome){
	r1 <- forecast[outcome==1,2]
	r0 <- forecast[outcome==0,2]

	dBar1 <- length(r1)/length(outcome)
	dBar0 <- length(r0)/length(outcome)

	r1Var <- dBar1*var(r1)
	if(is.na(r1Var)) r1Var <- 0

	r0Var <- dBar0*var(r0)
	if(is.na(r0Var)) r0Var <- 0

	varR <-  r1Var + r0Var

	r1Bias <- dBar1*(mean(r1)-1)^2
	if(is.na(r1Bias)) r1Bias <- 0
	r0Bias <- dBar0*(mean(r0)-0)^2
	if(is.na(r0Bias)) r0Bias <- 0
	biasR <- r1Bias + r0Bias

	brier <- varR + biasR

	outVec <- c(brier,varR,biasR)
	names(outVec) <- c('brier', 'var', 'bias')

	return(outVec)
}

# Weighted brier following the notation of [Young 2010]
brierWeighted <- function(forecast,outcome,weights){
	W <- sum(weights)
	brier <- 1/W * sum(weights*(forecast[,2]-outcome)^2)
	return(brier)
}

# function to calculate all of the Effective Interest Rate output
eir <- function(fcast,obs,fcastBaseline=NULL,landPropMat=NULL){
	lon  <- fcast$lon
	nlon <- length(lon)
	lat  <- fcast$lat
	nlat <- length(lat)

	if(is.null(landPropMat)) landPropMat <- landPropList$prop

	tropicsLats <- which(lat < 30 & lat > -30)


	nForecasts <- nrow(fcast$timeMap)


	# match the two time series
	tIndObs <- which(obs$timeMap[,3] %in% fcast$timeMap[,3])
	subObs <- obs$tercile[,,tIndObs]


	outcomeProbArr       <- array(NA, dim=c(nlon,nlat,nForecasts))
	baselineScoreArr <- array(NA, dim=c(nlon,nlat,nForecasts))

	for(tInd in 1:nForecasts){
		for(i in 1:nlon){
			for(j in 1:nlat){
				tempFcast <- fcast$forecast[i,j,tInd,1:3]
				tempObs <- subObs[i,j,tInd]

				if(!is.null(fcastBaseline)){
					tempBaseline <- fcastBaseline$forecast[i,j,tInd,1:3]
				} else{
					tempBaseline <- NULL
				}
			

				if(all(is.na(tempFcast)) | all(is.na(tempObs)) | 
					(!is.null(tempBaseline) & all(is.na(tempBaseline)))){
					next
				}

				outcomeProbArr[i,j,tInd] <- tempFcast[tempObs+2]
				
				# If we are using climatology, the outcome forecast
				# is always just 1/3, othwise use the baseline forecast
				if(is.null(fcastBaseline)){
					baselineScoreArr[i,j,tInd] <- 1/3	
				} else{
					baselineScoreArr[i,j,tInd] <- tempBaseline[tempObs+2]
				}
				
			}
		}
	}

	# we want to calculate four things
	# (1) global time series
	# (2) tropics time series
	# (3) global average EIR
	# (4) best season EIR

	# stuff for area weighting
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))
	weightMat <- cosMat*landPropMat

	# (1)/(2) Make the series calculations
	globalIgSeries  <- rep(NA, nForecasts)
	tropicsIgSeries <- rep(NA, nForecasts)

	globalBaselineSeries  <- rep(NA, nForecasts)
	tropicsBaselineSeries <- rep(NA, nForecasts)
	for(tInd in 1:nForecasts){
		# need to area and land area weight
		globalIgSeries[tInd]  <- log2Mean(c(outcomeProbArr[,,tInd]),c(weightMat))
		tropicsIgSeries[tInd] <- log2Mean(c(outcomeProbArr[,tropicsLats,tInd]),
									c(weightMat[,tropicsLats]))

		globalBaselineSeries[tInd]  <- log2Mean(c(baselineScoreArr[,,tInd]),c(weightMat))
		tropicsBaselineSeries[tInd] <- log2Mean(c(baselineScoreArr[,tropicsLats,tInd]),
									c(weightMat[,tropicsLats]))
	}

	ierSeriesGlobal  <- eirCalc(globalIgSeries,globalBaselineSeries)
	ierSeriesTropics <- eirCalc(tropicsIgSeries,tropicsBaselineSeries)

	# (3) Global average calculation

	igScore       <- apply(outcomeProbArr,c(1,2),log2Mean)
	baselineScore <- apply(baselineScoreArr,c(1,2),log2Mean)

	ierAnnualField <- eirCalc(igScore,baselineScore)

	# (4) best season field

	ierBestSeason <- ifelse(!is.na(ierAnnualField),-1,NA)
	whichSeason <- matrix(NA, nrow=nlon,ncol=nlat)
	for(s in 1:12){
		seasonInds <- which(fcast$timeMap[,2]==s)
		igTemp <- apply(outcomeProbArr[,,seasonInds],c(1,2),log2Mean)
		baselineTemp <- apply(baselineScoreArr[,,seasonInds],c(1,2),log2Mean)

		ierTemp <- eirCalc(igTemp,baselineTemp)

		# inefficient!
		for(i in 1:nlon){
			for(j in 1:nlat){
				if(!is.na(ierTemp[i,j]) & ierTemp[i,j]>ierBestSeason[i,j]){
					ierBestSeason[i,j] <- ierTemp[i,j]
					whichSeason[i,j] <- s
				}
			}
		}
	}

	global <- list(averageField = ierAnnualField,
				   bestField    = ierBestSeason,
				   whichSeason  = whichSeason,
				   series       = ierSeriesGlobal)

	tropics <- list(series = ierSeriesTropics)

	outList <- list(global=global,tropics=tropics)

	return(outList)
}

eirCalc <- function(ig,baseline){
	2^(baseline-ig) - 1
}

# Plot of the EIR field
eirFieldPlot <- function(lon,lat,z,plotdir,plotName,...){

	bins <- c(-1,0,0.01,0.025,0.05,0.1,0.2,0.3,1)
	pal  <- c('gray40', tim.colors(length(bins)-2))
	
	pdf(sprintf('%s/%s.pdf',plotdir,plotName),10,6)
	image(lon,lat,z,breaks=bins,col=pal,ylim=c(-60,90),...)
	world(add=T)
	dev.off()
}

# ajust a field to a given FDR, option for visualization plot
fdrFunction <- function(tempField, FDR = 0.25, plotName = NULL){
    goodTests <- which(!is.na(c(tempField)))

    dataVec <- c(1-tempField)[goodTests]

    iRank <- rank(dataVec)
    iOrder <- order(dataVec)

    mTotalTests <- length(goodTests)


    # (i/m)Q, where i is the rank, m is the total number 
    # of tests, and Q is the false discovery rate you choose.
    sigLevel <-  iRank * FDR /mTotalTests

    sigTests <- which(dataVec < sigLevel)

    # The largest P value that has P<(i/m)Q is significant, and all of the P values smaller than
    # it are also significant, even the ones that aren't less than their
    # Benjamini-Hochberg critical value.

    maxTest <- max(dataVec[sigTests])

    finalSig <- which(dataVec < maxTest)

    sigInds <- goodTests[finalSig]

    finalVec <- rep(NA,length(tempField))
    finalVec[sigInds] <- 1
    sigMat <- matrix(finalVec, length(lon), length(lat))

    if(!is.null(plotName)){
        pdf(plotName,6,6)
        plot(dataVec, sigLevel,xlim=c(0,1), ylim=c(0,1), xlab= 'P-value', ylab= 'Significance Cutoff',
            main=paste('FDR =',FDR,'Cutoff =', round(maxTest,4)))
        points(dataVec[finalSig], sigLevel[finalSig], col = 'blue')
        points(dataVec[sigTests], sigLevel[sigTests], col='red')
        abline(0,1)
        dev.off()
    }
    return(sigMat)
}

# F function for comparing two forecasts in the GROC analysis
Ffunction<- function(forecastA,forecastB){
	(forecastA[2]*forecastB[1] + forecastA[3]*forecastB[1] + forecastA[3]*forecastB[2])/
		(1 - sum(forecastA*forecastB))
}

# helper function in reliability calculation to bin and properly area weight
fieldReliability <- function(fcast,obs,terc,lat,binSize=0.05){
	# make appropriate weight array
	cosMat <- cos(matrix(lat,nrow=nrow(fcast),ncol=ncol(fcast),byrow=T)*(pi/180))
	weights <- array(cosMat,dim=dim(fcast))

	# make bins
	bins <- seq(0-(binSize)/2,1+(binSize)/2,by=binSize)
	centers <- bins[-1]-(binSize)/2


	# figure out which obs were correct
	obsInds <- which(obs==(terc-2))

	# determine the bin counts
	nk <- rep(NA, length(centers))
	fcastBar <- rep(NA,length(centers))
	wt <- rep(NA, length(centers))
	obsBar <- rep(NA, length(centers))


	# Do first bin outside loop to include zero
	tempInds <- which(fcast>=bins[1] & fcast<=bins[2])
	nk[1] <- length(tempInds)
	fcastBar[1] <- sum(fcast[tempInds]*weights[tempInds])/sum(weights[tempInds])
	obsBar[1] <- sum(weights[intersect(tempInds,obsInds)])/sum(weights[tempInds])


	# loop through the rest of the bins
	for(i in 2:length(nk)){
		tempInds <- which(fcast>bins[i] & fcast<=bins[i+1])
		nk[i] <- length(tempInds)
		fcastBar[i] <- sum(fcast[tempInds]*weights[tempInds])/sum(weights[tempInds])
		obsBar[i] <- sum(weights[intersect(tempInds,obsInds)])/sum(weights[tempInds])
	}

	return(cbind(centers,fcastBar,obsBar,nk))
}

# make a quick field plot over a lon,lat range
fieldSkillPlot <- function(lon,lat,z,plotdir,plotName,...){

	zr <- quantile(z,probs=c(0.005,.995),na.rm=T)
	zm <- max(abs(zr))

	palFull <- redBlue(512)
	colMap <- seq(-zm,zm,length=length(palFull))

	minColInd <- which.min(abs(colMap-zr[1]))
	maxColInd <- which.min(abs(colMap-zr[2]))

	pal <- palFull[minColInd:maxColInd]

	pdf(sprintf('%s/%s.pdf',plotdir,plotName),10,6)
	image.plot(lon,lat,z,zlim=zr,col=pal,...)
	world(add=T)
	dev.off()
}

# Extrat metadata from an IRI seasonal forecast
getMeta <- function(fileInd,fileList){
	year  <- as.numeric(substr(fileList[fileInd],4,7))
	season <- as.numeric(substr(fileList[fileInd],9,10))
	lead  <- as.numeric(substr(fileList[fileInd],2,2))
	
	out <- c(year,season,lead)	
	names(out) <- c('year','season','lead')

	return(out)
}

# A function that calculates the GROC at each location point over the entire time series, returning
# a field series of GROC over the forecast spatial domain
grocField <- function(fcast,obs,tropics=TRUE){
	require(fields)

	# extract objects from the forecast list
	lon  <- fcast$lon
	nlon <- length(lon)
	lat  <- fcast$lat
	nlat <- length(lat)

	tropicsLats <- which(lat < 30 & lat > -30)


	# the resulting series that is constructed in the loop over time
	grocField <- matrix(NA, nlon,nlat)

	gridList <- expand.grid(list(lon=1:nlon,lat=1:nlat))

	if(tropics){
		tropicsRows <- which(gridList[,2] %in% tropicsLats)
		gridList <- gridList[tropicsRows,]
	}

	# match the time series
	forecastTime <- fcast$timeMap[,3]
	obsTime      <- obs$timeMap[,3]
	forecastTimeInds <- which(forecastTime %in% obsTime)
	obsTimeInds  <- which(obsTime %in% forecastTime)


	pb   <- txtProgressBar(1, nrow(gridList), style=3)

	# locInd <- 5086 # with tropics FALSE!
	for(locInd in 1:nrow(gridList)){
		setTxtProgressBar(pb,locInd)

		lonInd <- gridList[locInd,1]
		latInd <- gridList[locInd,2]

		# pull the correct objects for the given location
		tempForecast <- fcast$forecast[lonInd,latInd,forecastTimeInds,1:3]
		tempObs  <- obs$tercile[lonInd,latInd,obsTimeInds]

		if(all(is.na(tempForecast)) | all(is.na(tempObs))) next

		denominatorD <- 0
		score        <- 0

		for(outcome in 1:0){
			outcomeInds <- which(tempObs == outcome)
			compInds    <- which(tempObs < outcome)

			outcomeForecastsFull <- tempForecast[outcomeInds,]
			compForecastsFull    <- tempForecast[compInds,]
			
			# clean the matricies to allow for not having a forecast every year
			if(is.matrix(outcomeForecastsFull)){
				if(all(is.na(outcomeForecastsFull))){
					next
				} else{
					goodRows <- which(!is.na(outcomeForecastsFull[,1]))
				}
			} else{
				if(all(is.na(outcomeForecastsFull))){
					next
				} else{
					goodRows <- 1
					outcomeForecastsFull <- matrix(outcomeForecastsFull,nrow=1,ncol=3)
				}
			}

			if(is.matrix(compForecastsFull)){
				if(all(is.na(compForecastsFull))){
					next
				} else{
					goodRowsComp <- which(!is.na(compForecastsFull[,1]))
				}
			} else{
				if(all(is.na(compForecastsFull))){
					next
				} else{
					goodRowsComp <- 1
					compForecastsFull <- matrix(compForecastsFull,nrow=1,ncol=3)
				}
			}
			

			outcomeForecasts       <- outcomeForecastsFull[goodRows,]
			compForecasts          <- compForecastsFull[goodRowsComp,]


			if(!is.matrix(outcomeForecasts)){
				outcomeForecasts <- matrix(outcomeForecasts,nrow=1,ncol=3)
			}
			if(!is.matrix(compForecasts)){
				compForecasts <- matrix(compForecasts,nrow=1,ncol=3)
			}

			denominatorD <- denominatorD + nrow(outcomeForecasts) * nrow(compForecasts)


			if(nrow(outcomeForecasts) > 0 & nrow(compForecasts) > 0){
				for(j in 1:nrow(outcomeForecasts)){
					forecastA <- outcomeForecasts[j,]
					for(k in 1:nrow(compForecasts)){
						forecastB <- compForecasts[k,]
						score <- score + Ifunction(forecastA,forecastB)
					}
				}
			}

		}

		grocField[lonInd,latInd] <- score/denominatorD
	}

	return(grocField)
}

grocFieldPlot <- function(lon,lat,z,plotdir,plotName,...){

	bins <- c(0,0.5,0.51,0.55,0.6,0.65,0.7,1)
	pal  <- c('gray40', tim.colors(length(bins)-2))

	pdf(sprintf('%s/%s.pdf',plotdir,plotName),10,6)
	image(lon,lat,z,breaks=bins,col=pal,ylim=c(-60,90),...)
	world(add=T)
	dev.off()
}

# A function that calculates the GROC at each time point over the field, returning
# a time series of GROC over the forecast time domain
grocSeries <- function(fcast,obs,tropics=TRUE){
	# extract objects from the forecast list
	lon  <- fcast$lon
	nlon <- length(lon)
	lat  <- fcast$lat
	nlat <- length(lat)

	nForecasts <- nrow(fcast$timeMap)

	tropicsLats <- which(lat < 30 & lat > -30)

	# the resulting series that are constructed in the loop over time
	scoreSeries <- rep(NA, nForecasts)
	denomSeries <- rep(NA, nForecasts)
	grocSeries  <- rep(NA, nForecasts)


	pb   <- txtProgressBar(1, nForecasts, style=3)

	for(tInd in 1:nForecasts){
		setTxtProgressBar(pb,tInd)

		# get the correct time series index in the observation
		outcomeInd <- which(obs$timeMap[,3] == fcast$timeMap[tInd,3])
		if(length(outcomeInd) <1) next

		# pull the correct objects for the given time step
		tempForecast <- fcast$forecast[,,tInd,]
		tempObs  <- obs$tercile[,,outcomeInd]

		if(tropics){
			tempForecast <- fcast$forecast[,tropicsLats,tInd,]
			tempObs  <- obs$tercile[,tropicsLats,outcomeInd]
		}

		denominatorD <- 0
		score        <- 0
		# loop over the three possible outcomes
		for(outcome in 1:0){
			outcomeMask <- which(tempObs == outcome)
			compMask    <- which(tempObs < outcome)

			outcomeForecastsFull   <- matrix(c(tempForecast),nrow=nrow(tempObs)*ncol(tempObs),ncol=3)[outcomeMask,]
			compForecastsFull      <- matrix(c(tempForecast),nrow=nrow(tempObs)*ncol(tempObs),ncol=3)[compMask,]
			goodRows 		       <- which(!is.na(outcomeForecastsFull[,1]))
			goodRowsComp   		   <- which(!is.na(compForecastsFull[,1]))

			outcomeForecasts       <- outcomeForecastsFull[goodRows,]
			compForecasts          <- compForecastsFull[goodRowsComp,]

			# increment the denomiantor by the number of pairs in this group
			denominatorD <- denominatorD + nrow(outcomeForecasts) * nrow(compForecasts)
			
			for(j in 1:nrow(outcomeForecasts)){
				forecastA <- outcomeForecasts[j,]
				for(k in 1:nrow(compForecasts)){
					forecastB <- compForecasts[k,]
					score <- score + Ifunction(forecastA,forecastB)
				}
			}
		}

		scoreSeries[tInd] <- score
		denomSeries[tInd] <- denominatorD

		grocSeries[tInd] <- score/denominatorD
	}

	totalGroc <- sum(scoreSeries)/sum(denomSeries)

	return(list(series=grocSeries,total=totalGroc))
}

# I function for the GROC analysis
Ifunction <- function(forecastA,forecastB){
	if(all(forecastA == .33)) forecastA <- rep(1/3,3)
	if(all(forecastB == .33)) forecastB <- rep(1/3,3)

	if(all(forecastA == forecastB)) return(0.5)

	Fval <- Ffunction(forecastA,forecastB)

	if(Fval == 0.5) return(0.5)
	else if(Fval < 0.5) return(0)
	else return(1)
}

# plot a seasonal Forecast that looks like the IRI forecast
iriPlot <- function(lon,lat,fcastArr,oceanMask,...){
	sigLow  <- fcastArr[,,1]
	sigLow[sigLow < 0.4]  <-  NA

	sigHigh  <- fcastArr[,,3]
	sigHigh[sigHigh < 0.4]  <-  NA

	probBreaks <- c(40,45,50,60,70,100)/100

	# NA all non ocean vals
	ocean <- oceanMask
	ocean[ocean==0] <- NA

	dry <- NA
	
	# NA all non -1 vals
	if(dim(fcastArr)[3]>3){
		dry <- fcastArr[,,4]
		dry[dry==0] <- NA
	}


	#grab all of the colors
	aboveColor   <- c(brewer.pal(7,'GnBu'))[3:7]
	belowColor   <- c(brewer.pal(7,'YlOrBr'))[3:7]
	dryColor     <- adjustcolor('darkred',alpha=0.15)
	oceanColor   <- adjustcolor('paleturquoise1',alpha=0.2)
	neutralColor <- 'white'

	# build the plot through adding layers
	image(lon,lat,sigLow, col=belowColor, breaks=probBreaks, xlab='',ylab='',...)
	image(lon,lat,sigHigh, col=aboveColor, breaks=probBreaks, add=T)
	image(lon,lat,ocean, col=oceanColor,add=T)
	if(!all(is.na(dry)) & dim(fcastArr)[3]>3){
		image(lon,lat,dry,col=dryColor,add=T)
	}

	world(add=T)
}

log2Mean <- function(x,weights=NULL){
	if(is.null(weights)) weights <- rep(1,length(x))
	goodInds <- !is.na(x)
	-weighted.mean(log(x[goodInds],base=2),weights[goodInds])
}

# Convert IRI raw ASCII data into an R matrix
makeTable <- function(fileInd, fileList, ddir){
	rawData <- readLines(sprintf('%s/%s',ddir,fileList[fileInd]))
	nxy <- length(rawData)

	datMat <- matrix(NA, nrow=nxy,ncol=5)
	colnames(datMat) <- c('lon','lat','low','mid','high')
	for(i in 1:nxy){
		lat    <- as.numeric(substr(rawData[i],1,6))
		lon    <- as.numeric(substr(rawData[i],7,13))
		low    <- as.numeric(substr(rawData[i],15,16))
		mid	   <- as.numeric(substr(rawData[i],19,20))
		high   <- as.numeric(substr(rawData[i],23,24))

		tempRow <- c(lon,lat,low,mid,high)

		datMat[i,] <- tempRow
	}

	return(datMat)
}

# Function to generate a nice diverging colorbar
redBlue <- function(n=256){
	require(fields)
	require(RColorBrewer)
	rev(designer.colors(n, brewer.pal(11,'RdBu')))
}

# Function that calculates the required output for a reliability diagram
reliability <- function(fcast,obs,binSize=0.05){
	# extract objects from the forecast list
	lat  <- fcast$lat
	tropicsLats <- which(lat < 30 & lat > -30)

	# match the two time series
	tIndObs <- which(obs$timeMap[,3] %in% fcast$timeMap[,3])
	subObs <- observedTercile$tercile[,,tIndObs]

	resultList <- list()
	resultListTropics <- list()


	for(terc in 1:3){
		tempForecast <- fcast$forecast[,,,terc]
		tempForecastTropics <- tempForecast[,tropicsLats,]

		resultList[[terc]] <- fieldReliability(tempForecast,subObs,terc,lat,binSize)
		resultListTropics[[terc]] <- fieldReliability(tempForecastTropics,subObs[,tropicsLats,],terc,lat[tropicsLats],binSize)

	}

	outList <- list(global=resultList,tropics=resultListTropics)
}

# function to plot a reliability diagram
reliabilityDiagram <- function(dat,yr=NULL,textSize=1.5,...){

	pal <- c('brown','darkgreen','blue')
	transPal <- adjustcolor(pal,alpha=0.3)
	transPal2 <- adjustcolor(pal,alpha=0.5)

	barArrange <- matrix(NA,nrow=nrow(dat[[1]]),ncol=3)
	for(i in 1:3) barArrange[,i] <- dat[[i]][,4]


	if(is.null(yr)){
		yr <- c(0,max(barArrange))*1.25
	}

	# plot the barplot of counts
	par(mar=c(5,5,4,5) + 0.1)
	barplot(t(barArrange),beside=T,axes=FALSE,ylim=yr,col=transPal)
	axis(4,cex.axis=textSize,cex.lab=textSize)
	
	par(new=TRUE) 

	# plot the reliability curves
	plot(x=NULL,xlim=c(0,1),ylim=c(0,1),
		xlab='Forecast Probability',ylab='Empirical Probability',
		cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize,...)
	

	mtext('# of Forecasts Issued', side=4, line=3,cex=1)

	grid()
	for(i in 1:3){
		rawDat <- dat[[i]]
		goodInds <- which(!is.na(rawDat[,2]))
		tempDat <- rawDat[goodInds,]

		points(tempDat[,2:3],type='b',col=pal[i],lwd=1.5)

		abline(lm(tempDat[,3] ~ tempDat[,2],weights = tempDat[,4]),lty=2,lwd=2,col=transPal2[i])

	}
	abline(0,1)


	legend('top',
		c('Dry Tercile', 'Normal Tercile', 'Wet Tercile'),
		lwd=2, col=pal,bg='white', cex = textSize, bty='n')

}

# Function that returns three objects where skill in measured by RPSS:
#	(1) A time series of global skill
#   (2) A time series of tropical skill
#	(3) A field of skill by location
rpss <- function(fcast1, obs, fcastBaseline=NULL,landPropMat=NULL){
	lon  <- fcast1$lon
	nlon <- length(lon)
	lat  <- fcast1$lat
	nlat <- length(lat)

	if(is.null(landPropMat)) landPropMat <- landPropList$prop

	nForecasts <- nrow(fcast1$timeMap)

	# store the RPS values at each point
	rpsForecast <- array(NA, dim=c(nlon,nlat,nForecasts))
	rpsBaseline <- array(NA, dim=c(nlon,nlat,nForecasts))

	#  lat-cosine scaling matrix to multiply before summing over the globe
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))

	#  lat-cosine scaling matrix to multiply before summing over the globe
	cosMat <- cos(matrix(lat,nrow=nlon,ncol=nlat,byrow=T)*(pi/180))

	# climatology forecast
	climatology <- rep(1/3,3)
	cumClimatology <- cumsum(climatology)

	# keep track of which alt forecasts used (use to debug later!)
	altIndVec <- c()

	# calculate how many forecasts happen at each grid cell
	nForecastMat <- apply(fcast1$forecast[,,,1],c(1,2), function(x) sum(!is.na(x)))

	# loop through  everything
	for(tInd in 1:nForecasts){
		outcomeInd <- which(obs$timeMap[,3] == fcast1$timeMap[tInd,3])
		tempForecast <- fcast1$forecast[,,tInd,]
		tempObs  <- obs$tercile[,,outcomeInd]

		if(!is.null(fcastBaseline)){
			altInd  <- which(fcastBaseline$timeMap[,3] == fcastBaseline$timeMap[tInd,3])
			altIndVec <- c(altIndVec,altInd)
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


	# calculate the total forecast and baseline values
	nForecastMat[nForecastMat ==0] <- NA
	rpsForecastTotal <- weighted.mean(rpsForecastField,w=nForecastMat * cosMat, na.rm=T)
	rpsBaselineTotal <- weighted.mean(rpsBaselineField,w=nForecastMat * cosMat, na.rm=T)

	rpsForecastTotalTropics <- weighted.mean(rpsForecastField[tropicsLats,],
									w=nForecastMat[tropicsLats,] * cosMat[tropicsLats,], na.rm=T)
	rpsBaselineTotalTropics <- weighted.mean(rpsBaselineField[tropicsLats,],
									w=nForecastMat[tropicsLats,] * cosMat[tropicsLats,], na.rm=T)


	rpssTotalGlobal <- 1 - rpsForecastTotal/rpsBaselineTotal
	rpssTotalTropics <- 1 - rpsForecastTotalTropics/rpsBaselineTotalTropics


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
	global  		<- list(field=rpssField,series=rpssGlobalSeries, total = rpssTotalGlobal)
	tropics 		<- list(series=rpssTropicsSeries, total = rpssTotalTropics)

	outList <- list(lon=lon,lat=lat,timeMap=fcast1$timeMap,
					global=global,tropics=tropics,
					bestMat = bestMat, 
					bestSeason = bestSeason,
					seasonalArr = seasonArr)

	return(outList)
}

# Function that decomposes the RPS in two different ways and calculates rpss
# and a reliability skill score
rpssDecompField <- function(fcast1,obs,tol=0.001,minForecasts=10){

	lon  <- fcast1$lon
	nlon <- length(lon)
	lat  <- fcast1$lat
	nlat <- length(lat)

	# climatology forecast
	climatology <- rep(1/3,3)
	cumClimatology <- cumsum(climatology)


	outcomeInds <- which(obs$timeMap[,3] %in% fcast1$timeMap[,3])

	rpsMat73  <- matrix(NA, nrow=nlon, ncol=nlat)
	rpsMat86  <- matrix(NA, nrow=nlon, ncol=nlat)

	relMat  <- matrix(NA, nrow=nlon, ncol=nlat)
	resMat  <- matrix(NA, nrow=nlon, ncol=nlat)
	uncMat  <- matrix(NA, nrow=nlon, ncol=nlat)

	varMat  <- matrix(NA, nrow=nlon, ncol=nlat)
	biasMat <- matrix(NA, nrow=nlon, ncol=nlat)

	rpssMat  <- matrix(NA, nrow=nlon, ncol=nlat)
	relssMat <- matrix(NA, nrow=nlon, ncol=nlat)

	badCells <- matrix(nrow=0,ncol=5)

	for(i in 1:nlon){
		for(j in 1:nlat){
			tempForecastFull <- fcast1$forecast[i,j,,1:3]
			tempOutcomeFull  <- obs$tercile[i,j,outcomeInds] +2

			badRows <- which(is.na(tempForecastFull[,1]) | is.na(tempOutcomeFull))

			if(length(badRows)>0){
				tempForecast <- tempForecastFull[-badRows,]
				tempOutcome  <- tempOutcomeFull[-badRows]			
			} else{
				tempForecast <- tempForecastFull
				tempOutcome  <- tempOutcomeFull	
			}

			if(length(tempOutcome) < minForecasts) next

			climForecast <- matrix(1/3,nrow=nrow(tempForecast),ncol=3)

			####################################################################
			# Decompose the Brier score into three parst following [Murphy 1973]
			# using the notation of [Stephenson et al. 2008]
			#
			# Note: this agrees with the RPS calculation in the RPSS function
			####################################################################
			# Break forecast into two brier calculations [Make outcome 0/1!!!!!]
			lowForecast <- cbind(tempForecast[,1],tempForecast[,2] + tempForecast[,3])
			lowClim     <- cbind(climForecast[,1],climForecast[,2] + climForecast[,3])
			lowOutcome  <- ifelse(tempOutcome==1,1,2)-1

			highForecast <- cbind(tempForecast[,1] + tempForecast[,2], tempForecast[,3])
			highClim     <- cbind(climForecast[,1] + climForecast[,2], climForecast[,3])
			highOutcome  <- ifelse(tempOutcome==3,2,1) -1

			lowDecomp     <- brierDecomp10WeightedBinned(lowForecast,lowOutcome,weights=rep(1,length(lowOutcome)))
			highDecomp    <- brierDecomp10WeightedBinned(highForecast,highOutcome,weights=rep(1,length(highOutcome)))

			lowClimDecomp <- brierDecomp10WeightedBinned(lowClim,lowOutcome,weights=rep(1,length(lowOutcome)))
			highClimDecomp <- brierDecomp10WeightedBinned(highClim,highOutcome,weights=rep(1,length(highOutcome)))

			testBrierL <- mean((lowForecast[,2]-lowOutcome)^2)
			testBrierH <- mean((highForecast[,2]-highOutcome)^2)

			rps0 <- lowDecomp + highDecomp
			rpsClim <- lowClimDecomp + highClimDecomp

			rpsMat73[i,j] <- rps0[1]
			relMat[i,j]   <- rps0[2]
			resMat[i,j]   <- rps0[3]
			uncMat[i,j]   <- rps0[4]

			rpssMat[i,j]  <- 1 - rps0[1]/rpsClim[1]
			relssMat[i,j] <- 1 - rps0[2]/rpsClim[2]
			####################################################################
			# Decompose the Brier score into two parts following [Murphy 1986]
			#
			# Note: this does not exactly agree with the RPS calculation 
			# avove or in the RPSS function. It appears to agree less when
			# sample size is small.

			# all of the variance terms are very small, need to think about if
			# is an accurate calculation in this case. My guess is yes because
			# we are issusing so many climatology forecasts.
			####################################################################
			lowBrier86 <- brierDecomp86(lowForecast,lowOutcome)
			highBrier86 <- brierDecomp86(highForecast,highOutcome)

			# total RPS
			rps <- lowBrier86 + highBrier86
			rpsClim <- lowClimDecomp + highClimDecomp

			if((rps[1] - rps0[1])>tol){
				badCells <- rbind(badCells,c(i,j,rps0[1],rps[1],length(tempOutcome)))
			}

			rpsMat86[i,j] <- rps[1]
			varMat[i,j]   <- rps[2]
			biasMat[i,j]  <- rps[3]


		}
	}

	decomp73 <- list(rps = rpsMat73, rel = relMat, res = resMat, unc = uncMat)
	decomp86 <- list(rps = rpsMat86, var = varMat, bias = biasMat)
	skill    <- list(rpss=rpssMat, relss=relssMat)
	info     <- list(badCells=badCells)

	outList <- list(lon=lon,lat=lat,decomp73=decomp73,decomp86=decomp86,
					skill=skill,info=info)

	return(outList)
}

# Function that decomposes the RPS in according to Murphy 1973 following
# the notation of Young 2010 and for a global or tropics series
rpssDecompSeries <- function(fcast1,obs,dryMask,tropics=FALSE){
	lon  <- fcast1$lon
	nlon <- length(lon)
	lat  <- fcast1$lat
	nlat <- length(lat)

	nForecasts <- nrow(fcast1$timeMap)


	#  lat-cosine scaling matrix to multiply before summing over the globe
	latMat <- matrix(lat,nlon,nlat,byrow=T)
	cosMat <- cos(latMat*(pi/180))

	# build the proper scaling vectors with the dry mask and cosMat
	cosVec <- c(cosMat)

	scalingMat <- matrix(nrow=length(cosVec),ncol=12)

	for(i in 1:12){
		scalingMat[,i] <- cosVec*c(dryMask[,,i])
	}
	# climatology forecast
	climatology <- rep(1/3,3)
	cumClimatology <- cumsum(climatology)


	outcomeInds <- which(obs$timeMap[,3] %in% fcast1$timeMap[,3])

	rpsSeries73  <- rep(NA, nForecasts)
	rpsSeries86  <- rep(NA, nForecasts)

	relSeries  <- rep(NA, nForecasts)
	resSeries  <- rep(NA, nForecasts)
	uncSeries  <- rep(NA, nForecasts)

	varSeries  <- rep(NA, nForecasts)
	biasSeries <- rep(NA, nForecasts)

	rpssSeries  <- rep(NA, nForecasts)
	relssSeries <- rep(NA, nForecasts)

	badCells <- matrix(nrow=0,ncol=5)

	# BEGIN LOOP IN TIME
	for(tInd in 1:nForecasts){

		tempSeason <- fcast1$timeMap[tInd,2]

		tempForecastFull <- fcast1$forecast[,,tInd,1:3]
		tempForecastMat <- cbind(c(tempForecastFull[,,1]),
								 c(tempForecastFull[,,2]),
								 c(tempForecastFull[,,3]))

		tempOutcomeFull  <- c(obs$tercile[,,outcomeInds[tInd]] +2)

		# turn into vector format for analysis using old code
		if(tropics){
			badRows <- which(is.na(tempForecastMat[,1]) | 
							is.na(tempOutcomeFull) |
							c(abs(latMat))>30 | 
							scalingMat[,tempSeason]==0)
		} else{
			badRows <- which(is.na(tempForecastMat[,1]) | 
							is.na(tempOutcomeFull) | 
							scalingMat[,tempSeason]==0)
		}
		
		if(length(badRows)>0){
			tempForecast <- tempForecastMat[-badRows,]
			tempOutcome  <- tempOutcomeFull[-badRows]		
			tempWeights  <- scalingMat[-badRows,tempSeason]
		} else{
			tempForecast <- tempForecastMat
			tempOutcome  <- tempOutcomeFull
			tempWeights  <- scalingMat[-badRows,tempSeason]	
		}

		climForecast <- matrix(1/3,nrow=nrow(tempForecast),ncol=3)

		####################################################################
		# Decompose the Brier score into three parst following [Murphy 1973]
		# using the notation of [Stephenson et al. 2008]
		#
		# Note: this agrees with the RPS calculation in the RPSS function
		####################################################################
		# Break forecast into two brier calculations [Make outcome 0/1!!!!!]
		lowForecast <- cbind(tempForecast[,1],tempForecast[,2] + tempForecast[,3])
		lowClim     <- cbind(climForecast[,1],climForecast[,2] + climForecast[,3])
		lowOutcome  <- ifelse(tempOutcome==1,1,2)-1

		highForecast <- cbind(tempForecast[,1] + tempForecast[,2], tempForecast[,3])
		highClim     <- cbind(climForecast[,1] + climForecast[,2], climForecast[,3])
		highOutcome  <- ifelse(tempOutcome==3,2,1) -1

		lowDecomp     <- brierDecomp10WeightedBinned(lowForecast,lowOutcome,tempWeights)
		highDecomp    <- brierDecomp10WeightedBinned(highForecast,highOutcome,tempWeights)	

		lowClimDecomp  <- brierDecomp10WeightedBinned(lowClim,lowOutcome,tempWeights)
		highClimDecomp <- brierDecomp10WeightedBinned(highClim,highOutcome,tempWeights)

		rps0 <- lowDecomp + highDecomp
		rpsClim <- lowClimDecomp + highClimDecomp

		rpsSeries73[tInd] <- rps0[1]
		relSeries[tInd]   <- rps0[2]
		resSeries[tInd]   <- rps0[3]
		uncSeries[tInd]  <- rps0[4]

		rpsSeries73[tInd]  <- 1 - rps0[1]/rpsClim[1]
		relssSeries[tInd] <- 1 - rps0[2]/rpsClim[2]
	}
	# END LOOP

	decomp73 <- list(rps = rpsSeries73, rel = relSeries, res = resSeries, unc = uncSeries)
	skill    <- list(rpss=rpsSeries73, relss=relssSeries)

	outList <- list(timeMap=fcast1$timeMap,decomp73=decomp73,skill=skill)

	return(outList)
}

rpssFieldPlot <- function(lon,lat,z,plotdir,plotName,...){

	bins <- c(-1,0,0.01,0.05,0.1,0.2,0.3,1)
	pal  <- c('gray40', tim.colors(length(bins)-2))

	pdf(sprintf('%s/%s.pdf',plotdir,plotName),10,6)
	image(lon,lat,z,breaks=bins,col=pal,ylim=c(-60,90),...)
	world(add=T)
	dev.off()
}

# calculate the total rpss for all forecasts
rpssTotalDecomp <- function(fcast1,obs,dryMask,tropics=FALSE){
	lon  <- fcast1$lon
	nlon <- length(lon)
	lat  <- fcast1$lat
	nlat <- length(lat)

	nForecasts <- nrow(fcast1$timeMap)


	#  lat-cosine scaling matrix to multiply before summing over the globe
	latMat <- matrix(lat,nlon,nlat,byrow=T)
	cosMat <- cos(latMat*(pi/180))

	# build the proper scaling vectors with the dry mask and cosMat
	cosVec <- c(cosMat)

	scalingMat <- matrix(nrow=length(cosVec),ncol=12)

	for(i in 1:12){
		scalingMat[,i] <- cosVec*c(dryMask[,,i])
	}
	# climatology forecast
	climatology <- rep(1/3,3)
	cumClimatology <- cumsum(climatology)


	outcomeInds <- which(obs$timeMap[,3] %in% fcast1$timeMap[,3])

	lowForecast <- matrix(nrow=0,ncol=2)
	lowClim     <- matrix(nrow=0,ncol=2)
	lowOutcome  <- numeric()

	highForecast <- matrix(nrow=0,ncol=2)
	highClim     <- matrix(nrow=0,ncol=2)
	highOutcome  <- numeric()

	weights <- numeric()

	# BEGIN LOOP IN TIME
	for(tInd in 1:nForecasts){

		tempSeason <- fcast1$timeMap[tInd,2]

		tempForecastFull <- fcast1$forecast[,,tInd,1:3]
		tempForecastMat <- cbind(c(tempForecastFull[,,1]),
								 c(tempForecastFull[,,2]),
								 c(tempForecastFull[,,3]))

		tempOutcomeFull  <- c(obs$tercile[,,outcomeInds[tInd]] +2)

		if(tropics){
			badRows <- which(is.na(tempForecastMat[,1]) | 
							is.na(tempOutcomeFull) |
							c(abs(latMat))>30 | 
							scalingMat[,tempSeason]==0)
		} else{
			badRows <- which(is.na(tempForecastMat[,1]) | 
							is.na(tempOutcomeFull) | 
							scalingMat[,tempSeason]==0)
		}
		
		if(length(badRows)>0){
			tempForecast <- tempForecastMat[-badRows,]
			tempOutcome  <- tempOutcomeFull[-badRows]		
			tempWeights  <- scalingMat[-badRows,tempSeason]
		} else{
			tempForecast <- tempForecastMat
			tempOutcome  <- tempOutcomeFull
			tempWeights  <- scalingMat[-badRows,tempSeason]	
		}



		climForecast <- matrix(1/3,nrow=nrow(tempForecast),ncol=3)

		tempLowFcast    <- cbind(tempForecast[,1],tempForecast[,2] + tempForecast[,3])
		tempLowClim     <- cbind(climForecast[,1],climForecast[,2] + climForecast[,3])
		tempLowOutcome  <- ifelse(tempOutcome==1,1,2)-1

		tempHighFcast    <- cbind(tempForecast[,1] + tempForecast[,2], tempForecast[,3])
		tempHighClim     <- cbind(climForecast[,1] + climForecast[,2], climForecast[,3])
		tempHighOutcome  <- ifelse(tempOutcome==3,2,1) -1


		# append the fcast/outcomes to the list
		lowForecast <- rbind(lowForecast,tempLowFcast)
		lowClim     <- rbind(lowClim,tempLowClim)
		lowOutcome  <- c(lowOutcome,tempLowOutcome)

		highForecast <- rbind(highForecast,tempHighFcast)
		highClim     <- rbind(highClim,tempHighClim)
		highOutcome  <- c(highOutcome,tempHighOutcome)

		weights <- c(weights,tempWeights)
	}

	# calculate the relevant decomposition and forecast statistics

	lowDecomp     <- brierDecomp10WeightedBinned(lowForecast,lowOutcome,weights)
	highDecomp    <- brierDecomp10WeightedBinned(highForecast,highOutcome,weights)	

	lowClimDecomp  <- brierDecomp10WeightedBinned(lowClim,lowOutcome,weights)
	highClimDecomp <- brierDecomp10WeightedBinned(highClim,highOutcome,weights)

	rps <- lowDecomp + highDecomp
	rpsClim <- lowClimDecomp + highClimDecomp

	rpss  <- 1 - rps[1]/rpsClim[1]

	return(list(rpss=rpss,rps=rps[1],rel=rps[2],res=rps[3],unc=rps[4]))
}

# averages a monthly time series into a seasonal time series
seasonalAverage12 <- function(ts,tYear){

	nMonth <- length(ts)
	nYear <- length(tYear)

	outMat <- matrix(NA, nrow=length(tYear)-1, ncol=12)
	
	rownames(outMat) <- tYear[2:nYear]
	colnames(outMat) <- c('DJF', 'JFM', 'FMA', 
						  'MAM', 'AMJ', 'MJJ',
						  'JJA', 'JAS', 'ASO',
						  'SON', 'OND', 'NDJ')

	for(i in 2:nYear){
		tsInd <- 1 + 12*(i-1)

		for(j in 1:12){
			seasonInds <- ( (j-2) : j)
			outMat[i-1,j] <- mean(ts[tsInd+seasonInds])
		}
	}



	return(outMat)
}

# applies a function to convert a monthly field into a seasonal field
seasonalField12 <- function(arr,tYear,FUN){
	nlon <- dim(arr)[1]
	nlat <- dim(arr)[2]

	outTime <- (length(tYear)-1) * 12 # not sure about this dimensionality rn

	outArray <- array(NA, dim = c(nlon,nlat,outTime))

	nYear <- length(tYear) - 1

	pb   <- txtProgressBar(2, nYear, style=3)

	for(i in 2:nYear){
		setTxtProgressBar(pb, i)

		tsInd <- 1 + 12*(i-1)

		for(j in 1:12){
			seasonInds <- ( (j-2) : j)
			if(max(tsInd+seasonInds) <= dim(arr)[3]){
				outArray[,,((i-2)*12 + j)] <- apply(arr[,,tsInd+seasonInds],c(1,2),FUN)
			}
		}
	}

	return(outArray)
}

# plot the skill of multiple
seriesSkillPlot <- function(x,s1,s2,s3,xTrim,s4,plotdir,plotName,ylimInc=NA,...){
	yr <- range(s1,s3,ylimInc,na.rm=T)
	pdf(sprintf('%s/%s.pdf',plotdir,plotName),10,6)
	plot(x,s1,type='l',col='red',ylim=yr,xlab='Year',...)
	points(x,s2,type='l',col='blue')
	points(x,s3,type='l',col='black',lwd=1.5)
	points(xTrim,s4,type='l',col='darkgreen',lwd=1)

	abline(h=c(0,0.5))
	legend('bottomleft',
		c('IRI Forecast','Probabilistic ENSO','Deterministic ENSO','Realistic ENSO'),
		col=c('black','blue','red','darkgreen'),lwd=1.5)
	dev.off()
}


# takes a forecast list object and pulls out specific months
subsetForecastList <- function(forecastList, inds){
	subInds <- which(forecastList$timeMap[,2] %in% inds)

	timeMap <- forecastList$timeMap[subInds,]
	fcast   <- forecastList$forecast[,,subInds,]

	listOut <- list(lon=forecastList$lon, lat=forecastList$lat,
					forecast=fcast, timeMap=timeMap)

	return(listOut)
}