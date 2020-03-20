# load in the processed enso series
load('Data/RawProcessed/ninaSeasonalOni.Rda')

# arrange in a matrix for more intuitive calculations
ninaMat <- matrix(c(ninaSeasonal$nina34), ncol=12,byrow=TRUE)

ninaIndicatorTemp <- matrix(0,nrow=nrow(ninaMat),ncol=ncol(ninaMat))

ninaIndicatorTemp[ninaMat >= cutoff]  <- 1
ninaIndicatorTemp[ninaMat <= -cutoff] <- -1


# Use the NOAA definition (need 5 consecutive months over cutoff)
# Use the JJA(7)-MJJ(6) as the year for determining max value
ninaIndVec  <- c(t(ninaIndicatorTemp))
oniIndVec   <- rep(0, length(ninaIndVec))
strengthVec <- rep(0, length(ninaIndVec))

monthVec <- ninaSeasonal$timeMap[,2]

count <- 1
state <- ninaIndVec[1]
maxVal <- abs(ninaSeasonal$nina34[1])

for(i in 2:length(ninaIndVec)){
	nextState  <- ninaIndVec[i]
	nextNino34 <- ninaSeasonal$nina34[i]

	if(nextState == state){
		count <- count + 1

		if(monthVec[i] > 6 && count > 9){
			maxVal <- 0
			count  <- 1
		}
		maxVal <- max(maxVal, abs(nextNino34))

	} else{
		count <- 1
		maxVal <- 0
	}

	if(count >= consecutiveMonths & state != 0){
		oniIndVec[(i-count+1):i] <- nextState
		strengthVec[(i-count+1):i] <- maxVal
	}

	state <- nextState
}


maxOniIndVec <- oniIndVec
maxOniIndVec[strengthVec<strengthCutoff] <- 0

# Test print to visually inspect the results to make sure it's all ok
# cbind(oniIndVec,maxOniIndVec,strengthVec,ninaSeasonal$nina34,
#		ninaSeasonal$timeMap[,2])

# package in the legacy matrix form
oniIndicator <- matrix(maxOniIndVec, ncol=12,byrow=TRUE)

# we do the calculation for each season (column)

ninoList <- list()
ninaList <- list()

# make sure we are working with the correct year range
tYearTemp <- min(ninaSeasonal$timeMap[,1]):max(ninaSeasonal$timeMap[,1])
timeInds  <- which(tYearTemp %in% startYear:endYear)

for(i in 1:ncol(oniIndicator)){

	ninoInds <- which(oniIndicator[timeInds,i]== 1)
	ninaInds <- which(oniIndicator[timeInds,i]== -1)

	# get the values of each
	ninoMatTemp <- cbind(tYearTemp[ninoInds],ninaMat[ninoInds,i])
	ninaMatTemp <- cbind(tYearTemp[ninaInds],ninaMat[ninaInds,i])

	ninoList[[i]] <- ninoMatTemp[order(ninoMatTemp[,2],decreasing=TRUE),]
	ninaList[[i]] <- ninaMatTemp[order(ninaMatTemp[,2]),]
}

# Needed for next steps of analysis:
# nina Indicator: a matrix with (-1,0,1) values for LN/N/EN [nYear,12]

# New year/strength tracking:
# nino/aList: a list with an element for each season with year and nino3.4
# values, sorted by strength [[12]]

# Old year/strength tracking:
# ensoArray: the ENSO statistic values of the top years [topYears,12,2]
# yearArray: the years of the top values for each month [topYears,12,2]

ninaIndicator <- oniIndicator[timeInds,]
rownames(ninaIndicator) <- startYear:endYear

save(ninaIndicator,ninoList,ninaList,file=sprintf('%s/ninaSeasonalInd.Rda',ddir))