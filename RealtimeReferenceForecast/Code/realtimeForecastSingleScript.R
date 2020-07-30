library(ncdf4)

###############################################################################
# Variables that may need to be changed
###############################################################################
url <- 'https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt'
ddir <- 'Data/Realtime'
dlMethod <- 'curl'

# ENSO event detection parameters
cutoff <- 0.5
strengthCutoff <- 1
consecutiveMonths <- 5

# forecast parameters
ddirImpacts <- 'Data/Impacts/Cru_2.5_1951_2016_Data'
forecastNinoCutoff <- 0.5
###############################################################################
# Pull and Clean Data from NOAA
###############################################################################

# download the latest file
download.file(url,sprintf('%s/realtimeONI.txt',ddir),method=dlMethod,quiet=TRUE)

# read the text file into R
oniTable  <- read.table(sprintf('%s/realtimeONI.txt', ddir),
					header=TRUE, stringsAsFactors=FALSE)

oniSeries <- oniTable[,4]

# make a timeMap object
timeMap <- cbind(oniTable[,2],NA,NA)
timeMap[,2] <- rep(1:12,length=nrow(timeMap))
timeMap[,3] <- timeMap[,1] + (timeMap[,2]-1)/12
colnames(timeMap) <- NULL



# figure out when events occur
ninaIndVec <- rep(0, length(oniSeries))
ninaIndVec[oniSeries >= cutoff]  <- 1
ninaIndVec[oniSeries <= -cutoff] <- -1

oniIndVec   <- rep(0, length(ninaIndVec))
strengthVec <- rep(0, length(ninaIndVec))

count <- 1
state <- ninaIndVec[1]
maxVal <- abs(oniSeries[1])

for(i in 2:length(ninaIndVec)){
	nextState  <- ninaIndVec[i]
	nextNino34 <- oniSeries[i]

	if(nextState == state){
		count <- count + 1

		if(timeMap[i,2] > 6 && count > 9){
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

# screen the event vector by strength of event
eventVec <- oniIndVec
eventVec[strengthVec<strengthCutoff] <- 0

###############################################################################
# Make the EBF
###############################################################################

# get the results list object from the impacts analysis
load(sprintf('%s/maskedAnalysisResults.Rda',ddirImpacts))

# unpack the results
lon 	    <- resultsList$lon
lat         <- resultsList$lat
sampleSize  <- resultsList$sampleSize
sigCounts   <- resultsList$sigCounts
ensoYears   <- resultsList$ensoYears
countsArray <- resultsList$countsArray


# array to fill
ensoForecastArray <- array(NA, dim=c(length(lon),
									 length(lat),
									 nrow(timeMap),3))

# need a different format of the eventVec
isNina <- eventVec
isNina[isNina== -1] <- 2

# optional cutoff for only strong el nino/la nina events
isStrongVec <- abs(oniSeries) >= forecastNinoCutoff
isNina[!isStrongVec] <- 0

# setup other things
climatology <- rep(1/3,3)
landMask <- ifelse(is.na(sampleSize[,,1]),NA,1)

# build a forecast for each time point
for(tInd in 1:length(isNina)){

	# determine the oracle ENSO state for the season of the forecast
	ninaState <- isNina[tInd]
	
	# build objects for the loop
	tempArray <- array(NA, dim=c(length(lon),length(lat),3))
	tempArrayProb <- tempArray

	if(ninaState==0){
		for(i in 1:3){
			tempArray[,,i] <- landMask  * climatology[i]
			tempArrayProb[,,i] <- landMask  * climatology[i]
		}
	} else{
		
		tempSeason <- timeMap[tInd,2]

		highLow <- sigCounts[,,tempSeason,ninaState,]
		highLow[highLow %in% c(-1,0)] <- NA

		countsTemp <- countsArray[,,tempSeason,ninaState,]
		numYearsTemp <- ensoYears[,,tempSeason,ninaState]

		for(i in 1:length(lon)){
			for(j in 1:length(lat)){
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
			tempArray[,,k] <- tempArray[,,k] * landMask 
			tempArrayProb[,,k] <- tempArrayProb[,,k] * landMask
		}
	}

	ensoForecastArray[,,tInd,] <- tempArrayProb
}


###############################################################################
# Save a netcdf file
###############################################################################
# need to build time dimension as months since 1960

position0 <- which(timeMap[,1]==1960 & timeMap[,2]==1)
timeVec <- -(position0-1):(nrow(timeMap)-position0)+0.5

# define dimensions
londim   <- ncdim_def("lon","degrees_east",as.double(lon)) 
latdim   <- ncdim_def("lat","degrees_north",as.double(lat)) 
timeDim   <- ncdim_def("season",'center month of 3-month season',as.double(timeVec))

terc    <- ncdim_def("tercile",'Above-Normal/Near-Normal/Below-Normal',as.double(1:3))

# define variables
fillvalue <- 1e32
dlname <- "Forecast Probability of Tercile"
fcast_def <- ncvar_def("fcast",'EBF Reference Tericle Precipitation Forecast',list(londim,latdim,timeDim,terc),
						fillvalue,dlname,prec="single")

# create netCDF file and put arrays
ncfname <- sprintf('%s/referenceEBF.nc',ddir)
ncout <- nc_create(ncfname,list(fcast_def),force_v4=T)

# put variables
ncvar_put(ncout,fcast_def,ensoForecastArray)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"season","axis","T")

# add global attributes
ncatt_put(ncout,0,"title",'Lenssen et al. 2020 Results')
ncatt_put(ncout,0,"institution",'IRI')
history <- paste("N. Lenssen", date(), sep=", ")
ncatt_put(ncout,0,"history",history)

# close the file, writing data to disk
nc_close(ncout)


# send the data to the IRI serve
system(sprintf('scp %s lenssen@orca01:/data/lenssen/realtimeAlternativeForecasts/',ncfname))








