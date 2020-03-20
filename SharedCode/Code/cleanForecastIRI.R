# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# pull all of the file names from the forecast dir. We are only working with the 2.5 deg
# forecasts for now
fullFileList <- system('ls Data/Raw/iriData',intern=TRUE)

highResFileInds <- grep('1deg',fullFileList,value=FALSE)
xfileInds        <- grep('*x', fullFileList, value= FALSE)
fileList <- fullFileList[-c(highResFileInds,xfileInds)]



# load in all of the forecasts

nFile <- length(fileList)

# build the array that stores the forecast data. We have an additional category to 
# track the location of the -1 flag
ForecastArray <- array(NA, dim=c(LON_SIZE,LAT_SIZE,nFile,FCAST_CATS+1))

metaMap <- matrix(NA, nrow=nFile, ncol=3)
colnames(metaMap) <- c('year', 'season', 'lead')


# The year and season is the target, the lead is how many months before
# the target date the forecast was issued. Therefore, we line up this year
# and season with the actual obs. We can also then plot it as a function of
# lead time.
badInds <- NULL


pb   <- txtProgressBar(1, nFile, style=3)

for(i in 1:nFile){
	setTxtProgressBar(pb, i)

	tempStr <- fileList[i]

	# put everything in a table (.25 sec)
	tab <- makeTable(i,fileList,'Data/Raw/iriData')

	# get the lon lat lists (2.5 x 2.5 degree)
	lon  <- sort(unique(tab[,1]))
	nlon <- length(lon)
	lat  <- sort(unique(tab[,2]))
	nlat <- length(lat)

	# if we don't have the correct number of lines, skip and save the ind
	if(nlon != LON_SIZE | nlat != LAT_SIZE){
		badInds <- c(badInds,i)
		next
	}
	
	if(i==1){
		# get the row inds to map the IRI tab to a fields tab (a few secs, only need once)
		targetGrid <- make.surface.grid(list(x=lon,y=lat))
		indMap <- apply(targetGrid,1,function(x) which(x[1] == tab[,1] & x[2] == tab[,2]))

		# get a ocean mask
		oceanMask <- matrix(tab[indMap,3],nrow=nlon,ncol=nlat) == -9
	}

	# make an matrix
	lowMat  <- matrix(tab[indMap,3],nrow=nlon,ncol=nlat)
	midMat  <- matrix(tab[indMap,4],nrow=nlon,ncol=nlat)
	highMat <- matrix(tab[indMap,5],nrow=nlon,ncol=nlat)

	# get the dry mat for the fcast
	dryMat  <- matrix(tab[indMap,3],nrow=nlon,ncol=nlat) == -1

	# combine into an array
	ForecastArray[,,i,] <- array(c(lowMat,midMat,highMat,dryMat),dim=c(nlon,nlat,4))

	# extract the metadata where the row of the matrix aligns with the
	# position on the list
	metaMap[i,] <- getMeta(i,fileList)

}

# Build the final time map matrix 
continuousTime <- metaMap[,1] + (metaMap[,2]-1)/12
timeMap <- cbind(metaMap[,1:2],continuousTime,metaMap[,3])

colnames(timeMap) <- c('year', 'season', 'continuousTime', 'lead')

# build the initial forecastList
forecastList <- list(lon=lon,lat=lat,forecast=ForecastArray,timeMap=timeMap)

###############################################################################
# Create a mask of forecasts in consistent locations to handle early hand drawn
# maps. We do this my grabbing locations that have greater than X forecasts
# issued over the forecast record
###############################################################################

lead <- 1

leadInds <- which(forecastList$timeMap[,4]==lead)


# make the correct regions in the forecast NA
forecastNA <- forecastList$forecast
forecastNA[forecastNA<0] <- NA

forecastLead1 <- list(lon=forecastList$lon,
					  lat=forecastList$lat,
					  forecast=forecastNA[,,leadInds,],
					  timeMap=forecastList$timeMap[leadInds,])

countMap <- matrix(NA, nrow=length(lon), ncol=length(lat))

for(i in 1:nrow(countMap)){
	for(j in 1:ncol(countMap)){
		series <- forecastLead1$forecast[i,j,,1]
		countMap[i,j] <- sum(series>0,na.rm=T)
	}
}

# build the mask according to the cutoff
cutoffMat <- countMap
cutoffMat[cutoffMat < forecastMaskCutoff] <- NA

# mask out antartica as well as it only recieves climatology forecasts
antarticaMask <- which(lat< -60)
cutoffMat[,antarticaMask] <- NA

# convert to a mask
masterForecastMask <- !is.na(cutoffMat)
masterForecastMask[masterForecastMask==0] <- NA


###############################################################################
# convert to probabilities and fix the climatology forecasts
###############################################################################

forecastProp <- ForecastArray
forecastProp[,,,1:3] <- ForecastArray[,,,1:3]/100

forecastPropFinal <- forecastProp
forecastPropFinal[forecastPropFinal==0.33] <- 1/3

# package up the final product and save
forecastList  <- list(lon=lon,
					  lat=lat,
					  forecast=forecastPropFinal,
					  mask=masterForecastMask,
					  timeMap=timeMap)


#save(forecastList, file='Data/RawProcessed/forecastObjIRI.Rda')