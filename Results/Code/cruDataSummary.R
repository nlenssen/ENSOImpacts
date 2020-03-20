source('Results/Namelists/namelist_FinalResults.Rnl')

load('Data/RawProcessed/landProp_0.5.Rda')

###############################################################################
# Load in the relevant count Data
###############################################################################

# load in the ncdf metadata
handle <- nc_open(sprintf('Data/Raw/%s',cruFile))

lon <- ncvar_get(handle,'lon')
lat <- ncvar_get(handle,'lat')

# set up all of the necessary time objects to properly subset the data
time <- ncvar_get(handle,'time')

tFull <- seq(1901,cruEndYear+1,length=length(time)+1)[1:length(time)]
year <- floor(tFull)
month <- rep(1:12,length(time)/12)

fullTimeMat <- cbind(year,month)

# subset based on starting year and make all of the necessary time objects
timeInds <- which(year >= cruStartYear)

timeMat <- fullTimeMat[timeInds,]
continuousTime <- tFull[timeInds]
# grab the ppt monthly time series for the single location

counts <- ncvar_get(handle, 'stn', start=c(1,1,timeInds[1]),
								count=c(-1,-1,length(timeInds)) )

nc_close(handle)


# Fill array with NAs

counts[counts==-999] <- NA

# make tYear
tYear <- (cruStartYear):cruEndYear

############################
# Convert counts to coverage
############################
coverageArray <- counts
coverageArray[coverageArray<3] <- NA
coverageArray[coverageArray>2] <- 1

noCoverageArray <- counts
noCoverageArray[noCoverageArray>2] <- NA
noCoverageArray[noCoverageArray<3] <- 1

####################################
# Get time series of global coverage
####################################
gl <- make.surface.grid(list(lon,lat))

landLandProp <- c(landPropList$prop)*cos(gl[,2]*(pi/180))

totalArea <- sum(landLandProp,na.rm=T)

# same caclulation for 30/30
inds50 <- which(lat > -50 & lat < 50)
lat50 <- lat[inds50]
gl50 <- make.surface.grid(list(lon,lat[inds50]))

landLandProp50 <- c(landPropList$prop[,inds50])*cos(gl50[,2]*(pi/180))

totalArea50 <- sum(landLandProp50,na.rm=TRUE)

# same caclulation for 30/30
inds30 <- which(lat > -30 & lat < 30)
lat30 <- lat[inds30]
gl30 <- make.surface.grid(list(lon,lat[inds30]))

landLandProp30 <- c(landPropList$prop[,inds30])*cos(gl30[,2]*(pi/180))

totalArea30 <- sum(landLandProp30,na.rm=T)

# calculate the global coverage time series
propCovered <- rep(NA, dim(counts)[3])
propCovered50 <- rep(NA, dim(counts)[3])
propCovered30 <- rep(NA, dim(counts)[3])

for(i in 1:length(propCovered)){
	tempSlice <- coverageArray[,,i]

	tempArea  <- sum(c(tempSlice)*landLandProp,na.rm=T)
	propCovered[i] <- tempArea/totalArea

	tempArea50 <- sum(c(tempSlice[,inds50])*landLandProp50,na.rm=T)
	propCovered50[i] <- tempArea50/totalArea50

	tempArea30 <- sum(c(tempSlice[,inds30])*landLandProp30,na.rm=T)
	propCovered30[i] <- tempArea30/totalArea30
}


# get every 12th to remove the seasonal cycle




####################################
# Save the output needed for plots
####################################
tSeason <- rep(1:12, length(propCovered)/12)
plotInds <- which(tSeason==3)

counts90 <- counts[,,plotInds[which(tYear==1990)]]
counts15 <- counts[,,plotInds[which(tYear==2015)]]

covered90 <- ifelse(counts90 >2,1,NA)
notCovered90 <- ifelse(counts90 <3,1,NA)
covered15 <- ifelse(counts15 >2,1,NA)
notCovered15 <- ifelse(counts15 <3 &covered90==1,1,NA)




tYearCru <- tYear
cruLon <- lon
cruLat <- lat

# save the output so the figure script can lift less
save(cruLon,cruLat,tYearCru,
	propCovered,propCovered30,propCovered50,
	covered90,covered15,notCovered90,notCovered15,
	file='Data/RawProcessed/cruDataSummary.Rda')



