# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# load in the data and do initial cleaning
rawLand <- t(as.matrix(read.table(filenameLand,skip=6,header=F,sep=" ")))
rawWater <- t(as.matrix(read.table(filenameWater,skip=6,header=F,sep=" ")))

rawLand[rawLand==-9999] <- 0
rawWater[rawWater==-9999] <- 0

totalLandRaw <- (rawLand + rawWater)[1:720,]

# rotate around and make the world arranged correctly
nlon <- nrow(totalLandRaw)
nlat <- ncol(totalLandRaw)

lon <- seq(-179.75,179.75,by=0.5)
lat <- seq(-59.75,84.75,by=0.5)

landRawFlip <- matrix(NA, nrow=nlon,ncol=nlat)

for(i in 1:nlat){
	landRawFlip[,i] <- totalLandRaw[,nlat-i+1]
}

landAreaFinal <- ifelse(landRawFlip==0,NA,landRawFlip)

# now, calculate the area of the gridcells depending on latitude
borderLat <- seq(-60,85,by=0.5) *  pi/180

R <- 6378.1
deltaLon <- 0.5

latArea <- rep(NA,nlat)

for(i in 1:(length(borderLat)-1)){
	latArea[i] <- (pi/180)*R^2*abs(sin(borderLat[i])-sin(borderLat[i+1])) * deltaLon
}

# now, apply the lat area to the matrix to determine proportion
landProp <- matrix(NA, nrow=nlon, ncol=nlat)
for(i in 1:nrow(landAreaFinal)){
	landProp[i,] <- landAreaFinal[i,] / latArea
}

# make all of the values that are close to one exactly one
landPropCorrected <- ifelse(landProp > (1 - tol),1,landProp)

# finally, add in the last latitudes to have the same dimensionality as the
# other global datasets
latFull <- seq(-89.75,89.75,by=0.5)
nlatFull <- length(latFull)

latInds <- which(latFull %in% lat)

landPropFinal <- matrix(NA, nrow=nlon,ncol=nlatFull)

landPropFinal[,latInds] <- landPropCorrected

landPropFinal[is.na(landPropFinal)] <- 0
# build the final list object
landPropList <- list(lon=lon,lat=latFull,prop=landPropFinal)

save(landPropList,file='Data/RawProcessed/landProp_0.5.Rda')

