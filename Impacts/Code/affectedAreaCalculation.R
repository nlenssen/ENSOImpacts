load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

percentageArr      <- array(NA, dim=c(12,2,2,length(alphaLevels)))
totalPercentageArr <- array(NA, dim=c(2,2,length(alphaLevels)))

for(l in 1:length(alphaLevels)){
# unpack the list (to facilitate multiple data formats)
lon <- resultsList$lon
lat <- resultsList$lat
tropicsLat <- which(lat > -30 & lat < 30)

timeMap <- resultsList$timeMap
sampleSize <- resultsList$sampleSize
pValArray <- resultsList$pValArray
dryMask   <- resultsList$dryMask

# carry out the analysis
landMask <- ifelse(is.na(sampleSize[,,1]),NA,1)

pValArray[pValArray == -1 ] <- NA

gl <- make.surface.grid(list(lon,lat))


totalArea <- rep(NA,4)
booleanDryMask <- array(NA, dim=dim(dryMask))

for(i in 1:12){
	dryInds <- which(c(dryMask[,,i])==1)
	
	tempMat <- landMask
	tempMat[dryInds] <- 0

	booleanDryMask[,,i] <- ifelse(tempMat==1,TRUE,FALSE)

	totalArea[i] <- sum(c(tempMat)*cos(gl[,2]*(pi/180)),na.rm=T)
}


# calculate max extent of non-dry reigons
maxExtent <- apply(booleanDryMask,c(1,2),any)
maxArea   <- sum(maxExtent*cos(gl[,2]*(pi/180)),na.rm=T)

# make the a boolean effected matrix
effectedArray <- pValArray
effectedArray[effectedArray==-1] <- NA
effectedArray <- ifelse(!is.na(effectedArray) & effectedArray > alphaLevels[l], 1,NA)


# calculate the percentage effected
naArray <- array(NA, dim=c(dim(dryMask),2,2))

percentage <- array(NA, dim=c(12,2,2))
totalPercentage <- matrix(NA, nrow=2,ncol=2)

for(el in 1:2){
	for(lh in 1:2){
		tempYearMat <- matrix(NA, nrow=length(lon),ncol=length(lat))

		for(s in 1:12){
			
			percentage[s,el,lh] <- sum(effectedArray[,,s,el,lh]*cos(gl[,2]*(pi/180)),na.rm=T)/totalArea[s]
			tempYearMat[!is.na(effectedArray[,,s,el,lh])] <- 1
		}

		yearArea <- sum(c(tempYearMat)*cos(gl[,2]*(pi/180)),na.rm=T)
		
		totalPercentage[el,lh] <- yearArea/maxArea
	}
}

percentageArr[,,,l]     <- percentage
totalPercentageArr[,,l] <- totalPercentage
}

save(percentageArr,totalPercentageArr,file=sprintf('%s/affectedArea.Rda',ddir))