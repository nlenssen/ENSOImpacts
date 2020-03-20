

# eventual parameters to add to some sort of namelist
impactsFDR <- 0.35

fdrDryMask <- TRUE

# load in the results. The object with the p-values are 
# provided in resultsList$pValArray
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

# pull out some useful objects
lon <- resultsList$lon
lat <- resultsList$lat

# for area calculations
landMask <- ifelse(is.na(resultsList$sampleSize[,,1]),NA,1)



# fdr adjusted fields to be stored here
fdrArray <- array(NA, dim=c(length(lon),length(lat), 12,2,2))
initialArea <- array(NA, dim=c(12,2,2))
fdrArea <- array(NA, dim=c(12,2,2))

# Eventual loop over all fields!

for(mn in 1:12){
for(el in 1:2){
for(ab in 1:2){


	tempFieldRaw <- resultsList$pValArrayFull[,,mn,el,ab]

	# figure out which zeros are due to low data
	cutoffValue <- 20
	lowDataMask <- ifelse(resultsList$sampleSize[,,mn] < cutoffValue,NA,1)
	

	# then mask for dry areas
	if(fdrDryMask){
		tempDryMask <- ifelse(!is.na(resultsList$dryMask[,,mn]),NA,1)
	} else{
		tempDryMask <- matrix(1,length(lon),length(lat))
	}

	# then take the resulting field to pass off to the FDR
	tempField <- tempFieldRaw * lowDataMask * tempDryMask

	# Now, tempField is the p-values of our multiple testing problem. Our first
	# attempt at handling multiple testing is throwing off-the-shelf BH-procedure
	# at is using the built-in R function

	# fdr calculation call
	fdrArray[,,mn,el,ab] <- fdrFunction(tempField, FDR = impactsFDR,
				   plotName = NULL)

	# Do some area calculations

	# build the base avaiable land area
	tempDryMask <- ifelse(!is.na(resultsList$dryMask[,,mn]),NA,1)

	cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))* tempDryMask * landMask

	base85 <- ifelse(tempField > 0.85, 1, NA)

	initialArea[mn,el,ab] <- sum(base85 * cosMat, na.rm=T) / sum(cosMat,na.rm=T)

	fdrArea[mn,el,ab] <- sum(fdrArray[,,mn,el,ab] * cosMat, na.rm=T)/sum(cosMat,na.rm=T)

}
}
}



# make plots to compare the current paper results with the FDR results
analysisList <- c('Data/Impacts/Cru_0.5_1951_2016_11_Data',
				  'Data/Impacts/Cru_2.5_1951_2016_11_Data',
				  'Data/Impacts/New_0.5_1951_1996_8_Data')

# corresponding to 85, 90, 95, 99
alphaInd <- 1

# first, make an affected area series in time
percentageOrg <- array(NA, dim=c(12,2,2,length(analysisList)))

for(i in 1:length(analysisList)){

	load(sprintf('%s/affectedArea.Rda',analysisList[i]))
	load(sprintf('%s/ninaSeasonalInd.Rda',analysisList[i]))

	# get at the 85% confidence level
	percentageOrg[,,,i] <- percentageArr[,,,alphaInd]
}


# Make line plots for four distinct seasons from the new analysis
textSize <- 1.5

seasonInds <- c(1,4,7,10)
seasonNames <- c('DJF', 'MAM', 'JJA', 'SON')
ninoNames <- c('EN','LN')
ninoLongNames <- c('El Niño', 'La Niña')


anomNames <- c('W','D')
anomLongNames <- c('Above Average', 'Below Average')

subInds <- c(1,4,7,10)


pdf('~/Desktop/scratchPlots/fdrWork/areaViz.pdf',14,7)
par(mfrow=c(1,2),mar=c(5, 5, 4, 3) + 0.1)
plot(1:4,NULL,xaxt = "n", xlab='', ylab='Global Proportion Affected',
	ylim=c(0,0.5), main=sprintf('El Niño'),
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)

axis(1, at=1:4, labels=seasonNames,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

points(1:4,initialArea[subInds,1,1],
	type='l',lwd=2,col='deepskyblue3')


points(1:4,initialArea[subInds,1,2],
	type='l',lwd=2,col='salmon4',lty=1)

points(1:4,fdrArea[subInds,1,1],
	type='l',lwd=2,col='deepskyblue3',lty=2)


points(1:4,fdrArea[subInds,1,2],
	type='l',lwd=2,col='salmon4',lty=2)

legend('topleft', c('Above Normal','Below Normal'),
col=c('deepskyblue3','salmon4'),lty=c(1,1),lwd=2.5,cex=1.5,horiz=F,bty='n')

plot(1:4,NULL,xaxt = "n", xlab='', ylab='',
ylim=c(0,0.5),
main=sprintf('La Niña'),
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

grid(lwd=1.5)

axis(1, at=1:4, labels=seasonNames,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

points(1:4,initialArea[subInds,2,1],
	type='l',lwd=2,col='deepskyblue3')


points(1:4,initialArea[subInds,2,2],
	type='l',lwd=2,col='salmon4',lty=1)

points(1:4,fdrArea[subInds,2,1],
	type='l',lwd=2,col='deepskyblue3',lty=2)


points(1:4,fdrArea[subInds,2,2],
	type='l',lwd=2,col='salmon4',lty=2)

dev.off()
# some strange warnings, investigate why all masked in certain fields
sigCounts <- apply(fdrArray,c(3,4,5),function(x) sum(!is.na(x)))

# [,,4,1,1] is a field with all NA
mn <- 7
el <- 2
ab <- 2

tempFieldRaw <- resultsList$pValArrayFull[,,mn,el,ab]

# figure out which zeros are due to low data
cutoffValue <- 20
lowDataMask <- ifelse(resultsList$sampleSize[,,mn] < cutoffValue,NA,1)

# then mask for dry areas
tempDryMask <- ifelse(!is.na(resultsList$dryMask[,,mn]),NA,1)
tempDryMask <- matrix(1,length(lon),length(lat))

# then take the resulting field to pass off to the FDR
tempField <- tempFieldRaw * lowDataMask * tempDryMask



result <- fdrFunction(tempField, FDR = impactsFDR,
				   plotName = '~/Desktop/scratchPlots/fdrWork/strange.pdf')







# # calculate the proportion of area from the 85% significance level with no correction
# base85 <- ifelse(tempField > 0.85, 1, NA)

# cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))

# initialArea <- sum(base85 * cosMat, na.rm=T) / sum(cosMat)

# area255 <- sum(fdr255 * cosMat, na.rm=T)/sum(cosMat) / initialArea
# area30  <- sum(fdr30 * cosMat, na.rm=T)/sum(cosMat) / initialArea
# area35  <- sum(fdr35 * cosMat, na.rm=T)/sum(cosMat) / initialArea
# area40  <- sum(fdr40 * cosMat, na.rm=T)/sum(cosMat) / initialArea
# area50  <- sum(fdr50 * cosMat, na.rm=T)/sum(cosMat) / initialArea
