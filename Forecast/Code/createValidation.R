load(sprintf('Data/RawProcessed/%s',validationDataName))


# create the climatology for each of the 12 months. Get both the means
# and the tercile values
tercileList <- list()
meanList    <- list()
dryMask <- array(NA, dim=c(nrow(observationList$observations),ncol(observationList$observations),12))
# 
climatologyYears <- 
	climatologyStartYear:(climatologyStartYear+climatologyLength-1)

for(i in 1:12){
	timeInds <- which(observationList$timeMap[,2]==i & 
					  observationList$timeMap[,1] %in% climatologyYears)
	
	tempTercile <- aperm(apply(observationList$observations[,,timeInds],
								c(1,2),quantile,probs=c(1/3,2/3),na.rm=T),
							c(2,3,1))

	tercileList[[i]] <- tempTercile

	meanList[[i]] <- apply(observationList$observations[,,timeInds],
					 	c(1,2),mean,na.rm=T)
	
}



# make the dry mask two different ways that appears to affect the results significatly

# first get the annual average mean (over climatology period)
if(fancyDryMask){
	annualTotalsArray <-  array(NA, dim=c(nrow(observationList$observations),
										  ncol(observationList$observations),
										  length(climatologyYears)))

	for(i in 1:length(climatologyYears)){
		# need to take non-overlapping seasons that sum up to a calendar year!
		yearInds <- which(observationList$timeMap[,1]==climatologyYears[i] &
						  observationList$timeMap[,2] %in% c(2,5,8,11) )

		annualTotalsArray[,,i] <- apply(observationList$observations[,,yearInds],c(1,2),sum)
	}

	meanAnnualTotal <- apply(annualTotalsArray,c(1,2),mean,na.rm=T)



	for(s in 1:12){
	for(i in 1:nrow(observationList$observations)){
	for(j in 1:ncol(observationList$observations)){

		isDry <- FALSE

		if(is.na(meanList[[s]][i,j])| is.na(meanAnnualTotal[i,j])) next

		if(meanList[[s]][i,j] <= 0.15 * meanAnnualTotal[i,j] & meanList[[s]][i,j] <= 50){
			isDry <- TRUE
		}else if(tercileList[[s]][i,j,1] < dryCutoff){
			isDry <- TRUE
		}

		dryMask[i,j,s] <- ifelse(isDry, 1, 0)
	}
	}
	}

} else{
	for(i in 1:12) dryMask[,,i] <- ifelse(tercileList[[i]][,,2] < dryCutoff,0,1)
}


# figure out which of the terciles was observed for the record
obsTercile <- array(0, dim=dim(observationList$observations))
obsTercile[is.na(observationList$observations[,,1])] <- NA

for(t in 1:dim(observationList$observations)[3]){
	seasonInd <- observationList$timeMap[t,2]

	# pull the correct climatology
	tempTercile <- tercileList[[seasonInd]]

	# Use (-1,0,1) as the indicators for the low avg and high
	obsTercile[,,t][observationList$observations[,,t] < tempTercile[,,1]] <- -1
	obsTercile[,,t][observationList$observations[,,t] > tempTercile[,,2]] <-  1
}

climatology <- list(lon=observationList$lon,lat=observationList$lat,
	mean = meanList, tercile=tercileList,timeMap=observationList$timeMap)

observedTercile <- list(lon=observationList$lon,lat=observationList$lat,
	tercile=obsTercile,timeMap=observationList$timeMap)

save(climatology,file=sprintf('%s/climatology.Rda',ddir))
save(observedTercile,file=sprintf('%s/obsTercile.Rda',ddir))
save(dryMask,file=sprintf('%s/dryMask.Rda',ddir))

