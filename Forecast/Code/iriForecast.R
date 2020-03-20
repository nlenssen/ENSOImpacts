load('Data/RawProcessed/forecastObjIRI.Rda')

# rebuild the iri forecast list with the desired lead time
leadIndsFull <- which(forecastList$timeMap[,4]==lead & forecastList$timeMap[,1] <= endYear)

# handle the inability of the impacts analysis to get the last two months of the final
# year of data
badLeadInds <- which(forecastList$timeMap[leadIndsFull,1] == endYear & 
	forecastList$timeMap[leadIndsFull,2] >10)
leadInds <- leadIndsFull[-badLeadInds]

# make the correct regions in the forecast NA
forecastNA <- forecastList$forecast
forecastNA[forecastNA<0] <- NA

# mask all of the forecasts to make sure we aren't evaluating strange ocean reigons
# due to old forecast maps that are hand drawn
forecastFinal <- forecastNA[,,leadInds,]

for(i in 1:length(leadInds)){
	for(k in 1:3){
		forecastFinal[,,i,k] <- forecastNA[,,leadInds[i],k] * forecastList$mask		
	}
}


# package up and save
iriForecastList <- list(lon=forecastList$lon,
					  lat=forecastList$lat,
					  forecast=forecastFinal,
					  timeMap=forecastList$timeMap[leadInds,])

save(iriForecastList,file=sprintf('%s/iriForecast.Rda',ddir))