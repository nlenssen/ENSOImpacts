# The goal of this script is to run the forecast analysis for leads of 1-4
# months over the record of the ENSO forecast. The approximate workflow
# should be
#
# (1) Load and subset the IRI forecast and the ENSO forecast at lead i
# (2) Create the appropriate realistic ENSO forecast at lead i
# (3) Verify the two forecasts at lead i
#		- For now, just calculate the RPSS of the total forecast and compare
#		over the record
#		- Allow easy drop of more stats as needed.

source('Forecast/Code/igFunctions.R')

# loop over the various leads
leads <- 1:4

# loop objs
iriGlobal <- matrix(NA,nrow=length(leads),ncol=3)
iriTropics <- matrix(NA,nrow=length(leads),ncol=3)
ensoGlobal <- matrix(NA,nrow=length(leads),ncol=3)
ensoTropics <- matrix(NA,nrow=length(leads),ncol=3)

iriGlobalIgn <- matrix(NA,nrow=length(leads),ncol=2)
iriTropicsIgn <- matrix(NA,nrow=length(leads),ncol=2)
ensoGlobalIgn <- matrix(NA,nrow=length(leads),ncol=2)
ensoTropicsIgn <- matrix(NA,nrow=length(leads),ncol=2)
##################
# (1) load in shit
##################

# load the nino 3.4 time series
load('Data/RawProcessed/ninaSeasonal.Rda')
load('Data/Impacts/Cru_0.5_1951_2016_Data/ninaSeasonalInd.Rda')


# load in the iri seasonal forecasts
load('Data/RawProcessed/forecastObjIRI.Rda')
load('Data/RawProcessed/landProp_2.5.Rda')

# load in the appropriate observational data
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

landPropMat <- landPropList$prop

# load the count array to get spatial signal
# reminder that countArray has dimensions
# lon x lat x season x el/la x high/low anom
load(sprintf('Data/Impacts/%s/maskedAnalysisResults.Rda',impactsAnalysisName))


for(i in 1:length(leads)){


	# load in the ENSO forecast for lead
	load(sprintf('Data/RawProcessed/ensoForecastIRILead%d.Rda',leads[i]))

	####################################
	# (2) create realistic ENSO forecast 
	####################################

	# build ENSO real forecast (function)
	ensoRealForecastListLead <- makeRealisticForecast(leads[i],ensoForecastList,
												 ninaIndicator,ninaSeasonal,resultsList,forecastList)

	# subset the iri forecast to the same time with the correct lead
	inds <- which(forecastList$timeMap[,3] %in% ensoRealForecastListLead$timeMap[,3] & 
				  forecastList$timeMap[,4] == leads[i])

	iriForecastList <- list(lon=forecastList$lon,
							lat=forecastList$lat,
							forecast=forecastList$forecast[,,inds,],
							mask=forecastList$mask,
							timeMap=forecastList$timeMap[inds,])


	# make the missing values NA so they are prooerly handled during the ign analysis
	fcastClean <- iriForecastList$forecast[,,,1:3]
	fcastClean[fcastClean < 0 ] <- NA

	iriForecastListClean <- list(lon=forecastList$lon,
							lat=forecastList$lat,
							forecast=fcastClean,
							mask=forecastList$mask,
							timeMap=forecastList$timeMap[inds,])

	#######################################
	# (3) verify the ENSO and IRI forecasts
	#######################################

	# RPS Decomposition
	iriTotalGlobal  <- rpssTotalDecomp(iriForecastList,observedTercile,dryMask,tropics=FALSE)
	iriTotalTropics <- rpssTotalDecomp(iriForecastList,observedTercile,dryMask,tropics=TRUE)

	ensoRealTotalGlobal  <- rpssTotalDecomp(ensoRealForecastListLead,observedTercile,dryMask,tropics=FALSE)
	ensoRealTotalTropics <- rpssTotalDecomp(ensoRealForecastListLead,observedTercile,dryMask,tropics=TRUE)


	iriGlobal[i,]  <- c(iriTotalGlobal$rpss,iriTotalGlobal$res,iriTotalGlobal$rel)
	iriTropics[i,] <- c(iriTotalTropics$rpss,iriTotalTropics$res,iriTotalTropics$rel)

	ensoGlobal[i,]  <- c(ensoRealTotalGlobal$rpss,ensoRealTotalGlobal$res,ensoRealTotalGlobal$rel)
	ensoTropics[i,] <- c(ensoRealTotalTropics$rpss,ensoRealTotalTropics$res,ensoRealTotalTropics$rel)


	# RIGN Decomposition
	yzIri  <- yzForecast(iriForecastListClean, observedTercile)
	yzEnso <- yzForecast(ensoRealForecastListLead, observedTercile)

	iriTotalGlobalEir     <-  rignTotal(yzIri, landPropMat, tropics = FALSE)$eir
	iriTotalGlobalIgnRes  <-  rignTotalDecomp(yzIri, landPropMat, tropics = FALSE)$res

	iriTotalTropicsEir    <- rignTotal(yzIri, landPropMat, tropics = TRUE)$eir
	iriTotalTropicsIgnRes <- rignTotalDecomp(yzIri, landPropMat, tropics = TRUE)$res

	ensoRealTotalGlobalEir     <- rignTotal(yzEnso, landPropMat, tropics = FALSE)$eir
	ensoRealTotalGlobalIgnRes  <- rignTotalDecomp(yzEnso, landPropMat, tropics = FALSE)$res

	ensoRealTotalTropicsEir    <- rignTotal(yzEnso, landPropMat, tropics = TRUE)$eir
	ensoRealTotalTropicsIgnRes <- rignTotalDecomp(yzEnso, landPropMat, tropics = TRUE)$res


	iriGlobalIgn[i,]  <- c(iriTotalGlobalEir,iriTotalGlobalIgnRes)
	iriTropicsIgn[i,] <- c(iriTotalTropicsEir,iriTotalTropicsIgnRes)

	ensoGlobalIgn[i,]  <- c(ensoRealTotalGlobalEir,ensoRealTotalGlobalIgnRes)
	ensoTropicsIgn[i,] <- c(ensoRealTotalTropicsEir,ensoRealTotalTropicsIgnRes)
}


save(iriGlobal,iriTropics,ensoGlobal,ensoTropics,file=sprintf('%s/multiLead.Rda',ddir))
