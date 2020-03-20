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

source('SharedCode/AllFunctions.R')
source('Forecast/Scratch/igFunctions.R')

# stuff from the namelist
impactsAnalysisName <- 'Cru_2.5_1951_1996_8_Data'

startYear <- 1997
endYear  <- 2016

ninoCutoff <- 1
ninaCutoff <- -1

FCAST_CATS <- 3

climatologyLength <- 30

dryCutoff <- 10

fancyDryMask <- FALSE

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

# load in the iri seasonal forecasts
load('Data/RawProcessed/forecastObjIRI.Rda')


# load in the appropriate observational data
load('Data/Forecast/CMAP_Data/obsTercile.Rda')
load('Data/Forecast/CMAP_Data/dryMask.Rda')
load('Data/RawProcessed/landProp_2.5.Rda')

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

	source('Forecast/Scratch/makeRealisticForecast.R')
	ensoRealForecastListLead <- makeRealisticForecast(leads[i],ensoForecastList,
												 ninaSeasonal,resultsList,forecastList)

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






# get the theoretical limit for the enso forecast over the same time period!
ddir <- 'Data/Forecast/CMAP_Data'
load(sprintf('%s/ensoForecast.Rda',ddir))

inds <- which(ensoProbForecastList$timeMap[,3] %in% ensoRealForecastListLead$timeMap[,3])

ensoProbForecastListTrim <- list(lon=ensoProbForecastList$lon,
								lat=ensoProbForecastList$lat,
								forecast=ensoProbForecastList$forecast[,,inds,],
								timeMap=ensoProbForecastList$timeMap[inds,])

ensoProbTotalGlobal  <- rpssTotalDecomp(ensoProbForecastListTrim,observedTercile,dryMask,tropics=FALSE)
ensoProbTotalTropics <- rpssTotalDecomp(ensoProbForecastListTrim,observedTercile,dryMask,tropics=TRUE)

textSize <- 1.5

# plot up the results
pdf('Figures/Results/LeadWork/leadTimeComp.pdf',14,7)
par(mfrow=c(1,2),mar=c(5, 5, 4, 3) + 0.1)

plot(1:4,iriTropics[,1],ylim=c(-0.01,0.04),type='b',lwd=2,col='blue',
	xlab='Lead Time (months)', ylab='Ranked Probability Skill Score',
	main='Forecast Skill (All Seasons 2003-2016)',
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

grid(lwd=1.5)
points(1:4,ensoTropics[,1],type='b',lwd=2,col='red')

points(1:4,iriGlobal[,1],type='b',lwd=2,col='blue',lty=3)
points(1:4,ensoGlobal[,1],type='b',lwd=2,col='red',lty=3)
abline(h=0)

legend(0.67,0.043, '(a)', bty='n',cex=textSize)

# legend('topright', c('IRI Forecast (Tropics)', 'IRI Forecast (Global)',
# 					 'ENSO Forecast (Tropics)', 'ENSO Forecast (Global)'),
# 			col=c('blue','blue', 'red', 'red'),lwd=2,lty=c(1,3,1,3),cex=textSize)


plot(1:4,iriTropics[,2],ylim=c(-0.01,0.04),type='b',lwd=2,col='blue',
	xlab='Lead Time (months)', ylab='RPS Resolution Score',
	main='Forecast Resolution (All Seasons 2003-2016)',
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

grid(lwd=1.5)
points(1:4,ensoTropics[,2],type='b',lwd=2,col='red')

points(1:4,iriGlobal[,2],type='b',lwd=2,col='blue',lty=3)
points(1:4,ensoGlobal[,2],type='b',lwd=2,col='red',lty=3)
abline(h=0)

legend('topright', c('IRI Forecast (Tropics)', 'IRI Forecast (Global)',
					 'ENSO Forecast (Tropics)', 'ENSO Forecast (Global)'),
			col=c('blue','blue', 'red', 'red'),lwd=2,lty=c(1,3,1,3),cex=textSize,bty='n')

legend(0.67,0.043, '(b)', bty='n',cex=textSize)

dev.off()

# ig plot
pdf('Figures/Results/LeadWork/leadTimeCompIgn.pdf',14,7)
par(mfrow=c(1,2))
plot(1:4,iriTropicsIgn[,1],ylim=c(-0.01,0.04),type='b',lwd=2,col='blue',
	xlab='Lead Time (months)', ylab='Effective Interest Rate',
	main='Forecast Skill (All Seasons 2003-2016)')
grid(lwd=1.5)
points(1:4,ensoTropicsIgn[,1],type='b',lwd=2,col='red')

points(1:4,iriGlobalIgn[,1],type='b',lwd=2,col='blue',lty=3)
points(1:4,ensoGlobalIgn[,1],type='b',lwd=2,col='red',lty=3)
abline(h=0)

# abline(h = ensoProbTotalTropics$rpss, col= 'red', lwd=1.5, lty=1)
# abline(h = ensoProbTotalGlobal$rpss, col= 'red', lwd=1.5, lty=3)

legend('topright', c('IRI Forecast (Tropics)', 'IRI Forecast (Global)',
					 'ENSO Forecast (Tropics)', 'ENSO Forecast (Global)'),
			col=c('blue','blue', 'red', 'red'),lwd=2,lty=c(1,3,1,3))


plot(1:4,iriTropicsIgn[,2],ylim=c(-0.01,0.04),type='b',lwd=2,col='blue',
	xlab='Lead Time (months)', ylab='IGN Resolution Score',
	main='Forecast Resolution (All Seasons 2003-2016)')
grid(lwd=1.5)
points(1:4,ensoTropicsIgn[,2],type='b',lwd=2,col='red')

points(1:4,iriGlobalIgn[,2],type='b',lwd=2,col='blue',lty=3)
points(1:4,ensoGlobalIgn[,2],type='b',lwd=2,col='red',lty=3)
abline(h=0)

legend('topright', c('IRI Forecast (Tropics)', 'IRI Forecast (Global)',
					 'ENSO Forecast (Tropics)', 'ENSO Forecast (Global)'),
			col=c('blue','blue', 'red', 'red'),lwd=2,lty=c(1,3,1,3))

dev.off()







