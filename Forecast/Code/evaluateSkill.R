# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')

field   <- rpssDecompField(iriForecastList,observedTercile)
global  <- rpssDecompSeries(iriForecastList,observedTercile,dryMask,tropics=FALSE)
tropics <- rpssDecompSeries(iriForecastList,observedTercile,dryMask,tropics=TRUE)

totalGlobal  <- rpssTotalDecomp(iriForecastList,observedTercile,dryMask,tropics=FALSE)
totalTropics <- rpssTotalDecomp(iriForecastList,observedTercile,dryMask,tropics=TRUE)

iriReliability <- reliability(iriForecastList,observedTercile,binSize=0.05)

iriGrocTropics <- grocSeries(iriForecastList,observedTercile,tropics=TRUE)
iriGrocField <- grocField(iriForecastList,observedTercile,tropics=FALSE)

# calculation for the det enso forecast
ensoField   <- rpssDecompField(ensoForecastList,observedTercile)
ensoGlobal  <- rpssDecompSeries(ensoForecastList,observedTercile,dryMask,tropics=FALSE)
ensoTropics <- rpssDecompSeries(ensoForecastList,observedTercile,dryMask,tropics=TRUE)

ensoTotalGlobal <- rpssTotalDecomp(ensoForecastList,observedTercile,dryMask,tropics=FALSE)
ensoTotalTropics <- rpssTotalDecomp(ensoForecastList,observedTercile,dryMask,tropics=TRUE)

ensoReliability <- reliability(ensoForecastList,observedTercile,binSize=0.05)

ensoDetGrocTropics <- grocSeries(ensoForecastList,observedTercile,tropics=TRUE)

# calculation for the prob enso forecast
ensoProbField   <- rpssDecompField(ensoProbForecastList,observedTercile)
ensoProbGlobal  <- rpssDecompSeries(ensoProbForecastList,observedTercile,dryMask,tropics=FALSE)
ensoProbTropics <- rpssDecompSeries(ensoProbForecastList,observedTercile,dryMask,tropics=TRUE)

ensoProbTotalGlobal  <- rpssTotalDecomp(ensoProbForecastList,observedTercile,dryMask,tropics=FALSE)
ensoProbTotalTropics <- rpssTotalDecomp(ensoProbForecastList,observedTercile,dryMask,tropics=TRUE)

ensoProbReliability <- reliability(ensoProbForecastList,observedTercile,binSize=0.05)

ensoProbGrocTropics <- grocSeries(ensoProbForecastList,observedTercile,tropics=TRUE)
ensoProbGrocField <- grocField(ensoProbForecastList,observedTercile,tropics=FALSE)

# calculation for the real enso forecast
realInds <- ensoRealForecastList$timeInds

ensoRealForecastListTrim <- list(lon=ensoRealForecastList$lon,
								 lat=ensoRealForecastList$lat,
						 		 forecast=ensoRealForecastList$forecast[,,realInds,],
						 		 timeMap=ensoRealForecastList$timeMap[realInds,])

ensoRealField   <- rpssDecompField(ensoRealForecastList,observedTercile)
ensoRealGlobal  <- rpssDecompSeries(ensoRealForecastList,observedTercile,dryMask,tropics=FALSE)
ensoRealTropics <- rpssDecompSeries(ensoRealForecastList,observedTercile,dryMask,tropics=TRUE)

ensoRealTotalGlobal <- rpssTotalDecomp(ensoRealForecastList,observedTercile,dryMask,tropics=FALSE)
ensoRealTotalTropics <- rpssTotalDecomp(ensoRealForecastList,observedTercile,dryMask,tropics=TRUE)

ensoRealReliability <- reliability(ensoRealForecastList,observedTercile,binSize=0.05)

ensoRealGrocTropics <- grocSeries(ensoRealForecastListTrim,observedTercile,tropics=TRUE)


# alt reference forecast rpss calculations
rpssComp <- rpss(iriForecastList,observedTercile,ensoProbForecastList)


# collect the total calcs for easier plotting
totalTabGlobal <- rbind(
totalGlobal,
ensoProbTotalGlobal,
ensoRealTotalGlobal,
ensoTotalGlobal)

totalTabTropics <- rbind(
totalTropics,
ensoProbTotalTropics,
ensoRealTotalTropics,
ensoTotalTropics)

totalGrocTropics <- c(iriGrocTropics$total,ensoProbGrocTropics$total,ensoRealGrocTropics$total,ensoDetGrocTropics$total)


#rm the uncesseary large objects and save the entire image
rm(iriForecastList,ensoForecastList,ensoProbForecastList,observedTercile)
save.image(file=sprintf('%s/skillScores.Rda',ddir))