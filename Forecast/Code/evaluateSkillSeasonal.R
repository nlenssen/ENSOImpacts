# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')

# seasons to run the analysis on
sInds <- c(2,5,8,11)

seasonNames <- c('DJF', 'JFM', 'FMA', 
				 'MAM', 'AMJ', 'MJJ',
				 'JJA', 'JAS', 'ASO',
				 'SON', 'OND', 'NDJ')

# loop through the key parts of the analysis for each season
for(i in 1:length(sInds)){
###
# subset the data appropriately 
###

# do the observations
observedTercileSub <- observedTercile
obsInds <- which(observedTercileSub$timeMap[,2]==sInds[i])

observedTercileSub$timeMap <- observedTercile$timeMap[obsInds,]
observedTercileSub$tercile <- observedTercile$tercile[,,obsInds]

# do the forecasts
iriForecastListSub      <- subsetForecastList(iriForecastList,sInds[i])
ensoForecastListSub     <- subsetForecastList(ensoForecastList,sInds[i])
ensoProbForecastListSub <- subsetForecastList(ensoProbForecastList,sInds[i])
ensoRealForecastListSub <- subsetForecastList(ensoRealForecastList,sInds[i])


###
# Run the verification on the subset
###
field   <- rpssDecompField(iriForecastListSub,observedTercileSub)
global  <- rpssDecompSeries(iriForecastListSub,observedTercileSub,dryMask,tropics=FALSE)
tropics <- rpssDecompSeries(iriForecastListSub,observedTercileSub,dryMask,tropics=TRUE)

totalGlobal  <- rpssTotalDecomp(iriForecastListSub,observedTercileSub,dryMask,tropics=FALSE)
totalTropics <- rpssTotalDecomp(iriForecastListSub,observedTercileSub,dryMask,tropics=TRUE)

iriReliability <- reliability(iriForecastListSub,observedTercileSub,binSize=0.05)

# calculation for the det enso forecast
ensoField   <- rpssDecompField(ensoForecastListSub,observedTercileSub)
ensoGlobal  <- rpssDecompSeries(ensoForecastListSub,observedTercileSub,dryMask,tropics=FALSE)
ensoTropics <- rpssDecompSeries(ensoForecastListSub,observedTercileSub,dryMask,tropics=TRUE)

ensoTotalGlobal <- rpssTotalDecomp(ensoForecastListSub,observedTercileSub,dryMask,tropics=FALSE)
ensoTotalTropics <- rpssTotalDecomp(ensoForecastListSub,observedTercileSub,dryMask,tropics=TRUE)

ensoReliability <- reliability(ensoForecastListSub,observedTercileSub,binSize=0.05)


# calculation for the prob enso forecast
ensoProbField   <- rpssDecompField(ensoProbForecastListSub,observedTercileSub)
ensoProbGlobal  <- rpssDecompSeries(ensoProbForecastListSub,observedTercileSub,dryMask,tropics=FALSE)
ensoProbTropics <- rpssDecompSeries(ensoProbForecastListSub,observedTercileSub,dryMask,tropics=TRUE)

ensoProbTotalGlobal  <- rpssTotalDecomp(ensoProbForecastListSub,observedTercileSub,dryMask,tropics=FALSE)
ensoProbTotalTropics <- rpssTotalDecomp(ensoProbForecastListSub,observedTercileSub,dryMask,tropics=TRUE)

ensoProbReliability <- reliability(ensoProbForecastListSub,observedTercileSub,binSize=0.05)

# collect the total calcs for easier plotting
totalTabGlobal <- rbind(
totalGlobal,
ensoProbTotalGlobal,
ensoTotalGlobal)

totalTabTropics <- rbind(
totalTropics,
ensoProbTotalTropics,
ensoTotalTropics)


###
# Make some plots
###

lon <- iriForecastList$lon
lat <- iriForecastList$lat

textSize <- 1.25


# resolution fields for the two of interest
maxCutoff <- 0.25

z <- array(NA, dim=c(length(lon),length(lat),3))

z[,,1]  <- field$decomp73$res
z[,,2]  <- ensoProbField$decomp73$res
z[,,3]  <- ensoField$decomp73$res

z[z>maxCutoff] <- maxCutoff
z[z<0.01] <- 0

zr <- c(0,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

mainVec <- c(sprintf("IRI Forecast %s Resolution",seasonNames[sInds[i]]),
			 sprintf("ENSO Probabilistic Forecast %s Resolution",seasonNames[sInds[i]]))

pdf(sprintf('%s/Seasonal/fieldRes%02d.pdf',plotdir,sInds[i]),10,12)
par(mfrow=c(2,1))

image.plot(lon,lat,z[,,1],col=pal,ylim=c(-60,90),main=mainVec[1],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)

image.plot(lon,lat,z[,,2],col=pal,ylim=c(-60,90),main=mainVec[2],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)
dev.off()


}