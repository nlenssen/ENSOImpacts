source('Forecast/Namelists/namelist_CMAP.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')

plotdir <- '/Users/lenssen/Desktop/scratchPlots/seasonalRpss'
# make a forecast mask
fcastCounts <- apply(iriForecastList$forecast[,,,1],c(1,2), function(x) sum(!is.na(x)))
fcastMask <- ifelse(fcastCounts < 50, NA, 1)

# rpss calculations on the annual mean (three figures in Fig 11)
field   <- rpssDecompField(iriForecastList,observedTercile)
ensoProbField   <- rpssDecompField(ensoProbForecastList,observedTercile)
rpssComp <- rpss(iriForecastList,observedTercile,ensoProbForecastList)

# get grid
lon  <- field$lon
lat  <- field$lat
nlon <- length(lon)
nlat <- length(lat)

# do the seasonal calculations

seasons <- c(2,5,8,11)

seasonPM <- c(-1,1)

# output objects
iriSeasonalField  <- array(NA, dim=c(nlon,nlat,length(seasons)))
ensoSeasonalField <- array(NA, dim=c(nlon,nlat,length(seasons)))
rpssSeasonalComp  <- array(NA, dim=c(nlon,nlat,length(seasons)))
rpssSeasonalComp2 <- array(NA, dim=c(nlon,nlat,length(seasons)))


for(i in 1:length(seasons)){
	tempSeason <- sort(unique(c(seasons[i],seasonPM + seasons[i])))


	# subset IRI
	iriInds <- which(iriForecastList$timeMap[,2] %in% tempSeason)

	subIriFcast <- iriForecastList
	subIriFcast$timeMap <- iriForecastList$timeMap[iriInds,]
	subIriFcast$forecast <- iriForecastList$forecast[,,iriInds,]

	# calculate iri field
	iriSeasonalField[,,i] <- rpss(subIriFcast,observedTercile)$global$field

	# subset ensoprob
	ensoIndsFull <- which(ensoProbForecastList$timeMap[,2] %in% tempSeason)
	ensoInds     <- ensoIndsFull[ensoIndsFull < 199]

	subEnsoFcast <- ensoProbForecastList
	subEnsoFcast$timeMap <- ensoProbForecastList$timeMap[ensoInds,]
	subEnsoFcast$forecast <- ensoProbForecastList$forecast[,,ensoInds,]

	# calculate ENSO field
	ensoSeasonalField[,,i] <- rpss(subEnsoFcast,observedTercile)$global$field


	# calculate the rpss relative to the other forecast
	# Q: WHY ARE THESE TWO CALCULATIONS DIFFERENT?!?
	#rpssSeasonalComp[,,i] <- rpss(subIriFcast,observedTercile,ensoProbForecastList)$global$field
	rpssSeasonalComp2[,,i] <- rpss(subIriFcast,observedTercile,subEnsoFcast)$global$field
}





# make some plots that look like the paper plots currently do

titleVec <- c('DJF - FMA', 'MAM - MJJ', 'JJA - ASO', 'SON - NDJ')
pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))
textSize <- 1.5

# IRI RPSS
zMax <- max(abs(iriSeasonalField),na.rm=T)
zForce <- 0.4
iriSeasonalField[iriSeasonalField > zForce] <- zForce
iriSeasonalField[iriSeasonalField < 0] <- 0
zr <- c(0,zForce)

pdf(sprintf('%s/iriSeasonal.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,iriSeasonalField[,,i] * fcastMask,zlim=zr, col=pal, 
		main=paste0('IRI (Clim. Baseline) [', titleVec[i],']'),xlab='',ylab='',ylim=c(-60,90))
	world(add=T)
}
dev.off()


# ENSO Prob RPSS
zMax <- max(abs(ensoSeasonalField),na.rm=T)
zForce <- 0.4
ensoSeasonalField[ensoSeasonalField > zForce] <- zForce
ensoSeasonalField[ensoSeasonalField < 0] <- 0
zr <- c(0,zForce)

pdf(sprintf('%s/ensoSeasonal.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,ensoSeasonalField[,,i] * fcastMask,zlim=zr, col=pal, 
		main=paste0('ENSO Prob. [', titleVec[i],']'),xlab='',ylab='',ylim=c(-60,90))
	world(add=T)
}
dev.off()


# Alt Baseline RPSS
zMax <- max(abs(rpssSeasonalComp),na.rm=T)
zForce <- 0.4
rpssSeasonalComp[rpssSeasonalComp > zForce] <- zForce
zr <- c(-zForce,zForce)

pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

pdf(sprintf('%s/iriCompSeasonal.pdf',plotdir),18,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,rpssSeasonalComp2[,,i] * fcastMask,zlim=zr, col=pal2, 
		main=paste0('IRI (ENSO Baseline) [', titleVec[i], ']'),xlab='',ylab='',ylim=c(-60,90),
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)
	world(add=T)
	abline(h=c(30,-30), lty=3)
	legend(-189, 92, paste0('(',letters[i],')'), cex = textSize, bty='n')
}
dev.off()


# Other plots as a sanity check

titleVec <- c('DJF - FMA', 'MAM - MJJ', 'JJA - ASO', 'SON - NDJ')

# IRI RPSS
zMax <- max(abs(iriSeasonalField),na.rm=T)
zForce <- 0.4
iriSeasonalField[iriSeasonalField > zForce] <- zForce
zr <- c(-zForce,zForce)

pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,iriSeasonalField[,,i],zlim=zr, col=pal2, 
		main=titleVec[i],xlab='',ylab='')
	world(add=T)
}



# ENSO Prob RPSS
zMax <- max(abs(ensoSeasonalField),na.rm=T)
zForce <- 0.4
ensoSeasonalField[ensoSeasonalField > zForce] <- zForce
zr <- c(-zForce,zForce)

pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,ensoSeasonalField[,,i],zlim=zr, col=pal2, 
		main=titleVec[i],xlab='',ylab='')
	world(add=T)
}


# Alt Baseline RPSS
zMax <- max(abs(rpssSeasonalComp),na.rm=T)
zForce <- 0.2
rpssSeasonalComp[rpssSeasonalComp > zForce] <- zForce
zr <- c(-zForce,zForce)

pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,rpssSeasonalComp2[,,i],zlim=zr, col=pal2, 
		main=titleVec[i],xlab='',ylab='')
	world(add=T)
}




