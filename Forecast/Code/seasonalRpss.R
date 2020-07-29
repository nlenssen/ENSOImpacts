skipCalc <- FALSE

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')

# make a forecast mask
fcastCounts <- apply(iriForecastList$forecast[,,,1],c(1,2), function(x) sum(!is.na(x)))
fcastMask <- ifelse(fcastCounts < 50, NA, 1)

# get grid
lon  <- iriForecastList$lon
lat  <- iriForecastList$lat
nlon <- length(lon)
nlat <- length(lat)

# do the seasonal calculations

seasons <- c(2,5,8,11)

seasonPM <- c(-1,1)

# output objects
if(!skipCalc){
iriSeasonalField  <- array(NA, dim=c(nlon,nlat,length(seasons)))
ensoSeasonalField <- array(NA, dim=c(nlon,nlat,length(seasons)))
rpssSeasonalComp  <- array(NA, dim=c(nlon,nlat,length(seasons)))
rpssSeasonalComp2 <- array(NA, dim=c(nlon,nlat,length(seasons)))

iriSeasonalGroc  <- array(NA, dim=c(nlon,nlat,length(seasons)))
ensoSeasonalGroc <- array(NA, dim=c(nlon,nlat,length(seasons)))

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
	rpssSeasonalComp2[,,i] <- rpss(subIriFcast,observedTercile,subEnsoFcast)$global$field


	# calculate the GROCs
	iriSeasonalGroc[,,i]  <- grocField(subIriFcast,observedTercile,tropics=FALSE)
	ensoSeasonalGroc[,,i] <- grocField(subEnsoFcast,observedTercile,tropics=FALSE)
}


}


# make some plots that look like the paper plots currently do

titleVec <- c('DJF - FMA', 'MAM - MJJ', 'JJA - ASO', 'SON - NDJ')
pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))
textSize <- 1.5
dryColor <- adjustcolor('black', alpha=0.15)

# IRI RPSS
zMax <- max(abs(iriSeasonalField),na.rm=T)
zForce <- 0.4
iriSeasonalField[iriSeasonalField > zForce] <- zForce
iriSeasonalField[iriSeasonalField < 0] <- 0
zr <- c(0,zForce)

pdf(sprintf('%s/Seasonal/iriSeasonal.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,iriSeasonalField[,,i] * fcastMask * dryMask[,,seasons[i]],zlim=zr, col=pal, 
		main=paste0('IRI Forecast RPSS (Climatology Baseline) [', titleVec[i],']'),xlab='',ylab='',ylim=c(-60,90))
	image(lon,lat,ifelse(dryMask[,,seasons[i]]==0,1,NA),add=T,col=dryColor)
	world(add=T)
	abline(h=c(30,-30), lty=3)
	legend(-189, 92, paste0('(',letters[i],')'), cex = textSize, bty='n')
}
dev.off()


# ENSO Prob RPSS
zMax <- max(abs(ensoSeasonalField),na.rm=T)
zForce <- 0.4
ensoSeasonalField[ensoSeasonalField > zForce] <- zForce
ensoSeasonalField[ensoSeasonalField < 0] <- 0
zr <- c(0,zForce)

pdf(sprintf('%s/Seasonal/ensoSeasonal.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,ensoSeasonalField[,,i] * fcastMask * dryMask[,,seasons[i]],zlim=zr, col=pal, 
		main=paste0('Probabilistic Known-ENSO EBF RPSS [', titleVec[i],']'),xlab='',ylab='',ylim=c(-60,90))
	image(lon,lat,ifelse(dryMask[,,seasons[i]]==0,1,NA),add=T,col=dryColor)
	world(add=T)
	abline(h=c(30,-30), lty=3)
	legend(-189, 92, paste0('(',letters[i],')'), cex = textSize, bty='n')
}
dev.off()



# Alt Baseline RPSS
zMax <- max(abs(rpssSeasonalComp2),na.rm=T)
zForce <- 0.4
rpssSeasonalComp2[rpssSeasonalComp2 > zForce] <- zForce
zr <- c(-zForce,zForce)

pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

pdf(sprintf('%s/Seasonal/iriCompSeasonal.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,rpssSeasonalComp2[,,i] * fcastMask * dryMask[,,seasons[i]],zlim=zr, col=pal2, 
		main=paste0('IRI Forecast RPSS (EBF Baseline) [', titleVec[i], ']'),xlab='',ylab='',ylim=c(-60,90))
	image(lon,lat,ifelse(dryMask[,,seasons[i]]==0,1,NA),add=T,col=dryColor)
	world(add=T)
	abline(h=c(30,-30), lty=3)
	legend(-189, 92, paste0('(',letters[i],')'), cex = textSize, bty='n')
}
dev.off()

# IRI Forecast GROC
z <- iriSeasonalGroc
zMax <- 1
minCutoff <- 0.5
z[z<minCutoff] <- minCutoff
z[z>zMax] <- zMax

zr <- c(minCutoff,zMax)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

pdf(sprintf('%s/Seasonal/iriSeasonalGroc.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,z[,,i] * fcastMask * dryMask[,,seasons[i]],zlim=zr, col=pal, 
		main=paste0('IRI Forecast GROC [', titleVec[i],']'),xlab='',ylab='',ylim=c(-60,90))
	image(lon,lat,ifelse(dryMask[,,seasons[i]]==0,1,NA),add=T,col=dryColor)
	world(add=T)
	abline(h=c(30,-30), lty=3)
	legend(-189, 92, paste0('(',letters[i],')'), cex = textSize, bty='n')
}
dev.off()


# ENSO Prob GROC
z <- ensoSeasonalGroc
zMax <- 1
minCutoff <- 0.5
z[z<minCutoff] <- minCutoff
z[z>zMax] <- zMax

zr <- c(minCutoff,zMax)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

pdf(sprintf('%s/Seasonal/ensoSeasonalGroc.pdf',plotdir),16,10)
set.panel(2,2)
for(i in 1:4){
	image.plot(lon,lat,z[,,i] * fcastMask * dryMask[,,seasons[i]],zlim=zr, col=pal, 
		main=paste0('Probabilistic Known-ENSO EBF GROC [', titleVec[i],']'),xlab='',ylab='',ylim=c(-60,90))
	image(lon,lat,ifelse(dryMask[,,seasons[i]]==0,1,NA),add=T,col=dryColor)
	world(add=T)
	abline(h=c(30,-30), lty=3)
	legend(-189, 92, paste0('(',letters[i],')'), cex = textSize, bty='n')
}
dev.off()




############

# Other plots as a sanity check

# titleVec <- c('DJF - FMA', 'MAM - MJJ', 'JJA - ASO', 'SON - NDJ')

# # IRI RPSS
# zMax <- max(abs(iriSeasonalField),na.rm=T)
# zForce <- 0.4
# iriSeasonalField[iriSeasonalField > zForce] <- zForce
# zr <- c(-zForce,zForce)

# pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

# set.panel(2,2)
# for(i in 1:4){
# 	image.plot(lon,lat,iriSeasonalField[,,i],zlim=zr, col=pal2, 
# 		main=titleVec[i],xlab='',ylab='')
# 	world(add=T)
# }



# # ENSO Prob RPSS
# zMax <- max(abs(ensoSeasonalField),na.rm=T)
# zForce <- 0.4
# ensoSeasonalField[ensoSeasonalField > zForce] <- zForce
# zr <- c(-zForce,zForce)

# pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

# set.panel(2,2)
# for(i in 1:4){
# 	image.plot(lon,lat,ensoSeasonalField[,,i],zlim=zr, col=pal2, 
# 		main=titleVec[i],xlab='',ylab='')
# 	world(add=T)
# }


# # Alt Baseline RPSS
# zMax <- max(abs(rpssSeasonalComp),na.rm=T)
# zForce <- 0.2
# rpssSeasonalComp[rpssSeasonalComp > zForce] <- zForce
# zr <- c(-zForce,zForce)

# pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

# set.panel(2,2)
# for(i in 1:4){
# 	image.plot(lon,lat,rpssSeasonalComp2[,,i],zlim=zr, col=pal2, 
# 		main=titleVec[i],xlab='',ylab='')
# 	world(add=T)
# }




