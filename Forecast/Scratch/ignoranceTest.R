source('Forecast/Namelists/namelist_CMAP.Rnl')

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


# calculation for the det enso forecast
ensoField   <- rpssDecompField(ensoForecastList,observedTercile)
ensoGlobal  <- rpssDecompSeries(ensoForecastList,observedTercile,dryMask,tropics=FALSE)
ensoTropics <- rpssDecompSeries(ensoForecastList,observedTercile,dryMask,tropics=TRUE)

ensoTotalGlobal <- rpssTotalDecomp(ensoForecastList,observedTercile,dryMask,tropics=FALSE)
ensoTotalTropics <- rpssTotalDecomp(ensoForecastList,observedTercile,dryMask,tropics=TRUE)

# calculation for the prob enso forecast
ensoProbField   <- rpssDecompField(ensoProbForecastList,observedTercile)
ensoProbGlobal  <- rpssDecompSeries(ensoProbForecastList,observedTercile,dryMask,tropics=FALSE)
ensoProbTropics <- rpssDecompSeries(ensoProbForecastList,observedTercile,dryMask,tropics=TRUE)

ensoProbTotalGlobal  <- rpssTotalDecomp(ensoProbForecastList,observedTercile,dryMask,tropics=FALSE)
ensoProbTotalTropics <- rpssTotalDecomp(ensoProbForecastList,observedTercile,dryMask,tropics=TRUE)


# calculate the EIR (ignorance) calculations

iriEir      <- eir(iriForecastList,observedTercile)
ensoEir     <- eir(ensoForecastList,observedTercile)
ensoProbEir <- eir(ensoProbForecastList,observedTercile)
ensoRealEir <- eir(ensoRealForecastList,observedTercile)

# series plot
x  <- iriForecastList$timeMap[,3]

y  <- iriEir$tropics$series
y2 <- ensoEir$tropics$series[1:198]
y3 <- ensoProbEir$tropics$series[1:198]

y4 <- ensoRealEir$tropics$series

pal <- brewer.pal(4,'BrBG')

plot(x,y,type='l',lwd=2)
grid(lwd=1.5)
points(x, y2,type='l',lwd=2,col=pal[2])
points(x, y3,type='l',lwd=2,col=pal[3])
points(x, y4, type='l',lwd=2,col=pal[4])
abline(h=0)


# mock up of the field figure I would use in the paper
lon <- iriForecastList$lon
lat <- iriForecastList$lat

par(mfrow=c(2,1))

z <- iriEir$global$averageField
zMax <- 0.20

noSkill <- which(z<0)
greyMat <- matrix(NA, nrow=length(lon),ncol=length(lat))
greyMat[noSkill] <- 1
z[noSkill] <- NA

z[z>zMax] <- zMax

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))
image.plot(lon,lat,z,ylim=c(-60,90),col=pal)
image(lon,lat,greyMat,col='grey80',add=T)
world(add=T)

z <- ensoProbEir$global$averageField
zMax <- 0.20

noSkill <- which(z<0)
greyMat <- matrix(NA, nrow=length(lon),ncol=length(lat))
greyMat[noSkill] <- 1
z[noSkill] <- NA

z[z>zMax] <- zMax

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))
image.plot(lon,lat,z,ylim=c(-60,90),col=pal)
image(lon,lat,greyMat,col='grey80',add=T)
world(add=T)

