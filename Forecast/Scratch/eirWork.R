# source('Forecast/Namelists/namelist_CMAP.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')


######################
# first create the ignorance field calculation (unweighted, no decomp)
######################

# inputs
fcast       <- iriForecastList
obs         <- observedTercile
landPropMat <- landPropList$prop

source('Forecast/Scratch/igFunctions.R')

iriEir <- eirCalcFull(iriForecastList,observedTercile,landPropMat)
ensoDetEir  <- eirCalcFull(ensoForecastList,observedTercile,landPropMat)
ensoProbEir <- eirCalcFull(ensoProbForecastList,observedTercile,landPropMat)

forceYR <- c(0)
pal <- brewer.pal(4,'BrBG')
lw <- 1.75
textSize <- 1.5

x  <- iriForecastList$timeMap[,3]
y  <- iriEir$globalSeries
y2 <- ensoProbEir$globalSeries[1:198]
y3 <- ensoDetEir$globalSeries[1:198]

totalVec <- c(iriEir$globalTotal, ensoProbEir$globalTotal, ensoDetEir$globalTotal)


plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Effective Interest Rate',main='Average Effective Interest Rate (Tropics)',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y2,type='l',col=pal[4],lwd=lw)
points(x,y3,type='l',col=pal[1],lwd=lw)
abline(h=0)
legend('topleft',as.character(round(totalVec,3)),text.col=c('black',pal[c(4,1)]), bty='n')



# calculate the IRI EIR the RIGN way to compare the results (field should not be the same)
yzIri  <- yzForecast(iriForecastList,observedTercile)
fIri  <- rignField(yzIri)

lon <- fIri$lon
lat <- fIri$lat

zMax <- 0.4

z1 <- iriEir$field
z2 <- fIri$eir
d1 <- z1-z2


z1[z1<0] <- NA
z2[z2<0] <- NA

z1[z1>zMax] <- zMax
z2[z2>zMax] <- zMax

zr <- c(0,zMax)

zdmax <- max(abs(d1),na.rm=T)
zr2 <- c(-zdmax, zdmax)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))
pal2 <- rev(designer.colors(256,brewer.pal(9,'RdBu')))

pdf('~/Downloads/testEIR.pdf',30,7)
set.panel(1,3)

image.plot(lon,lat,z1,zlim=zr, col=pal, main = '(1) Proper EIR Calculation')
world(add=T)

image.plot(lon,lat,z2,zlim=zr, col=pal, main = '(2) Hacky REIR Calculation')
world(add=T)

image.plot(lon,lat,d1,zlim=zr2, col=pal2, main = "(1) - (2)")
dev.off()