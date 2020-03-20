source('Forecast/Namelists/namelist_Test.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))


plotdir <- '/Users/lenssen/Dropbox/DEES/LisaWork/MastersMeeting/Paper/Figures/CodeFigures'

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')

calculate <- TRUE
# calculations for the IRI Forecast
if(calculate){
field   <- rpssDecompField(iriForecastList,observedTercile)
global  <- rpssDecompSeries(iriForecastList,observedTercile,dryMask,tropics=FALSE)
tropics <- rpssDecompSeries(iriForecastList,observedTercile,dryMask,tropics=TRUE)

totalGlobal  <- rpssTotalDecomp(iriForecastList,observedTercile,dryMask,tropics=FALSE)
totalTropics <- rpssTotalDecomp(iriForecastList,observedTercile,dryMask,tropics=TRUE)

# calculation for the enso forecast
ensoField   <- rpssDecompField(ensoForecastList,observedTercile)
ensoGlobal  <- rpssDecompSeries(ensoForecastList,observedTercile,dryMask,tropics=FALSE)
ensoTropics <- rpssDecompSeries(ensoForecastList,observedTercile,dryMask,tropics=TRUE)

ensoTotalGlobal <- rpssTotalDecomp(ensoForecastList,observedTercile,dryMask,tropics=FALSE)
ensoTotalTropics <- rpssTotalDecomp(ensoForecastList,observedTercile,dryMask,tropics=TRUE)

# calculation for the enso forecast
ensoProbField   <- rpssDecompField(ensoProbForecastList,observedTercile)
ensoProbGlobal  <- rpssDecompSeries(ensoProbForecastList,observedTercile,dryMask,tropics=FALSE)
ensoProbTropics <- rpssDecompSeries(ensoProbForecastList,observedTercile,dryMask,tropics=TRUE)

ensoProbTotalGlobal  <- rpssTotalDecomp(ensoProbForecastList,observedTercile,dryMask,tropics=FALSE)
ensoProbTotalTropics <- rpssTotalDecomp(ensoProbForecastList,observedTercile,dryMask,tropics=TRUE)

# calculation for the enso forecast
ensoRealField   <- rpssDecompField(ensoRealForecastList,observedTercile)
ensoRealGlobal  <- rpssDecompSeries(ensoRealForecastList,observedTercile,dryMask,tropics=FALSE)
ensoRealTropics <- rpssDecompSeries(ensoRealForecastList,observedTercile,dryMask,tropics=TRUE)

ensoRealTotalGlobal <- rpssTotalDecomp(ensoRealForecastList,observedTercile,dryMask,tropics=FALSE)
ensoRealTotalTropics <- rpssTotalDecomp(ensoRealForecastList,observedTercile,dryMask,tropics=TRUE)

# calculate for the realistic det forecast list
ensoRealDetField    <- rpssDecompField(ensoRealDetForecastList,observedTercile)
ensoRealDetGlobal   <- rpssDecompSeries(ensoRealDetForecastList,observedTercile,dryMask,tropics=FALSE)
ensoRealDetTropics  <- rpssDecompSeries(ensoRealDetForecastList,observedTercile,dryMask,tropics=TRUE)

ensoRealDetTotalGlobal <- rpssTotalDecomp(ensoRealDetForecastList,observedTercile,dryMask,tropics=FALSE)
ensoRealDetTotalTropics <- rpssTotalDecomp(ensoRealDetForecastList,observedTercile,dryMask,tropics=TRUE)

}
######################################
# compare the tropics rpss series
######################################
pal <- brewer.pal(4,'BrBG')
lw <- 1.75
textSize <- 1.5

x  <- tropics$timeMap[,3]
y  <- tropics$skill$rpss
y2 <- ensoTropics$skill$rpss[1:198]
y3 <- ensoProbTropics$skill$rpss[1:198]
y4 <- ensoRealTropics$skill$rpss[1:198]
y5 <- ensoRealDetTropics$skill$rpss[1:198]

# with det included
pdf(sprintf('%s/tropicsRpss.pdf',plotdir),12,6.5)
plot(x,y,ylim=range(y,y2,y3,y4,y5,na.rm=T),type='l',lwd=lw, lty=3,
	xlab='',ylab='Ranked Probability Skill Score',main='Average Tropics Skill',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
points(x,y4,type='l',col=pal[4],lwd=lw)
points(x,y3,type='l',col=pal[3],lwd=lw)
points(x,y2,type='l',col=pal[2],lwd=lw)
points(x,y5,type='l',col=pal[1],lwd=lw)
abline(h=0)
# legend('topleft', 
# 	c('IRI Forecast (Lead 1)',
# 	  'ENSO Probabilistic Forecast (Oracle)',
# 	  'ENSO Probabilistic Forecast (Realistic)', 
# 	  'ENSO Deterministic Forecast (Oracle)',
# 	  'ENSO Deterministic Forecast (Realistic)'),
# 	col=c('black',pal[c(3,4,2,1)]),lwd=lw,lty=c(3,1,1,1,1))
dev.off()


######################################
# compare the tropics resolution series
######################################
pal <- brewer.pal(4,'BrBG')
lw <- 1.75

x  <- tropics$timeMap[,3]
y  <- tropics$decomp73$res
y2 <- ensoTropics$decomp73$res[1:198]
y3 <- ensoProbTropics$decomp73$res[1:198]
y4 <- ensoRealTropics$decomp73$res[1:198]
y5 <- ensoRealDetTropics$decomp73$res[1:198]


pdf(sprintf('%s/tropicsRes.pdf',plotdir),12,6.5)
plot(x,y,ylim=range(y,y2,y3,y4,y5,na.rm=T),type='l',lwd=lw, lty=3,
	xlab='',ylab='Total Resolution',main='Average Tropics Resolution',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
points(x,y4,type='l',col=pal[4],lwd=lw)
points(x,y3,type='l',col=pal[3],lwd=lw)
points(x,y2,type='l',col=pal[2],lwd=lw)
points(x,y5,type='l',col=pal[1],lwd=lw)
abline(h=0)
# legend('topleft', 
# 	c('IRI Forecast (Lead 1)',
# 	  'ENSO Probabilistic Forecast (Oracle)',
# 	  'ENSO Probabilistic Forecast (Realistic)', 
# 	  'ENSO Deterministic Forecast (Oracle)',
# 	  'ENSO Deterministic Forecast (Realistic)'),
# 	col=c('black',pal[c(3,4,2,1)]),lwd=lw,lty=c(3,1,1,1,1))

dev.off()

######################################
# compare the tropics reliability series
######################################
pal <- brewer.pal(4,'BrBG')
lw <- 1.75

x  <- tropics$timeMap[,3]
y  <- -tropics$decomp73$rel
y2 <- -ensoTropics$decomp73$rel[1:198]
y3 <- -ensoProbTropics$decomp73$rel[1:198]
y4 <- -ensoRealTropics$decomp73$rel[1:198]
y5 <- -ensoRealDetTropics$decomp73$rel[1:198]


pdf(sprintf('%s/tropicsRel.pdf',plotdir),12,6.5)
plot(x,y,ylim=range(y,y2,y3,y4,y5,na.rm=T),type='l',lwd=lw, lty=3,
	xlab='',ylab='Negative Total Reliability',main='Average Tropics Reliability',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
points(x,y4,type='l',col=pal[4],lwd=lw)
points(x,y3,type='l',col=pal[3],lwd=lw)
points(x,y2,type='l',col=pal[2],lwd=lw)
points(x,y5,type='l',col=pal[1],lwd=lw)
abline(h=0)
# legend('topleft', 
# 	c('IRI Forecast (Lead 1)',
# 	  'ENSO Probabilistic Forecast (Oracle)',
# 	  'ENSO Probabilistic Forecast (Realistic)', 
# 	  'ENSO Deterministic Forecast (Oracle)',
# 	  'ENSO Deterministic Forecast (Realistic)'),
# 	col=c('black',pal[c(3,4,2,1)]),lwd=lw,lty=c(3,1,1,1,1))

dev.off()





######################################
# RPSS Field
######################################
lon <- field$lon
lat <- field$lat

maxCutoff <- 0.25

z   <- ensoProbField$skill$rpss

z[z>maxCutoff] <- maxCutoff
z[z<0] <- 0

zr <- c(0,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

pdf(sprintf('%s/fieldRpssProb.pdf',plotdir),10,6)
image.plot(lon,lat,z,col=pal,ylim=c(-60,90),main="ENSO Oracle Probabilistic Forecast RPSS",xlab='',ylab='')
world(add=T)
dev.off()

######################################
# Resolution Field
######################################
lon <- field$lon
lat <- field$lat

maxCutoff <- 0.15

z   <- ensoProbField$decomp73$res

z[z>maxCutoff] <- maxCutoff

zr <- c(0,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

pdf(sprintf('%s/fieldResProb.pdf',plotdir),10,6)
image.plot(lon,lat,z,col=pal,ylim=c(-60,90),zlim=zr,main="ENSO Oracle Probabilistic Forecast Annual Total Resolution",xlab='',ylab='')
world(add=T)
dev.off()


######################################
# Plot some example Forecasts
######################################

# 2015.000          5.000
# FIND OUT WHAT THE ENSO FORECAST WAS THEN!

iriInd <- 180
iriForecastList$timeMap[iriInd,]
oceanMask <- ifelse(landPropList$prop==0,1,NA)
pdf(sprintf('%s/fcastDetOracle.pdf',plotdir),10,6)
iriPlot(lon,lat,ensoForecastList$forecast[,,iriInd,],oceanMask,ylim=c(-60,90))
dev.off()

pdf(sprintf('%s/fcastProbOracle.pdf',plotdir),10,6)
iriPlot(lon,lat,ensoProbForecastList$forecast[,,iriInd,],oceanMask,ylim=c(-60,90))
dev.off()

pdf(sprintf('%s/fcastDetReal.pdf',plotdir),10,6)
iriPlot(lon,lat,ensoRealDetForecastList$forecast[,,iriInd,],oceanMask,ylim=c(-60,90))
dev.off()

pdf(sprintf('%s/fcastProbReal.pdf',plotdir),10,6)
iriPlot(lon,lat,ensoRealForecastList$forecast[,,iriInd,],oceanMask,ylim=c(-60,90))
dev.off()


###########################################################################
# Make a table of the total performance of the 5 different forecasts
###########################################################################

# in the order of the legend
totalTab <- rbind(
totalTropics,
ensoProbTotalTropics,
ensoRealTotalTropics,
ensoTotalTropics,
ensoRealDetTotalTropics)

