source('Forecast/Namelists/namelist_Test.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')

plotdir <- '/Users/lenssen/Dropbox/DEES/LisaWork/Meetings/CompIssues'
analysis <- FALSE


# get the rpss with the ENSO-based alt forecast
rpssComp <- rpss(iriForecastList,observedTercile,ensoProbForecastList)


# check to make sure this function is doing what I think it is...
checkSeries <- rpssDecompSeries(iriForecastList, observedTercile, dryMask)


# get the field of iri Forecast counts
fcastCounts <- apply(iriForecastList$forecast[,,,1],c(1,2), function(x) sum(!is.na(x)))
fcastMask <- ifelse(fcastCounts < 50, NA, 1)

# just plot the field of rpssComp
lon <- iriForecastList$lon
lat <- iriForecastList$lat

z <- rpssComp$global$field * fcastMask

zMax <- max(abs(z), na.rm=T)

zCap <- 0.2
z[z>zCap] <- zCap

zr <- c(-zCap, zCap)

pal <- designer.colors(256, brewer.pal(11,'PiYG'))

pdf('~/Desktop/scratchPlots/ensoRefField.pdf',10,7)
image.plot(lon, lat, z, zlim=zr, col=pal, 
	xlab='', ylab='', main='IRI RPSS (Probabilstic ENSO reference)')
world(add=T)
abline(h=c(30,-30), lty=3)
dev.off()

## Ranked Probability Skill Score

# draw the field plot


difference <- rpssComp$global$field - rpssClim$global$field

zr <- c(-0.3,0.1)
zSeq <- seq(-0.3,0.3,length=501)
pal <- redBlue(length(zSeq))
zInd <- which(zSeq >= zr[1] & zSeq <= zr[2])

pdf(sprintf('%s/rpssAverageFieldComp.pdf',plotdir),10,7)
image.plot(lon,lat,difference,col=pal[zInd],ylim=c(-60,90),zlim=zr,
	main='Annual Average RPSS Increase Due to Impacts Validation',
	xlab='',ylab='')
world(add=T)
dev.off()

# best season plot

lon <- iriForecastList$lon
lat <- iriForecastList$lat

difference <- rpssComp$bestMat - rpssClim$bestMat

pdf(sprintf('%s/rpssBestFieldComp.pdf',plotdir),10,7)
image.plot(lon,lat,difference,col=pal[zInd],ylim=c(-60,90),zlim=c(-0.25,0.25),
	main='Best Season RPSS Increase Due to Impacts Validation',
	xlab='',ylab='')
world(add=T)
dev.off()

# Determine the increase in skill in the IRI forecast when both the IRI and the 
# ENSO forecast are verified relative to climatology

# annual average
lon <- iriForecastList$lon
lat <- iriForecastList$lat

zr <- c(-0.1,0.1)

zmax <- max(abs(zr))
zSeq <- seq(-zmax,zmax,length=501)
pal <- redBlue(length(zSeq))
zInd <- which(zSeq >= zr[1] & zSeq <= zr[2])

pdf(sprintf('%s/rpssAverageFieldComp2.pdf',plotdir),10,7)
difference <- rpssClim$global$field - rpssComp2$global$field
image.plot(lon,lat,difference,col=pal[zInd],zlim=zr,ylim=c(-60,90),
	main='IRI-ENSO Impacts Annual Average RPSS Skill (Climatology as Reference)',
	xlab='',ylab='')
world(add=T)
dev.off()


# best season difference
lon <- iriForecastList$lon
lat <- iriForecastList$lat

zr <- c(-0.3,0.3)

zmax <- max(abs(zr))
zSeq <- seq(-zmax,zmax,length=501)
pal <- redBlue(length(zSeq))
zInd <- which(zSeq >= zr[1] & zSeq <= zr[2])

pdf(sprintf('%s/rpssBestFieldComp2.pdf',plotdir),10,7)
difference <- rpssClim$bestMat- rpssComp2$bestMat
image.plot(lon,lat,difference,col=pal[zInd],zlim=zr,ylim=c(-60,90),
	main='IRI-ENSO Impacts Best Season RPSS Skill (Climatology as Reference)',
	xlab='',ylab='')
world(add=T)
dev.off()

# plot the global result


pdf(sprintf('%s/rpssGlobalComp.pdf',plotdir),10,14)
set.panel(2,1)
yr <- c(-0.1,0.15)
time <- rpssComp$timeMap[,3]
plot(time,rpssComp$global$series,type='l',ylim=yr,col='darkgreen',lwd=1.5,
	xlab='Year', ylab='RPSS',
	main='Global Forecast Skill with Climatology and Impacts Reference Forecasts')
points(time,rpssClim$global$series,type='l',lwd=1.5)
points(time,rpssComp2$global$series[1:198],type='l',lwd=1.5,col='blue',lty=3)
abline(h=0)

legend('bottomright', c('IRI Forecast (ENSO Reference)', 
						'IRI Forecast (Climatology Reference)',
						'ENSO Forecast (Climatology Reference)'),
	lwd=2,col=c('darkgreen','black','blue'),lty=c(1,1,3))


difference <- rpssComp$global$series - rpssClim$global$series
difference2 <- rpssClim$global$series - rpssComp2$global$series[1:198]
plot(time,difference,type='l',col='firebrick4',ylim=yr,
	xlab='Year', ylab='RPSS Gain due to IRI Forecast')
abline(h=0)
points(time,difference,type='l',col='firebrick4',lwd=1.5)
points(time,difference2,type='l',col='firebrick4',lwd=1.5,lty=3)

legend('bottomright', c('IRI (ENSO Reference) - IRI (Climatology Reference)', 
						'IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)'),
	lwd=2,col=c('firebrick4','firebrick4'),lty=c(1,3))

dev.off()

# plot the tropical result

pdf(sprintf('%s/rpssTropicsComp.pdf',plotdir),10,14)
set.panel(2,1)
yr <- c(-0.15,0.2)
time <- rpssComp$timeMap[,3]
plot(time,rpssComp$tropics$series,type='l',ylim=yr,col='darkgreen',lwd=1.5,
	xlab='Year', ylab='RPSS',
	main='Tropics Forecast Skill with Climatology and Impacts Reference Forecasts')
points(time,rpssClim$tropics$series,type='l',lwd=1.5)
points(time,rpssComp2$tropics$series[1:198],type='l',lwd=1.5,col='blue',lty=3)
abline(h=0)

legend('bottomright', c('IRI Forecast (ENSO Reference)', 
						'IRI Forecast (Climatology Reference)',
						'ENSO Forecast (Climatology Reference)'),
	lwd=2,col=c('darkgreen','black','blue'),lty=c(1,1,3))


difference <- rpssComp$tropics$series - rpssClim$tropics$series
difference2 <- rpssClim$tropics$series - rpssComp2$tropics$series[1:198]
plot(time,difference,type='l',col='firebrick4',ylim=yr,
	xlab='Year', ylab='RPSS Gain due to IRI Forecast')
abline(h=0)
points(time,difference,type='l',col='firebrick4',lwd=1.5)
points(time,difference2,type='l',col='firebrick4',lwd=1.5,lty=3)

legend('bottomright', c('IRI (ENSO Reference) - IRI (Climatology Reference)', 
						'IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)'),
	lwd=2,col=c('firebrick4','firebrick4'),lty=c(1,3))

dev.off()

## Effective Interest Rate

# draw the field difference Plotplot
pdf(sprintf('%s/eirAverageFieldComp.pdf',plotdir),10,7)
lon <- iriForecastList$lon
lat <- iriForecastList$lat

difference <- eirComp$global$averageField - eirClim$global$averageField
image.plot(lon,lat,difference,col=redBlue(),zlim=c(-0.25,0.25),ylim=c(-60,90),
	main='[IRI Skill (ENSO Reference) - IRI Skill (Climatology Reference)] Annual Average EIR',
	xlab='',ylab='')
world(add=T)
dev.off()

pdf(sprintf('%s/eirBestFieldComp.pdf',plotdir),10,7)
lon <- iriForecastList$lon
lat <- iriForecastList$lat

zr <- c(-0.5,1)
zSeq <- seq(-1,1,length=501)
pal <- redBlue(length(zSeq))
zInd <- which(zSeq >= zr[1] & zSeq <= zr[2])

difference <- eirComp$global$bestField - eirClim$global$bestField
image.plot(lon,lat,difference,col=pal[zInd],zlim=zr,ylim=c(-60,90),
	main='[IRI (ENSO Reference) - IRI (Climatology Reference)] Best Seasosn EIR',
	xlab='',ylab='')
world(add=T)
dev.off()

# Determine the increase in skill in the IRI forecast when both the IRI and the 
# ENSO forecast are verified relative to climatology

# annual average
lon <- iriForecastList$lon
lat <- iriForecastList$lat

zr <- c(-0.25,0.25)

zmax <- max(abs(zr))
zSeq <- seq(-zmax,zmax,length=501)
pal <- redBlue(length(zSeq))
zInd <- which(zSeq >= zr[1] & zSeq <= zr[2])

pdf(sprintf('%s/eirAverageFieldComp2.pdf',plotdir),10,7)
difference <- eirClim$global$averageField - eirComp2$global$averageField
image.plot(lon,lat,difference,col=pal[zInd],zlim=zr,ylim=c(-60,90),
	main='[IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)] Annual Average EIR',
	xlab='',ylab='')
world(add=T)
dev.off()

# best season difference
lon <- iriForecastList$lon
lat <- iriForecastList$lat

zr <- c(-0.5,1)

zmax <- max(abs(zr))
zSeq <- seq(-zmax,zmax,length=501)
pal <- redBlue(length(zSeq))
zInd <- which(zSeq >= zr[1] & zSeq <= zr[2])

pdf(sprintf('%s/eirBestFieldComp2.pdf',plotdir),10,7)
difference <- eirClim$global$bestField - eirComp2$global$bestField
image.plot(lon,lat,difference,col=pal[zInd],zlim=zr,ylim=c(-60,90),
	main='[IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)] Best Season EIR',
	xlab='',ylab='')
world(add=T)
dev.off()



# plot the global result
pdf(sprintf('%s/eirGlobalComp.pdf',plotdir),10,14)
set.panel(2,1)
yr <- c(-0.15,0.15)
time <- rpssComp$timeMap[,3]
plot(time,eirComp$global$series,type='l',ylim=yr,col='darkgreen',lwd=1.5,
	xlab='Year', ylab='Effective Interest Rate',
	main='Global Forecast Skill with Climatology and Impacts Reference Forecasts')
points(time,eirClim$global$series,type='l',lwd=1.5)
points(time,eirComp2$global$series[1:198],type='l',lwd=1.5,lty=3,col='blue')
abline(h=0)

legend('bottomright', c('IRI Forecast (ENSO Reference)', 
						'IRI Forecast (Climatology Reference)',
						'ENSO Forecast (Climatology Reference)'),
	lwd=2,col=c('darkgreen','black','blue'),lty=c(1,1,3))


difference <- eirComp$global$series - eirClim$global$series
difference2 <- eirClim$global$series - eirComp2$global$series[1:198]
plot(time,difference,type='l',col='firebrick4',ylim=yr,
	xlab='Year', ylab='EIR Gain by IRI Forecast')
abline(h=0)
points(time,difference,type='l',col='firebrick4',lwd=1.5)
points(time,difference2,type='l',col='firebrick4',lwd=1.5,lty=3)

legend('bottomright', c('IRI (ENSO Reference) - IRI (Climatology Reference)', 
						'IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)'),
	lwd=2,col=c('firebrick4','firebrick4'),lty=c(1,3))

dev.off()

# plot the tropical result

pdf(sprintf('%s/eirTropicsComp.pdf',plotdir),10,14)
set.panel(2,1)
yr <- c(-0.15,0.15)
time <- rpssComp$timeMap[,3]
plot(time,eirComp$tropics$series,type='l',ylim=yr,col='darkgreen',lwd=1.5,
	xlab='Year', ylab='Effective Interest Rate',
	main='Tropics Forecast Skill with Climatology and Impacts Reference Forecasts')
points(time,eirClim$tropics$series,type='l',lwd=1.5)
points(time,eirComp2$tropics$series[1:198],type='l',lwd=1.5,lty=3,col='blue')
abline(h=0)

legend('bottomright', c('IRI Forecast (ENSO Reference)', 
						'IRI Forecast (Climatology Reference)',
						'ENSO Forecast (Climatology Reference)'),
	lwd=2,col=c('darkgreen','black','blue'),lty=c(1,1,3))


difference <- eirComp$tropics$series - eirClim$global$series
difference2 <- eirClim$tropics$series - eirComp2$global$series[1:198]
plot(time,difference,type='l',col='firebrick4',ylim=yr,
	xlab='Year', ylab='EIR Gain by IRI Forecast')
abline(h=0)
points(time,difference,type='l',col='firebrick4',lwd=1.5)
points(time,difference2,type='l',col='firebrick4',lwd=1.5,lty=3)

legend('bottomright', c('IRI (ENSO Reference) - IRI (Climatology Reference)', 
						'IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)'),
	lwd=2,col=c('firebrick4','firebrick4'),lty=c(1,3))

dev.off()








