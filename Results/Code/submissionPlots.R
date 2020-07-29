# Load needed final data products (skill scores must be loaded in run script)
load('Data/RawProcessed/ninaSeasonal.Rda')
load('Data/RawProcessed/landProp_0.5.Rda')
load('Data/RawProcessed/cruDataSummary.Rda')
load('Data/Impacts/Cru_0.5_1951_2016_Data/ninaSeasonalInd.Rda')

# load in the forecast info
load(sprintf('Data/Forecast/%s_Data/iriForecast.Rda',fcastProject))
load(sprintf('Data/Forecast/%s_Data/ensoForecast.Rda',fcastProject))
load(sprintf('Data/Forecast/%s_Data/obsTercile.Rda',fcastProject))
load(sprintf('Data/Forecast/%s_Data/dryMask.Rda',fcastProject))




###############################################################################
# get all of the forecast verification things loaded/run
###############################################################################

# make a mask for field plots (mostly to exclude sahara and strange coast points)
fcastCounts <- apply(iriForecastList$forecast[,,,1],c(1,2), function(x) sum(!is.na(x)))
fcastMask <- ifelse(fcastCounts < 50, NA, 1)



###############################################################################
# (1) Cartoon
###############################################################################

###############################################################################
# (2) Cru Data summary plots
###############################################################################

# series
tSeason <- rep(1:12, length(propCovered)/12)
plotInds <- which(tSeason==3)

fullCol <- adjustcolor('blue',alpha=0.7)
partCol <- adjustcolor('red',alpha=0.4)
noCol   <- adjustcolor('black',alpha=0.2)

x  <- tYearCru
y1 <- propCovered[plotInds]
y2 <- propCovered50[plotInds]
y3 <- propCovered30[plotInds]

yr <- c(0,1)

pdf(sprintf('%s/cruGlobalCoverageSeries.pdf',plotdir),10,4.5)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)
plot(x,y1,type='l',col='black',lwd=2,ylim=yr,
	xlab='Year',ylab='Coverage Proportion',
	cex.lab=textSize, cex.axis=1.25, cex.main=textSize, cex.sub=textSize)
points(x,y2,type='l',col='blue',lwd=2)
points(x,y3,type='l',col='red',lwd=2)
grid()
legend(1950,0.4,c('Global','50N-50S','30N-30S'),
		col=c('black','blue','red'),lwd=2,bty='n')
legend('topleft','(b)', bty = "n", )
dev.off()

# field

pdf(sprintf('%s/cruGlobalCoverageField.pdf',plotdir),10,6)
image(cruLon,cruLat,covered15, col=fullCol,ylim=c(-60,90),
				xlab='',ylab='',
				main='CRU TS 4.01 Precipitation Coverage',
	cex.lab=textSize, cex.axis=1.25, cex.main=textSize, cex.sub=textSize)
image(cruLon,cruLat,notCovered90,col=noCol,add=T)
image(cruLon,cruLat,notCovered15,col=partCol,add=T)
world(add=T)
legend('bottomleft', c('No Coverage', 'Coverage in 1990', 'Coverage in 2015'),
	col=c(noCol,partCol,fullCol),lwd=8,bty = "n")
abline(h=c(30,-30), lty=3)
legend('topleft','(a)', bty = "n", )
dev.off()


###############################################################################
# (3) Example Maps
###############################################################################

# Auto-generated from running Impacts analysis
system(sprintf('cp Figures/Impacts/%s_Figures/empiricalProbs/01NinaAboveDJFProb.pdf %s/03a_NinaAboveDJFProb.pdf',
	impactsProject, plotdir))

system(sprintf('cp Figures/Impacts/%s_Figures/empiricalProbs/01NinaBelowDJFProb.pdf %s/03b_NinaBelowDJFProb.pdf',
	impactsProject, plotdir))

###############################################################################
# (5) Tropics Resolution/Discrimination
###############################################################################
pal <- brewer.pal(4,'BrBG')

resNino <- -0.0005
disNino <- 0.455

pdf(sprintf('%s/04_resolutionDiscrimination.pdf',plotdir),10,14)
par(mfrow=c(2,1),mar=c(5, 5, 4, 3) + 0.1)

forceYR <- c(0)

x  <- global$timeMap[,3]
y  <- global$decomp73$res
y2 <- ensoGlobal$decomp73$res[1:198]
y3 <- ensoProbGlobal$decomp73$res[1:198]

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Resolution Score',main='Mean Tropics Resolution',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y3,type='l',col=pal[4],lwd=lw)
points(x,y2,type='l',col=pal[1],lwd=lw)
abline(h=0)
legend('topright',as.character(round(as.numeric(totalTabTropics[c(1,2,4),4]),3)),text.col=c('black',pal[c(4,1)]), bty='n')
legend('topleft', '(a)', cex = textSize, bty='n')

addNino(y=c(-100,resNino,resNino,-100))

# discrimination plot
forceYR <- c(0.46)

x  <- global$timeMap[,3]

y  <- iriGrocTropics$series
y2 <- ensoDetGrocTropics$series[1:198]
y3 <- ensoProbGrocTropics$series[1:198]


plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Generalized ROC Score',main='Mean Tropics Discrimination',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y3,type='l',col=pal[4],lwd=lw)
points(x,y2,type='l',col=pal[1],lwd=lw)
points(x,y,lwd=lw,type='l')
abline(h=0.5)

legend(2000.5,0.75, 
	c('IRI Forecast',
	  'Probabilistic Known-ENSO EBF',
	  'Deterministic Known-ENSO EBF'),
	col=c('black',pal[c(4,1)]),lwd=lw,lty=c(1,1,1,1),bty='n')

legend('topright',as.character(round(as.numeric(totalGrocTropics[c(1,2,4)]),4)),text.col=c('black',pal[c(4,1)]), bty='n')
legend('topleft', '(b)', cex = textSize, bty='n')

addNino(y=c(-100,disNino,disNino,-100))

dev.off()

###############################################################################
# (6) Resolution Field of the two forecasts of interest
###############################################################################

lon <- field$lon
lat <- field$lat

maxCutoff <- 0.15

z <- array(NA, dim=c(length(lon),length(lat),2))

z[,,1]  <- field$decomp73$res * fcastMask
z[,,2]  <- ensoProbField$decomp73$res * fcastMask

z[z>maxCutoff] <- maxCutoff
z[z<0.01] <- 0

zr <- c(0,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

mainVec <- c("IRI Forecast Resolution",
			 "ENSO Realistic Forecast Resolution",
			 "Probabilistic Known-ENSO EBF Resolution",
			 "Deterministic Known-ENSO EBF Resolution")

pdf(sprintf('%s/05_fieldRes.pdf',plotdir),10,12)
par(mfrow=c(2,1))

image.plot(lon,lat,z[,,1],col=pal,ylim=c(-60,90),main=mainVec[1],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)	
abline(h=c(30,-30), lty=3)
legend(-197, 95, '(a)', cex = textSize, bty='n')


image.plot(lon,lat,z[,,2],col=pal,ylim=c(-60,90),main=mainVec[3],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)
abline(h=c(30,-30), lty=3)
legend(-197, 95, '(b)', cex = textSize, bty='n')

dev.off()

###############################################################################
# (6) Geographic Discrimination of the the IRI and prob ENSO
###############################################################################

lon <- field$lon
lat <- field$lat

minCutoff <- 0.5
zMax <- 0.85

z <- array(NA, dim=c(length(lon),length(lat),2))

z[,,1]  <- iriGrocField * fcastMask
z[,,2]  <- ensoProbGrocField * fcastMask


z[z<minCutoff] <- minCutoff
z[z>zMax] <- zMax

zr <- c(minCutoff,zMax)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

mainVec <- c("IRI Forecast Discrimination (GROC)",
			 "Probabilistic Known-ENSO EBF Discrimination (GROC)")

pdf(sprintf('%s/06_fieldGroc.pdf',plotdir),10,12)
par(mfrow=c(2,1))

image.plot(lon,lat,z[,,1],col=pal,ylim=c(-60,90),main=mainVec[1],xlab='',ylab='',zlim=zr,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)	
abline(h=c(30,-30), lty=3)
legend(-197, 95, '(a)', cex = textSize, bty='n')

image.plot(lon,lat,z[,,2],col=pal,ylim=c(-60,90),main=mainVec[2],xlab='',ylab='',zlim=zr,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
world(add=T)	
abline(h=c(30,-30), lty=3)
legend(-197, 95, '(b)', cex = textSize, bty='n')

dev.off()

###############################################################################
# (7) Reliability Time Series
###############################################################################

pal <- brewer.pal(4,'BrBG')
relNino <- 0.011

pdf(sprintf('%s/07_relSeries.pdf',plotdir),10,7)
par(mfrow=c(1,1),mar=c(5, 5, 4, 3) + 0.1)

forceYR <- c(0.01)

x  <- tropics$timeMap[,3]
y  <- -tropics$decomp73$rel
y2 <- -ensoTropics$decomp73$rel[1:198]
y3 <- -ensoProbTropics$decomp73$rel[1:198]
plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Negative Reliability Score',main='Mean Tropics Reliability',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y3,type='l',col=pal[4],lwd=lw)
points(x,y2,type='l',col=pal[2],lwd=lw)
abline(h=0)

legend('bottomleft', 
	c('IRI Forecast',
	  'Probabilistic Known-ENSO EBF',
	  'Deterministic Known-ENSO EBF'),
	col=c('black',pal[c(4,2)]),lwd=lw,lty=c(1,1,1,1),bty='n')

legend('bottomright',as.character(round(as.numeric(totalTabTropics[c(1,2,4),3]),4)),text.col=c('black',pal[c(4,2)]), bty='n')

addNino(y=c(100,relNino,relNino,100))

dev.off()


###############################################################################
# (8) Reliability Diagrams
###############################################################################

# tropics without realistic included
pdf(sprintf('%s/08_reliabilityDiagramsTropicsTrim.pdf',plotdir),18,5.5)

par(mfrow=c(1,3))

yr <- c(0,245000)

reliabilityDiagram(iriReliability$tropics,yr,main='IRI Forecast Tropics Reliability')
legend(-0.09, 1.06, '(a)', bty='n', cex=textSize)
#reliabilityDiagram(ensoRealReliability$tropics,yr,main='ENSO Realistic Tropics Reliability')
reliabilityDiagram(ensoProbReliability$tropics,yr,main='Probabilistic Known-ENSO EBF Tropics Reliability')
legend(-0.09, 1.06, '(b)', bty='n', cex=textSize)

reliabilityDiagram(ensoReliability$tropics,yr,main='Deterministic Known-ENSO EBF Tropics Reliability')
legend(-0.09, 1.06, '(c)', bty='n', cex=textSize)

dev.off()

###############################################################################
# (9) RPSS Time Series for global and tropics
###############################################################################
pal <- brewer.pal(4,'BrBG')
forceYR <- c(0.2,-0.75)
rpsNino <- -0.75

pdf(sprintf('%s/09_rpssSeriesTrim.pdf',plotdir),10,14)
par(mfrow=c(2,1), mar=c(5, 5, 4, 3) + 0.1)

x  <- global$timeMap[,3]
y  <- global$skill$rpss
y2 <- ensoGlobal$skill$rpss[1:198]
y3 <- ensoProbGlobal$skill$rpss[1:198]

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Ranked Probability Skill Score',main='Mean Global Skill (RPSS, Climatology Reference)',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y3,type='l',col=pal[4],lwd=lw)
points(x,y2,type='l',col=pal[1],lwd=lw)
abline(h=0)
legend('topleft' , '(a)', cex = textSize, bty='n')
legend(2015,-0.44,as.character(round(as.numeric(totalTabGlobal[c(1,2,4),1]),4)),text.col=c('black',pal[c(4,1)]), bty='n')

addNino(y=c(-100,rpsNino,rpsNino,-100))


# add the alt to the plot

x  <- tropics$timeMap[,3]
y  <- tropics$skill$rpss
y2 <- ensoTropics$skill$rpss[1:198]
y3 <- ensoProbTropics$skill$rpss[1:198]

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Ranked Probability Skill Score',main='Mean Tropics Skill (RPSS, Climatology Reference)',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)

points(x,y3,type='l',col=pal[4],lwd=lw)
points(x,y2,type='l',col=pal[1],lwd=lw)
abline(h=0)

legend(2001,0.25,
	c('IRI Forecast',
	  'Probabilistic Known-ENSO EBF',
	  'Deterministic Known-ENSO EBF'),
	col=c('black',pal[c(4,1)]),lwd=lw,lty=c(1,1,1,1),bty='n')
legend('topleft', '(b)', cex = textSize, bty='n')
legend(2015,-0.44,as.character(round(as.numeric(totalTabTropics[c(1,2,4),1]),4)),text.col=c('black',pal[c(4,1)]), bty='n')

addNino(y=c(-100,rpsNino,rpsNino,-100))


dev.off()

###############################################################################
# (10) RPSS Fields
###############################################################################


lon <- field$lon
lat <- field$lat

maxCutoff <- 0.25 #0.25

z <- array(NA, dim=c(length(lon),length(lat),2))

z[,,1]  <- field$skill$rpss * fcastMask
z[,,2]  <- ensoProbField$skill$rpss * fcastMask

z[z>maxCutoff] <- maxCutoff
z[z<0] <- 0

zr <- c(0,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))


z2 <- rpssComp$global$field * fcastMask

zMax2 <- max(abs(z), na.rm=T)

zCap2 <- 0.2 #0.2
z[z>zCap2] <- zCap2

zr2 <- c(-zCap2, zCap2)

pal2 <- designer.colors(256, brewer.pal(11,'PiYG'))

mainVec <- c("IRI Forecast RPSS",
			 "Probabilistic Known-ENSO EBF RPSS")

pdf(sprintf('%s/10_fieldRpssProb.pdf',plotdir),10,18)
par(mfrow=c(3,1))

image.plot(lon,lat,z[,,1],col=pal,ylim=c(-60,90), zlim=zr,main=mainVec[1],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)
abline(h=c(30,-30), lty=3)	
legend(-189, 92, '(a)', cex = textSize, bty='n')


image.plot(lon,lat,z[,,2],col=pal,ylim=c(-60,90), zlim=zr,main=mainVec[2],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)
abline(h=c(30,-30), lty=3)
legend(-189, 92, '(b)', cex = textSize, bty='n')


image.plot(lon, lat, z2, ylim=c(-60,90), zlim=zr2, col=pal2, 
	xlab='', ylab='', main='IRI Forecast RPSS (EBF reference)',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)
world(add=T)
abline(h=c(30,-30), lty=3)
legend(-189, 92, '(c)', cex = textSize, bty='n')

dev.off()

###############################################################################
# (11) RPSS alternative series comparison
###############################################################################

forceYR <- c(-0.22,0)
altCol <- 'red2'
aNino <- -0.22

pdf(sprintf('%s/11_rpssSeriesComp.pdf',plotdir),10,14)
par(mfrow=c(2,1), mar=c(5, 5, 4, 3) + 0.1)

x  <- tropics$timeMap[,3]
y  <- global$skill$rpss
y2 <- rpssComp$global$series[1:198]
y3 <- tropics$skill$rpss
y4 <- rpssComp$tropics$series[1:198]

totalVec <- c(as.numeric(totalTabGlobal[1,1]),rpssComp$global$total,
			  as.numeric(totalTabTropics[1,1]), rpssComp$tropics$total)

plot(x,y,ylim=range(y,y2,y3,y4,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Ranked Probability Skill Score',main='Mean Global Skill of IRI Forecast (RPSS)',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)

points(x,y2,type='l',col=altCol,lwd=lw)
abline(h=0)
legend('topleft', '(a)', cex = textSize, bty='n')

legend(2015.5,-0.01,as.character(signif(as.numeric(totalVec[1:2]),3)),text.col=c('black',altCol), bty='n')

legend(2001,0.22, 
	c('Climatology Reference',
	  'EBF Reference'),cex=textSize,
	col=c('black',altCol),lwd=lw,lty=c(1,1),bty='n')

addNino(y=c(-100,aNino,aNino,-100))


plot(x,y3,ylim=range(y,y2,y3,y4,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Ranked Probability Skill Score',main='Mean Tropics Skill of IRI Forecast (RPSS)',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)

points(x,y4,type='l',col=altCol,lwd=lw)
abline(h=0)
legend('topleft', '(b)', cex = textSize, bty='n')

legend(2015.5,-0.01,as.character(signif(as.numeric(totalVec[3:4]),3)),
	text.col=c('black',altCol), bty='n')

legend(2001,0.22, 
	c('Climatology Reference','EBF Reference'),cex=textSize,
	col=c('black',altCol),lwd=lw,lty=c(1,1),bty='n')

addNino(y=c(-100,aNino,aNino,-100))


dev.off()

###############################################################################
# (12) Lead time plot
###############################################################################

# source the results from the multi lead analysis
load(sprintf('%s/multiLead.Rda',ddir))

# plot up the results
pdf(sprintf('%s/12_leadTimeComp.pdf',plotdir),14,7)
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
					 'Forecast-ENSO EBF (Tropics)', 'Forecast-ENSO EBF (Global)'),
			col=c('blue','blue', 'red', 'red'),lwd=2,lty=c(1,3,1,3),cex=textSize,bty='n')

legend(0.67,0.043, '(b)', bty='n',cex=textSize)

dev.off()
