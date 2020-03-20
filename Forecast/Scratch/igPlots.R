source('Forecast/Namelists/namelist_CMAP.Rnl')
source('Forecast/Scratch/igFunctions.R')

load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))
load(sprintf('%s/dryMask.Rda',ddir))

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')
landPropMat <- landPropList$prop


plotdir <- '/Users/lenssen/Dropbox/DEES/LisaWork/ensoPaper/workingFigures/igPlots'
runAnalysis <- TRUE

## Final Figures to make (All RANKED Ignorance Scores)

if(runAnalysis){
# get the yz form of the forecasts
yzIri  <- yzForecast(iriForecastList,observedTercile)
yzDet  <- yzForecast(ensoForecastList,observedTercile)
yzProb <- yzForecast(ensoProbForecastList,observedTercile)

# hyperparameters for the analysis

binSize <- 0.05

# field calculations
fIri  <- rignField(yzIri)
fDet  <- rignField(yzDet)
fProb <- rignField(yzProb)

# series calculations
sIri  <- rignSeries(yzIri, landPropMat, tropics = FALSE)
sDet  <- rignSeries(yzDet, landPropMat, tropics = FALSE)
sProb <- rignSeries(yzProb, landPropMat, tropics = FALSE)

sIriTropics  <- rignSeries(yzIri, landPropMat, tropics = TRUE)
sDetTropics  <- rignSeries(yzDet, landPropMat, tropics = TRUE)
sProbTropics <- rignSeries(yzProb, landPropMat, tropics = TRUE)

# field decomp calculations
fdIri  <- rignFieldDecomp(yzIri, binSize)
fdDet  <- rignFieldDecomp(yzDet, binSize)
fdProb <- rignFieldDecomp(yzProb, binSize)

# series decomp calculations
sdIri  <- rignSeriesDecomp(yzIri, landPropMat, binSize, tropics=FALSE)
sdDet  <- rignSeriesDecomp(yzDet, landPropMat, binSize, tropics=FALSE)
sdProb <- rignSeriesDecomp(yzProb, landPropMat, binSize, tropics=FALSE) 

sdIriTropics  <- rignSeriesDecomp(yzIri, landPropMat, binSize, tropics=TRUE)
sdDetTropics  <- rignSeriesDecomp(yzDet, landPropMat, binSize, tropics=TRUE)
sdProbTropics <- rignSeriesDecomp(yzProb, landPropMat, binSize, tropics=TRUE)

# Total calculations
tIri  <- rignTotal(yzIri, landPropMat, tropics = FALSE)
tDet  <- rignTotal(yzDet, landPropMat, tropics = FALSE)
tProb <- rignTotal(yzProb, landPropMat, tropics = FALSE)

tIriTropics  <- rignTotal(yzIri, landPropMat, tropics = TRUE)
tDetTropics  <- rignTotal(yzDet, landPropMat, tropics = TRUE)
tProbTropics <- rignTotal(yzProb, landPropMat, tropics = TRUE)

tdIri  <- rignTotalDecomp(yzIri,landPropMat, binSize, tropics=FALSE)
tdDet  <- rignTotalDecomp(yzDet,landPropMat, binSize, tropics=FALSE)
tdProb <- rignTotalDecomp(yzProb,landPropMat, binSize, tropics=FALSE)

tdIriTropics  <- rignTotalDecomp(yzIri,landPropMat, binSize, tropics=TRUE)
tdDetTropics  <- rignTotalDecomp(yzDet,landPropMat, binSize, tropics=TRUE)
tdProbTropics <- rignTotalDecomp(yzProb,landPropMat, binSize, tropics=TRUE)

}
##########
# Figure 6
##########
# (a) tropics resolution (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob

# (b) tropics reliability (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob

forceYR <- c(0)
pal <- brewer.pal(4,'BrBG')
lw <- 1.75
textSize <- 1.5

pdf(sprintf('%s/relResSeries.pdf',plotdir),10,14)
par(mfrow=c(2,1))



x  <- iriForecastList$timeMap[,3]
y  <- sdIriTropics$res
y2 <- sdDetTropics$res[1:198]
y3 <- sdProbTropics$res[1:198]

totalVec <- c(tdIriTropics$res, tdProbTropics$res, tdDetTropics$res)

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Ignorance Resolution Score',main='Mean Tropics Resolution',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y2,type='l',col=pal[1],lwd=lw)
points(x,y3,type='l',col=pal[4],lwd=lw)
abline(h=0)

legend('topright',as.character(round(totalVec,3)),text.col=c('black',pal[c(4,1)]), bty='n')
legend('topleft', '(b)', cex = textSize, bty='n')




x  <- iriForecastList$timeMap[,3]
y  <- -sdIriTropics$rel
y2 <- -sdDetTropics$rel[1:198]
y3 <- -sdProbTropics$rel[1:198]

totalVec <- c(tdIriTropics$rel, tdProbTropics$rel, tdDetTropics$rel)

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Ignorance Reliability Score',main='Mean Tropics Reliability',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y2,type='l',col=pal[1],lwd=lw)
points(x,y3,type='l',col=pal[4],lwd=lw)
abline(h=0)

legend('bottomright',as.character(round(totalVec,3)),text.col=c('black',pal[c(4,1)]), bty='n')
legend('topleft', '(b)', cex = textSize, bty='n')

dev.off()


##########
# Figure 7
##########
# resolution fields for
#	- IRI
#	- ENSO Prop

lon <- iriForecastList$lon
lat <- iriForecastList$lat

maxCutoff <- 0.15

z <- array(NA, dim=c(length(lon),length(lat),3))

z[,,1]  <- fdIri$res
z[,,2]  <- fdDet$res
z[,,3]  <- fdProb$res

z[z>maxCutoff] <- maxCutoff
z[z<0.01] <- 0

zr <- c(0,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

mainVec <- c("IRI Resolution",
			 "ENSO Deterministic Forecast Resolution",
			 "ENSO Probabilistic Forecast Resolution")



pdf(sprintf('%s/fieldRes.pdf',plotdir),10,12)
par(mfrow=c(2,1))

image.plot(lon,lat,z[,,1],col=pal,ylim=c(-60,90),main=mainVec[1],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)	

image.plot(lon,lat,z[,,3],col=pal,ylim=c(-60,90),main=mainVec[3],xlab='',ylab='',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)

dev.off()


##########
# Figure 9
##########

# EIR Global Series (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob



forceYR <- c(0)
pal <- brewer.pal(4,'BrBG')
lw <- 1.75
textSize <- 1.5

pdf(sprintf('%s/eirSeries.pdf',plotdir),10,14)
par(mfrow=c(2,1))



x  <- iriForecastList$timeMap[,3]
y  <- sIri$eir
y2 <- sDet$eir[1:198]
y3 <- sProb$eir[1:198]

totalVec <- c(tIri$eir, tProb$eir, tDet$eir)

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Effective Interest Rate',main='Mean Global Skill',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y2,type='l',col=pal[1],lwd=lw)
points(x,y3,type='l',col=pal[4],lwd=lw)
abline(h=0)

legend('topright',as.character(round(totalVec,3)),text.col=c('black',pal[c(4,1)]), bty='n',
	text.width = strwidth("-1.000"), xjust = 1, yjust = 1)
legend('topleft', '(b)', cex = textSize, bty='n')


# EIR Tropics Series (With Totals)
# 	- IRI
# 	- ENSO Det
# 	- ENSO Prob

x  <- iriForecastList$timeMap[,3]
y  <- sIriTropics$eir
y2 <- sDetTropics$eir[1:198]
y3 <- sProbTropics$eir[1:198]

totalVec <- c(tIriTropics$eir, tProbTropics$eir, tDetTropics$eir)

plot(x,y,ylim=range(y,y2,y3,forceYR,na.rm=T),type='l',lwd=lw, lty=1,
	xlab='',ylab='Effective Interest Rate',main='Mean Tropics Skill',
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
grid(lwd=1.5)
points(x,y2,type='l',col=pal[1],lwd=lw)
points(x,y3,type='l',col=pal[4],lwd=lw)
abline(h=0)

legend('topright',as.character(round(totalVec,3)),text.col=c('black',pal[c(4,1)]), bty='n',
	text.width = strwidth("-1.000"), xjust = 1, yjust = 1)
legend('topleft', '(b)', cex = textSize, bty='n')

dev.off()

###########
# Figure 10
###########

# EIR global fields for
#	- IRI
#	- ENSO Prop

lon <- iriForecastList$lon
lat <- iriForecastList$lat

minCutoff <- 0
maxCutoff <- 0.25

z <- array(NA, dim=c(length(lon),length(lat),2))

z[,,1]  <- fIri$eir
z[,,2]  <- fProb$eir

z[z<minCutoff] <- minCutoff
z[z>maxCutoff] <- maxCutoff

zr <- c(minCutoff,maxCutoff)

pal <- c('white',designer.colors(256,brewer.pal(9,'YlGnBu')))

mainVec <- c("IRI Forecast Skill (EIR) ",
			 "ENSO Probabilistic Forecast Skill (EIR)")

pdf(sprintf('%s/fieldEir.pdf',plotdir),10,12)
par(mfrow=c(2,1))

image.plot(lon,lat,z[,,1],col=pal,ylim=c(-60,90),main=mainVec[1],xlab='',ylab='',zlim=zr,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)	
world(add=T)	

image.plot(lon,lat,z[,,2],col=pal,ylim=c(-60,90),main=mainVec[2],xlab='',ylab='',zlim=zr,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)
world(add=T)	

dev.off()


