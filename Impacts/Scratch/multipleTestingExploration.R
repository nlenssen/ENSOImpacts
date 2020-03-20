source('Impacts/Namelists/namelist_Cru_2.5_1951_2016_11.Rnl')

# load in the results. The object with the p-values are 
# provided in resultsList$pValArray
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

# pull out some useful objects
lon <- resultsList$lon
lat <- resultsList$lat

# grab out one field to experiment with multiple testing
mn <- 1
el <- 1
ab <- 1

tempFieldRaw <- resultsList$pValArrayFull[,,mn,el,ab]

# figure out which zeros are due to low data
cutoffValue <- 10
lowDataMask <- ifelse(resultsList$sampleSize[,,mn] < cutoffValue,NA,1)

tempField <- tempFieldRaw * lowDataMask

# Now, tempField is the p-values of our multiple testing problem. Our first
# attempt at handling multiple testing is throwing off-the-shelf BH-procedure
# at is using the builtin R function

bhTest <- matrix(1-p.adjust(c(1-tempField), method='BH'),length(lon),length(lat))



plotdir <- '/Users/lenssen/Desktop/scratchPlots/fdrWork'

# fdr calculation call
fdr255 <- fdrFunction(tempField, 0.255, plotName = sprintf('%s/fViz26.pdf',plotdir))
fdr30 <- fdrFunction(tempField, 0.30, plotName = sprintf('%s/fViz30.pdf',plotdir))
fdr35 <- fdrFunction(tempField, 0.35, plotName = sprintf('%s/fViz35.pdf',plotdir))
fdr40 <- fdrFunction(tempField, 0.40, plotName = sprintf('%s/fViz40.pdf',plotdir) )
fdr50 <- fdrFunction(tempField, 0.50, plotName = sprintf('%s/fViz50.pdf',plotdir) )

# calculate the proportion of area from the 85% significance level with no correction
base85 <- ifelse(tempField > 0.85, 1, NA)

cosMat <- cos(matrix(lat,nrow=length(lon),ncol=length(lat),byrow=T)*(pi/180))

initialArea <- sum(base85 * cosMat, na.rm=T) / sum(cosMat)

area255 <- sum(fdr255 * cosMat, na.rm=T)/sum(cosMat) / initialArea
area30  <- sum(fdr30 * cosMat, na.rm=T)/sum(cosMat) / initialArea
area35  <- sum(fdr35 * cosMat, na.rm=T)/sum(cosMat) / initialArea
area40  <- sum(fdr40 * cosMat, na.rm=T)/sum(cosMat) / initialArea
area50  <- sum(fdr50 * cosMat, na.rm=T)/sum(cosMat) / initialArea
# the least stringent confidence is a significance of 0.748 (or a p-value
# of just over 0.25).
hist(bhTest)

# make some plots
zr <- c(0,0.5)

pValBreaks <- sort(1-c(0.65,0.75,0.85,0.9,0.95,0.985,1))

probColors <- rev(tim.colors(7)[1:6])

# make a series of plots showing what FDR is doing

# first plot the raw p-values
pdf(sprintf('%s/raw.pdf',plotdir), 10,6)
image.plot(lon,lat,1 - tempField, breaks = pValBreaks, col = probColors, ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño)')
world(add=T)
dev.off()

# then plot just the significant regions
pdf(sprintf('%s/sigMask.pdf',plotdir), 10,6)
image.plot(lon,lat,base85 , ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño), alpha=0.15')
world(add=T)
dev.off()

# then plot the regions with 
pdf(sprintf('%s/fdr26.pdf',plotdir), 10,6)
image.plot(lon,lat,fdr255, ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño), FDR=0.255')
world(add=T)
legend('bottomleft', paste('Sig. area prop:',round(area255,3)))
dev.off()

pdf(sprintf('%s/fdr30.pdf',plotdir), 10,6)
image.plot(lon,lat,fdr30, ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño), FDR=0.30')
world(add=T)
legend('bottomleft', paste('Sig. area prop:',round(area30,3)))
dev.off()

pdf(sprintf('%s/fdr35.pdf',plotdir), 10,6)
image.plot(lon,lat,fdr35, ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño), FDR=0.35')
world(add=T)
legend('bottomleft', paste('Sig. area prop:',round(area35,3)))
dev.off()

pdf(sprintf('%s/fdr40.pdf',plotdir), 10,6)
image.plot(lon,lat,fdr40, ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño), FDR=0.40')
world(add=T)
legend('bottomleft', paste('Sig. area prop:',round(area40,3)))
dev.off()

pdf(sprintf('%s/fdr50.pdf',plotdir), 10,6)
image.plot(lon,lat,fdr50, ylim=c(-60,90),
	xlab='',ylab='',main = 'DJF Above Normal Precipitation (El Niño), FDR=0.50')
world(add=T)
legend('bottomleft', paste('Sig. area prop:',round(area50,3)))
dev.off()
