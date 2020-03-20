#source('Impacts/Namelists/namelist_Cru_0.5_1951_2016_11_alpha90.Rnl')
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

# unpack the list (to facilitate multiple data formats)
lon <- resultsList$lon
lat <- resultsList$lat
timeMap <- resultsList$timeMap
sampleSize <- resultsList$sampleSize
pValArray <- resultsList$pValArray
dryMask   <- resultsList$dryMask
sigCounts <- resultsList$sigCounts
ensoYears <- resultsList$ensoYears

# create the names needed to loop through
seasonInds <- c(1,4,7,10)
seasonInds <- 1:12
seasonName <- c('DJF', 'JFM', 'FMA', 
				'MAM', 'AMJ', 'MJJ',
				'JJA', 'JAS', 'ASO',
				'SON', 'OND', 'NDJ')

seasonNameFull <- c('December-February','January-March','February-April',
					'March-May','April-June','May-July',
					'June-August','July-September','August-October',
					'September-November','October-December','November-January')

nName      <- c('Nino', 'Nina')
ellaName   <- c('El Niño', 'La Niña')

anomName   <- c('Above', 'Below')

textSize <- 1.5

# colors to plot with
aboveColor   <- c(brewer.pal(7,'GnBu'))[3:7]
belowColor   <- c(brewer.pal(7,'YlOrBr'))[3:7]
hlColor <- rbind(aboveColor,belowColor)

nonSigColor  <- adjustcolor('black',alpha=0.15)
oceanColor   <- adjustcolor('paleturquoise1',alpha=0.3)
dryColor     <- adjustcolor('darkred',alpha=0.15)

# hard coded for now, think about how I may want to change this
pValBreaks <- c(0.85,0.9,0.95,0.98,0.99,1)

probColors <- tim.colors(6)[1:5]
probBreaks <- c(1/3,0.5,2/3,0.8,0.9,1)

probBreaksCont <- seq(1/3,1,length=257)
probColorsCont <- designer.colors(256,brewer.pal(9,'YlGnBu'))



# rearrange a single plot
nlon <- length(lon)
lonBreak <- min(which(lon > -30))
newLon <- c(lon[lonBreak:nlon],lon[1:(lonBreak-1)]+360)

lonInds <- c(which(lon > -30), which(lon < -30))

for(s in 1:length(seasonInds)){
for(el in 1:2){
for(hl in 1:2){
tempField <- pValArray[,,seasonInds[s],el,hl]

# extrat the non-sig regions
nonSig <- ifelse(tempField==-1,1,NA)

# set the sig reigons to NA in the plotting field
tempField[tempField==-1] <- NA

plotName <- sprintf('%02d%s%s%s',s,nName[el],anomName[hl],seasonName[seasonInds[s]])

plotTitle <- sprintf('Probability of %s Normal Precipitation (%s, %s)',
	anomName[hl], ellaName[el], seasonNameFull[seasonInds[s]])


tempCounts <- sigCounts[,,seasonInds[s],el,hl]
tempCounts[tempCounts==-1] <- NA

tempProb <- tempCounts/ensoYears[,,seasonInds[s],el]

png(sprintf('%s/empiricalProbsShift/png%sProb.png',plotdir,plotName),1250,600)
image(newLon,lat,tempProb[lonInds,],col=hlColor[hl,],breaks=probBreaks,
	main=plotTitle,xlab='',ylab='',ylim=c(-60,90),
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)
image(newLon,lat,nonSig[lonInds,],col=nonSigColor,add=TRUE)
image(newLon,lat,dryMask[lonInds,,seasonInds[s]],col=dryColor,add=TRUE)
world(add=TRUE)
abline(h=c(30,-30), lty=3)
dev.off()
}}}
