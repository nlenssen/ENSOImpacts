# source('Impacts/Namelists/namelist_Cru_0.5_1951_2016_alpha90.Rnl')
# plotdir <- '/Users/lenssen/Desktop/scratchPlots/testNewMaps'

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

oceanMask <- ifelse(is.na(ensoYears[,,1,1]),1,NA)

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

nName      <- c('Nino', 'Nina', 'Neutral')
ellaName   <- c('El Niño', 'La Niña', 'Neutral')

anomName   <- c('Above', 'Below')

textSize <- 1.25

# colors to plot with
aboveColor   <- c(brewer.pal(7,'GnBu'))[3:7]
aboveColor <- c(rgb(218,246,207,maxColorValue=255),
				rgb(191,244,168,maxColorValue=255),
				rgb(134,184,117,maxColorValue=255),
				rgb(87,148,201,maxColorValue=255),
				rgb(5,65,237,maxColorValue=255))
belowColor   <- c(brewer.pal(7,'YlOrBr'))[3:7]
belowColor <- c(rgb(251,249,76,maxColorValue=255),
				rgb(255,186,76,maxColorValue=255),
				rgb(197,133,65,maxColorValue=255),
				rgb(158,78,40,maxColorValue=255),
				rgb(112,55,14,maxColorValue=255))

hlColor <- rbind(aboveColor,belowColor)

nonSigColor  <- adjustcolor('black', alpha=0.15)
oceanColor   <- adjustcolor('paleturquoise1', alpha=0.15)
oceanColor   <- rgb(242,250,255,maxColorValue=255)

dryColor     <- adjustcolor('red1', alpha=0.25)
dryColor     <- adjustcolor(rgb(248,211,209,maxColorValue=255),alpha=0.5)


# hard coded for now, think about how I may want to change this
pValBreaks <- c(0.85,0.9,0.95,0.98,0.99,1)

probColors <- tim.colors(6)[1:5]
probBreaks <- c(1/3,0.5,2/3,0.8,0.9,1)

probBreaksCont <- seq(1/3,1, length=257)
probColorsCont <- designer.colors(256, brewer.pal(9,'YlGnBu'))

# lon x lat x season x el/la x high/low anom

for(s in 1:length(seasonInds)){
for(el in 1:2){

plotName <- sprintf('%02d%s%s',s,nName[el],seasonName[seasonInds[s]])

plotTitle <- sprintf('%s Historical Precipitation Anomalies (%s)', ellaName[el], seasonNameFull[seasonInds[s]])


# high field
hl <- 1
tempCounts <- sigCounts[,,seasonInds[s],el,hl]
tempCounts[tempCounts==-1] <- NA

tempProbHigh <- tempCounts/ensoYears[,,seasonInds[s],el]

# low field
hl <- 2
tempCounts <- sigCounts[,,seasonInds[s],el,hl]
tempCounts[tempCounts==-1] <- NA

tempProbLow <- tempCounts/ensoYears[,,seasonInds[s],el]

pdf(sprintf('%s/empiricalProbs/%sProb.pdf',plotdir,plotName),10,6)
image(lon,lat,dryMask[,,seasonInds[s]],col=dryColor,
	main=plotTitle,xlab='',ylab='',ylim=c(-60,90),
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

image(lon,lat,tempProbHigh,col=hlColor[1,],breaks=probBreaks, add=TRUE)
image(lon,lat,tempProbLow,col=hlColor[2,],breaks=probBreaks, add=TRUE)
#image(lon,lat,nonSig,col=nonSigColor,add=TRUE)

image(lon,lat,oceanMask,col=oceanColor,add=TRUE)
world(add=TRUE)
abline(h=c(30,-30), lty=3)

# optional labeling for paper plots
if(el==1){
	legend(-195, 95, '(a)', cex = textSize, bty='n')
} else{
	legend(-195, 95, '(b)', cex = textSize, bty='n')
}

dev.off()


}}


paperdir <- '/Users/lenssen/Dropbox/DEES/LisaWork/ensoPaper/Revision02/Supplement'
system(sprintf('cp %s/empiricalProbs/01NinoDJFProb.pdf %s/03a_NinoDJFProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/01NinaDJFProb.pdf %s/03b_NinaDJFProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/04NinoMAMProb.pdf %s/03a_NinoMAMProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/04NinaMAMProb.pdf %s/03b_NinaMAMProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/07NinoJJAProb.pdf %s/03a_NinoJJAProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/07NinaJJAProb.pdf %s/03b_NinaJJAProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/10NinoSONProb.pdf %s/03a_NinoSONProb.pdf',plotdir, paperdir))
system(sprintf('cp %s/empiricalProbs/10NinaSONProb.pdf %s/03b_NinaSONProb.pdf',plotdir, paperdir))