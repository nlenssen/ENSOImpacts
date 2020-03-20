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

nName      <- c('Nino', 'Nina', 'Neutral')
ellaName   <- c('El Niño', 'La Niña', 'Neutral')

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

# lon x lat x season x el/la x high/low anom

for(s in 1:length(seasonInds)){
for(el in 1:3){
for(hl in 1:2){
tempField <- pValArray[,,seasonInds[s],el,hl]

# extrat the non-sig regions
nonSig <- ifelse(tempField==-1,1,NA)

# set the sig reigons to NA in the plotting field
tempField[tempField==-1] <- NA

plotName <- sprintf('%02d%s%s%s',s,nName[el],anomName[hl],seasonName[seasonInds[s]])

plotTitle <- sprintf('Probability of %s Normal Precipitation (%s, %s)',
	anomName[hl], ellaName[el], seasonNameFull[seasonInds[s]])

# plot the p-value masked by significance
pdf(sprintf('%s/pVals/%sPval.pdf',plotdir,plotName),10,6)
image.plot(lon,lat,tempField,col=hlColor[hl,],breaks=pValBreaks,
	main=plotTitle,xlab='',ylab='')
image(lon,lat,nonSig,col=nonSigColor,add=TRUE)
image(lon,lat,dryMask[,,seasonInds[s]],col=dryColor,add=TRUE)
world(add=TRUE)
dev.off()


tempCounts <- sigCounts[,,seasonInds[s],el,hl]
tempCounts[tempCounts==-1] <- NA

tempProb <- tempCounts/ensoYears[,,seasonInds[s],el]

pdf(sprintf('%s/empiricalProbs/%sProb.pdf',plotdir,plotName),10,6)
image.plot(lon,lat,tempProb,col=hlColor[hl,],breaks=probBreaks,
	main=plotTitle,xlab='',ylab='',ylim=c(-60,90),
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize, legend.cex=textSize)
image(lon,lat,nonSig,col=nonSigColor,add=TRUE)
image(lon,lat,dryMask[,,seasonInds[s]],col=dryColor,add=TRUE)
world(add=TRUE)
abline(h=c(30,-30), lty=3)

# optional labeling for paper plots
# if(hl==1){
# 	legend(-197, 95, '(a)', cex = textSize, bty='n')
# } else{
# 	legend(-197, 95, '(b)', cex = textSize, bty='n')
# }

dev.off()

# pdf(sprintf('%s/empiricalProbsCont/%sProb.pdf',plotdir,plotName),10,6)
# image.plot(lon,lat,tempProb,breaks=probBreaksCont,col=probColorsCont,
# 	main=plotTitle,xlab='',ylab='')
# image(lon,lat,nonSig,col=nonSigColor,add=TRUE)
# image(lon,lat,dryMask[,,seasonInds[s]],col=dryColor,add=TRUE)
# world(add=TRUE)
# dev.off()

}}}