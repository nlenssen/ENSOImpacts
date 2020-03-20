source('Results/Namelists/namelist_Test.Rnl')

alphaLevels <- c(0.85,0.9,0.95,0.99)

alphaInd <- 1

analysisList <- system('ls -d Data/Impacts/Cru_0.5*',intern=TRUE)

# first, make an affected area series in time
percentageOrg <- array(NA, dim=c(12,2,2,length(analysisList)))

ninoMax <- ninaMax <- rep(NA, length(analysisList))

M <- c(11,8,11,11)
for(i in 1:length(analysisList)){

	load(sprintf('%s/affectedArea.Rda',analysisList[i]))
	load(sprintf('%s/ninaSeasonalInd.Rda',analysisList[i]))

	percentageOrg[,,,i] <- percentageArr[,,,alphaInd]

	ninoMax[i] <- sum(ensoArray[,,1]>1)/(12*M[i])
	ninaMax[i] <- sum(ensoArray[,,2]< -1)/(12*M[i])
}

# Do the same again but with the 2.5 data
analysisList25 <- system('ls -d Data/Impacts/Cru_2.5*',intern=TRUE)

# first, make an affected area series in time
percentageOrg25 <- array(NA, dim=c(12,2,2,length(analysisList25)))

ninoMax25 <- ninaMax25 <- rep(NA, length(analysisList))

for(i in 1:length(analysisList25)){

	load(sprintf('%s/affectedArea.Rda',analysisList25[i]))
	load(sprintf('%s/ninaSeasonalInd.Rda',analysisList25[i]))

	percentageOrg25[,,,i] <- percentageArr[,,,alphaInd]

	ninoMax25[i] <- sum(ensoArray[,,1]>1)/(12*M[i])
	ninaMax25[i] <- sum(ensoArray[,,2]< -1)/(12*M[i])
}


# Make line plots for four distinct seasons
plotdirTemp <- '/Users/lenssen/Dropbox/DEES/LisaWork/Meetings/November/AreaPlots'

analysisInds <- c(2,4)
analysisCol <- c('red', 'black')
seasonInds <- c(1,4,7,10)
seasonNames <- c('DJF', 'MAM', 'JJA', 'SON')
ninoNames <- c('EN','LN')
ninoLongNames <- c('El Nino', 'La Nina')


anomNames <- c('W','D')
anomLongNames <- c('Above Average', 'Below Average')

pdf(sprintf('%s/cruArea0.5.pdf',plotdirTemp),14,14)
set.panel(2,2)
for(el in 1:2){
	for(wd in 1:2){
		plot(1:4,NULL,xaxt = "n", xlab='', ylab='Global Area Affected',
		ylim=c(0,0.5),
		main=sprintf('%s, %s (0.5 degree)',ninoLongNames[el],anomLongNames[wd]))
		axis(1, at=1:4, labels=seasonNames)
		for(i in 1:2){
			points(1:4,percentageOrg[seasonInds,el,wd,analysisInds[i]],
				type='b',lwd=2,col=analysisCol[i])
		}
	}
}

legend('bottomright', c('1951-2016 Analysis (M=11)','1951-1996 Analysis (M=8)'),
	col=rev(analysisCol),lwd=2)
dev.off()

pdf(sprintf('%s/cruArea2.5.pdf',plotdirTemp),14,14)
set.panel(2,2)
for(el in 1:2){
	for(wd in 1:2){
		plot(1:4,NULL,xaxt = "n", xlab='', ylab='Global Area Affected',
		ylim=c(0,0.5),
		main=sprintf('%s, %s (2.5 degree)',ninoLongNames[el],anomLongNames[wd]))
		axis(1, at=1:4, labels=seasonNames)
		for(i in 1:2){
			points(1:4,percentageOrg25[seasonInds,el,wd,analysisInds[i]],
				type='b',lwd=2,col=analysisCol[i])
		}
	}
}

legend('bottomright', c('1951-2016 Analysis (M=11)','1951-1996 Analysis (M=8)'),
	col=rev(analysisCol),lwd=2)
dev.off()



el <- 1
wd <- 2
textSize <- 1.5
analysisInd <- 2
pdf(sprintf('%s/cruArea2.5EW.pdf',plotdirTemp),10,6)
set.panel(1,2)
plot(1:4,NULL,xaxt = "n", xlab='', ylab='Global Proportion Affected',
ylim=c(0,0.5),
main=sprintf('El Nino'),
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

axis(1, at=1:4, labels=seasonNames,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

points(1:4,percentageOrg25[seasonInds,1,1,4],
	type='l',lwd=2,col='deepskyblue3')
points(1:4,percentageOrg25[seasonInds,1,2,4],
	type='l',lwd=2,col='salmon4',lty=1)

plot(1:4,NULL,xaxt = "n", xlab='', ylab='',
ylim=c(0,0.5),
main=sprintf('La Nina'),
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

axis(1, at=1:4, labels=seasonNames,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

points(1:4,percentageOrg25[seasonInds,2,1,4],
	type='l',lwd=2,col='deepskyblue3',lty=1)
points(1:4,percentageOrg25[seasonInds,2,2,4],
	type='l',lwd=2,col='salmon4',lty=1)

dev.off()



points(1:4,percentageOrg25[seasonInds,1,1,4],
	type='b',lwd=2,col='deepskyblue3')
points(1:4,percentageOrg25[seasonInds,1,2,4],
	type='b',lwd=2,col='salmon4',lty=1)
points(1:4,percentageOrg25[seasonInds,2,1,4],
	type='b',lwd=2,col='deepskyblue3',lty=3)
points(1:4,percentageOrg25[seasonInds,2,2,4],
	type='b',lwd=2,col='salmon4',lty=3)


# make the legend
pdf(sprintf('%s/cruArea2.5EWLegend.pdf',plotdirTemp),10,6)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)

legend('bottom', c('Above Normal','Below Normal'),
col=c('deepskyblue3','salmon4'),lty=c(1,1),lwd=2.5,cex=1.5,horiz=T)

dev.off()