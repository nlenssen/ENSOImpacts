source('Results/Namelists/namelist_Test.Rnl')


plotdirTemp <- '/Users/lenssen/Desktop/test'




analysisList <- c('Data/Impacts/Cru_0.5_1951_2016_11_Data',
				  'Data/Impacts/New_0.5_1951_1996_8_Data')

# first, make an affected area series in time
percentageOrg <- array(NA, dim=c(12,2,2,length(analysisList)))

ninoMax <- ninaMax <- rep(NA, length(analysisList))

M <- c(11,8,11,11)
for(i in 1:length(analysisList)){

	load(sprintf('%s/affectedArea.Rda',analysisList[i]))
	load(sprintf('%s/ninaSeasonalInd.Rda',analysisList[i]))

	# get at the 85% confidence level
	percentageOrg[,,,i] <- percentageArr[,,,1]

	ninoMax[i] <- sum(ensoArray[,,1]>1)/(12*M[i])
	ninaMax[i] <- sum(ensoArray[,,2]< -1)/(12*M[i])
}


# Make line plots for four distinct seasons from the new analysis
el <- 1
wd <- 2
textSize <- 1.5
analysisInd <- 2
pdf(sprintf('%s/cruArea0.5EW.pdf',plotdirTemp),10,6)
set.panel(1,2)
plot(1:4,NULL,xaxt = "n", xlab='', ylab='Global Proportion Affected',
ylim=c(0,0.5),
main=sprintf('El Nino'),
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

axis(1, at=1:4, labels=seasonNames,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

points(1:4,percentageOrg[seasonInds,1,1,1],
	type='l',lwd=2,col='deepskyblue3')
points(1:4,percentageOrg[seasonInds,1,2,1],
	type='l',lwd=2,col='salmon4',lty=1)

plot(1:4,NULL,xaxt = "n", xlab='', ylab='',
ylim=c(0,0.5),
main=sprintf('La Nina'),
  cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

axis(1, at=1:4, labels=seasonNames,
	cex.lab=textSize, cex.axis=textSize, cex.main=textSize, cex.sub=textSize)

points(1:4,percentageOrg[seasonInds,2,1,1],
	type='l',lwd=2,col='deepskyblue3',lty=1)
points(1:4,percentageOrg[seasonInds,2,2,1],
	type='l',lwd=2,col='salmon4',lty=1)

legend('bottomright', c('Above Normal','Below Normal'),
col=c('deepskyblue3','salmon4'),lty=c(1,1),lwd=2.5,cex=1.5,horiz=F)

dev.off()

