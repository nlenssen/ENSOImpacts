source('Forecast/Namelists/namelist_Test.Rnl')

# load in the data needed
load(sprintf('%s/iriForecast.Rda',ddir))
load(sprintf('%s/ensoForecast.Rda',ddir))
load(sprintf('%s/obsTercile.Rda',ddir))

load('Data/RawProcessed/ninaSeasonal.Rda')

# land area dataset
load('Data/RawProcessed/landProp_2.5.Rda')


plotdir <- '/Users/lenssen/Dropbox/DEES/LisaWork/Meetings/CompIssues'
# work only with RPSS for now

rpssClim <- rpss(iriForecastList,observedTercile)
rpssComp <- rpss(iriForecastList,observedTercile,ensoProbForecastList)
rpssComp2 <- rpss(ensoProbForecastList,observedTercile)


eirClim <- eir(iriForecastList,observedTercile)
eirComp <- eir(iriForecastList,observedTercile,ensoProbForecastList)
eirComp2 <- eir(ensoProbForecastList,observedTercile)


# get non-overlapping seasons
subSeasons <- c(2,5,8,11)
rpssClim$timeMap[,2]

subInds <- which(rpssClim$timeMap[,2] %in% subSeasons)

subTimeMat <- rpssClim$timeMap[subInds,3]

# 
ninaTimeInds <- which(ninaSeasonal$timeMap[,3] %in% subTimeMat)

subNina <- ninaSeasonal$nina34[ninaTimeInds]

ninoYears <- which(subNina>1)
ninaYears <- which(subNina < -1)


# make a scatter of the nina index and the difference

# just look tropics for now
difference <- rpssComp$tropics$series[subInds] - rpssClim$tropics$series[subInds]
difference2 <- rpssClim$tropics$series[subInds] - rpssComp2$tropics$series[subInds]

pdf(sprintf('%s/rpssDiffScatter.pdf',plotdir),10,7)
plot(abs(subNina),difference,xlim=c(0,3),ylim=range(difference),
	xlab='Absolute Nino34 Standardized Anomaly',
	ylab='IRI (ENSO Reference) - IRI (Climatology Reference)',
	main='Tropics RPSS Skill Added by IRI Forecast')
points(abs(subNina[ninaYears]),difference[ninaYears],col='blue')
points(abs(subNina[ninoYears]),difference[ninoYears],col='red')
abline(h=0)
abline(v=1,lty=3)

dev.off()

pdf(sprintf('%s/rpssDiffScatter2.pdf',plotdir),10,7)
plot(abs(subNina),difference2,xlim=c(0,3),ylim=range(difference2),
	xlab='Absolute Nino34 Standardized Anomaly',
	ylab='IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)',
	main='Tropics RPSS Skill Added by IRI Forecast')
points(abs(subNina[ninaYears]),difference2[ninaYears],col='blue')
points(abs(subNina[ninoYears]),difference2[ninoYears],col='red')
abline(h=0)
abline(v=1,lty=3)
dev.off()


# Same plots for EIR
difference <- eirComp$tropics$series[subInds] - eirClim$tropics$series[subInds]
difference2 <- eirClim$tropics$series[subInds] - eirComp2$tropics$series[subInds]

pdf(sprintf('%s/eirDiffScatter.pdf',plotdir),10,7)
plot(abs(subNina),difference,xlim=c(0,3),ylim=range(difference),
	xlab='Absolute Nino34 Standardized Anomaly',
	ylab='IRI (ENSO Reference) - IRI (Climatology Reference)',
	main='Tropics EIR Skill Added by IRI Forecast')
points(abs(subNina[ninaYears]),difference[ninaYears],col='blue')
points(abs(subNina[ninoYears]),difference[ninoYears],col='red')
abline(h=0)
abline(v=1,lty=3)

dev.off()

pdf(sprintf('%s/eirDiffScatter2.pdf',plotdir),10,7)
plot(abs(subNina),difference2,xlim=c(0,3),ylim=range(difference2),
	xlab='Absolute Nino34 Standardized Anomaly',
	ylab='IRI (Climatology Reference) - ENSO Impacts (Climatology Reference)',
	main='Tropics EIR Skill Added by IRI Forecast')
points(abs(subNina[ninaYears]),difference2[ninaYears],col='blue')
points(abs(subNina[ninoYears]),difference2[ninoYears],col='red')
abline(h=0)
abline(v=1,lty=3)
dev.off()