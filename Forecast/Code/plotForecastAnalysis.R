#load(sprintf('%s/skillScores.Rda',ddir))

# compare global rpss series
seriesSkillPlot(x=iriRpss$timeMap[,3], 
				s1=naiveEnsoRpss$global$series[1:198],
				s2=probEnsoRpss$global$series[1:198],
				s3=iriRpss$global$series,
				xTrim=realEnsoRpss$timeMap[,3],
				s4=realEnsoRpss$global$series,
				plotdir = plotdir,
				plotName='rpssSeriesGlobal',
				ylimInc=-0.5,
				ylab='RPSS',
				main=sprintf('Global RPSS Skill (%s Verification)',validationName))

# compare tropical rpss series
seriesSkillPlot(x=iriRpss$timeMap[,3], 
				s1=naiveEnsoRpss$tropics$series[1:198],
				s2=probEnsoRpss$tropics$series[1:198],
				s3=iriRpss$tropics$series,
				xTrim=realEnsoRpss$timeMap[,3],
				s4=realEnsoRpss$tropics$series,
				plotdir = plotdir,
				plotName='rpssSeriesTropics',
				ylimInc=0,
				ylab='RPSS',
				main=sprintf('Tropics RPSS Skill (%s Verification)',validationName))

# compare tropical groc series
seriesSkillPlot(x=iriRpss$timeMap[,3],
				s1=grocSeriesNaiveEnso[1:198],
				s2=grocSeriesProbEnso[1:198],
				s3=grocSeriesIRI,
				xTrim=realEnsoRpss$timeMap[,3],
				s4=grocSeriesRealEnso,
				plotdir = plotdir,
				plotName='grocSeries',
				ylimInc=0.35,
				ylab='GROC',
				main=sprintf('Tropics GROC Skill (%s Verification)',validationName))

# compare tropical EIR series
seriesSkillPlot(x=iriRpss$timeMap[,3], 
				s1=rep(NA,198),
				s2=probEnsoEir$tropics$series[1:198],
				s3=iriEir$tropics$series,
				xTrim=realEnsoRpss$timeMap[,3],
				s4=realEnsoEir$tropics$series,
				plotdir = plotdir,
				plotName='eirSeriesTropics',
				ylimInc=-0.2,
				ylab='Effective Interest Rate',
				main=sprintf('Tropics EIR Skill (%s Verification)',validationName))


# make rpss global fields
lon <- iriRpss$lon
lat <- iriRpss$lat

# Plot the annual average RPSS field
rpssFieldPlot(lon,lat,iriRpss$global$field,plotdir,'iriFieldRpss',
	main=sprintf('Annual Mean RPSS of IRI Forecast (%s Verification)',validationName))

rpssFieldPlot(lon,lat,naiveEnsoRpss$global$field,plotdir,'naiveEnsoFieldRpss',
	main=sprintf('Annual Mean RPSS of Deterministic ENSO Forecast (%s Verification)',validationName))

rpssFieldPlot(lon,lat,probEnsoRpss$global$field,plotdir,'probEnsoFieldRpss',
	main=sprintf('Annual Mean RPSS of Probabilistic ENSO Forecast (%s Verification)',validationName))

rpssFieldPlot(lon,lat,realEnsoRpss$global$field,plotdir,'realEnsoFieldRpss',
	main=sprintf('Annual Mean RPSS of Realistic ENSO Forecast (%s Verification)',validationName))

# Plot the best season RPSS field
rpssFieldPlot(lon,lat,iriRpss$bestMat,plotdir,'iriFieldRpssBest',
	main=sprintf('Best Season RPSS of IRI Forecast (%s Verification)',validationName))

rpssFieldPlot(lon,lat,naiveEnsoRpss$bestMat,plotdir,'naiveEnsoFieldRpssBest',
	main=sprintf('Best Season RPSS of Deterministic ENSO Forecast (%s Verification)',validationName))

rpssFieldPlot(lon,lat,probEnsoRpss$bestMat,plotdir,'probEnsoFieldRpssBest',
	main=sprintf('Best Season RPSS of Probabilistic ENSO Forecast (%s Verification)',validationName))

rpssFieldPlot(lon,lat,realEnsoRpss$bestMat,plotdir,'realEnsoFieldRpssBest',
	main=sprintf('Best Season RPSS of Realistic ENSO Forecast (%s Verification)',validationName))

# Plot the annual average EIR field
eirFieldPlot(lon,lat,iriEir$global$averageField,plotdir,'iriFieldEir',
	main=sprintf('Annual Mean Effective Interest Rate of IRI Forecast (%s Verification)',validationName))

eirFieldPlot(lon,lat,naiveEnsoEir$global$averageField,plotdir,'naiveEnsoFieldEir',
	main=sprintf('Annual Mean Effective Interest Rate of Deterministic ENSO Forecast (%s Verification)',validationName))

eirFieldPlot(lon,lat,probEnsoEir$global$averageField,plotdir,'probEnsoFieldEir',
	main=sprintf('Annual Mean Effective Interest Rate of Probabilistic ENSO Forecast (%s Verification)',validationName))

eirFieldPlot(lon,lat,realEnsoEir$global$averageField,plotdir,'realEnsoFieldEir',
	main=sprintf('Annual Mean Effective Interest Rate of Realistic ENSO Forecast (%s Verification)',validationName))

# Plot the best season EIR field
eirFieldPlot(lon,lat,iriEir$global$bestField,plotdir,'iriFieldEirBest',
	main=sprintf('Best Season Effective Interest Rate of IRI Forecast (%s Verification)',validationName))

eirFieldPlot(lon,lat,naiveEnsoEir$global$bestField,plotdir,'naiveEnsoFieldEirBest',
	main=sprintf('Best Season Effective Interest Rate of Deterministic ENSO Forecast (%s Verification)',validationName))

eirFieldPlot(lon,lat,probEnsoEir$global$bestField,plotdir,'probEnsoFieldEirBest',
	main=sprintf('Best Season Effective Interest Rate of Probabilistic ENSO Forecast (%s Verification)',validationName))

eirFieldPlot(lon,lat,realEnsoEir$global$bestField,plotdir,'realEnsoFieldEirBest',
	main=sprintf('Best Season Effective Interest Rate of Realistic ENSO Forecast (%s Verification)',validationName))



# Plot reliability Diagrams
reliabilityDiagram(iriReliability,plotdir,'iriReliability',
	main=sprintf('IRI Forecast Reliability (%s Verification)',validationName))

reliabilityDiagram(naiveEnsoReliability,plotdir,'naiveReliability',
	main=sprintf('Naive ENSO Forecast Reliability (%s Verification)',validationName))

reliabilityDiagram(probEnsoReliability,plotdir,'probReliability',
	main=sprintf('Probabilistic ENSO Forecast Reliability (%s Verification)',validationName))

reliabilityDiagram(realEnsoReliability,plotdir,'realReliability',
	main=sprintf('Realistic ENSO Forecast Reliability (%s Verification)',validationName))

# GROC field plots
grocFieldPlot(lon,lat,iriField,plotdir,'iriFieldGroc',
	main=sprintf('GROC of IRI Forecast (%s Verification)',validationName))

grocFieldPlot(lon,lat,ensoField,plotdir,'naiveEnsoFieldGroc',
	main=sprintf('GROC of Deterministic ENSO Forecast (%s Verification)',validationName))

grocFieldPlot(lon,lat,ensoProbField,plotdir,'probEnsoFieldGroc',
	main=sprintf('GROC of Probabilistic ENSO Forecast (%s Verification)',validationName))

grocFieldPlot(lon,lat,ensoRealField,plotdir,'realEnsoFieldGroc',
	main=sprintf('GROC of Realistic ENSO Forecast (%s Verification)',validationName))
