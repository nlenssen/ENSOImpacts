source('Forecast/Namelists/namelist_FinalResults.Rnl')

# Build validation dataset
source('Forecast/Code/createValidation.R')

# Build the necessary forecast lists
source('Forecast/Code/ensoAlternativeForecast.R')
source('Forecast/Code/iriForecast.R')

# Verify the forecasts

if(verification){
	rm(list=ls())
	gc()
	source('Forecast/Namelists/namelist_FinalResults.Rnl')
	source('Forecast/Code/evaluateSkill.R')
}

# Option to verify skill at seasonal level

if(seasonalVerification){
	rm(list=ls())
	gc()
	source('Forecast/Namelists/namelist_FinalResults.Rnl')
	source('Forecast/Code/evaluateSkillSeasonal.R')
}