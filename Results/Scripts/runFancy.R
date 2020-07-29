source('Results/Namelists/namelist_FancyDryMask.Rnl')
load(sprintf('Data/Forecast/%s_Data/skillScores.Rda',fcastProject))
source('Results/Namelists/namelist_FancyDryMask.Rnl')

# run the CRU data summary
source('Results/Code/submissionPlots.R')