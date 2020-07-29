# setup the workspace properly
source('Results/Namelists/namelist_FinalResults.Rnl')
load(sprintf('Data/Forecast/%s_Data/skillScores.Rda',fcastProject))
source('Results/Namelists/namelist_FinalResults.Rnl')

# run the CRU Data summary
source('Results/Code/cruDataSummary.R')

# run the plotting script
source('Results/Code/submissionPlots.R')