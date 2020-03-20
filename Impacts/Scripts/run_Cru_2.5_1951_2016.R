# Load in the run parameters from the namelist
source(sprintf('Impacts/Namelists/namelist_Cru_2.5_1951_2016.Rnl'))

# Detect enso years
source('Impacts/Code/identifyEnso.R')

# Run the bulk of the analysis
source('Impacts/Code/maskedImpactsAnalysis.R')

# Make plots if requested
if(plotting) source('Impacts/Code/plotImpactsAnalysis.R')