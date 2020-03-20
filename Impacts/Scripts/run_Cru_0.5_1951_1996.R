# Load in the run parameters from the namelist
source(sprintf('Impacts/Namelists/namelist_Cru_0.5_1951_1996.Rnl'))

# Detect enso years
source('Impacts/Code/identifyEnso.R')

# Run the bulk of the analysis
source('Impacts/Code/maskedImpactsAnalysis.R')

# Make plots if requested
if(plotting) source('Impacts/Code/plotImpactsAnalysis.R')