# Make the ENSO seasonal series
source('SharedCode/Code/cleanDataENSO.R')

# Make the CRU seasonal series
source('SharedCode/Code/cleanDataCRU.R')

# Aggregate CRU to the 2.5 degree grid
source('SharedCode/Code/aggregateCRU.R')

# Clean the New et al Data
source('SharedCode/Code/cleanDataNew.R')

# Make the GPCC seasonal series at both resolutions using the script to handle
# loading multiple nc files for each analysis
source('SharedCode/Scripts/cleanGPCC.R')

# Make the GPCP seasonal series
source('SharedCode/Code/cleanGPCP.R')

# Clean the historical IRI seasonal forecast data
source('SharedCode/Code/cleanForecastIRI.R')

# Clean the historical IRI ENSO forecast data
source('SharedCode/Code/cleanForecastENSO.R')