# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# clean the 0.5 degree data
gpccPptFile <- 'gpcc7_0.5_1950_2013_ppt.nc'
gpccStnFile <- 'gpcc7_0.5_1950_2013_stn.nc'
gpccOutFileName <- 'gpccSeasonal_0.5.Rda'
source('SharedCode/Code/cleanDataGPCC.R')

# clean the 2.5 degree data
gpccPptFile <- 'gpcc7_2.5_1950_2013_ppt.nc'
gpccStnFile <- 'gpcc7_2.5_1950_2013_stn.nc'
gpccOutFileName <- 'gpccSeasonal_2.5.Rda'
source('SharedCode/Code/cleanDataGPCC.R')