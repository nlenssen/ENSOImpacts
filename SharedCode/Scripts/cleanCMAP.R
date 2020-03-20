# Clear workspace and garbage collect
rm(list = ls())
gc()

# Load namelist
source('SharedCode/Namelists/cleanData.Rnl')

# clean the 0.5 degree data
cmapFile <- 'cmap07.2018v1_2.5_1979_2017.nc'
cmapOutFileName <- 'cmapSeasonalv1_2.5.Rda'
source('SharedCode/Code/cleanDataCMAP.R')

# clean the 2.5 degree data
cmapFile <- 'cmap07.2018v2_2.5_1979_2017.nc'
cmapOutFileName <- 'cmapSeasonalv2_2.5.Rda'
source('SharedCode/Code/cleanDataCMAP.R')