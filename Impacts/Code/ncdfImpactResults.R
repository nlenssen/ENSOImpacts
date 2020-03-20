# source the namelist, load the data, and save the necessary field to the array
source('Impacts/Namelists/namelist_Cru_0.5_1951_2016.Rnl')
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

# do some additional setup in the first data to make the array the correct size

tempResults <- resultsList$propArray
netCdfArray <- array(NA, dim=c(dim(tempResults),length(alphaLevels)))

# put the data in the array for easy R -> ncdf work
netCdfArray[,,,,,1] <- tempResults

rm(tempResults)
gc()

source('Impacts/Namelists/namelist_Cru_0.5_1951_2016_alpha90.Rnl')
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

netCdfArray[,,,,,2] <- resultsList$propArray

source('Impacts/Namelists/namelist_Cru_0.5_1951_2016_alpha95.Rnl')
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

netCdfArray[,,,,,3] <- resultsList$propArray


source('Impacts/Namelists/namelist_Cru_0.5_1951_2016_alpha99.Rnl')
load(sprintf('%s/maskedAnalysisResults.Rda',ddir))

netCdfArray[,,,,,4] <- resultsList$propArray

# lon x lat x season x el/la x high/low anom


# write the actual netcdf file
library(ncdf4)

# define dimensions
londim   <- ncdim_def("lon","degrees_east",as.double(resultsList$lon)) 
latdim   <- ncdim_def("lat","degrees_north",as.double(resultsList$lat)) 
mondim   <- ncdim_def("time",'month',as.double(1:12))
eldim    <- ncdim_def("ENSO",'EN/LN',as.double(1:2))
andim    <- ncdim_def("anomaly",'AN/BN',as.double(1:2))
alphadim <- ncdim_def("significance",'alpha',as.double(alphaLevels))

# define variables
fillvalue <- 1e32
dlname <- "Probability of Anomaly"
prob_def <- ncvar_def("prob",'',list(londim,latdim,mondim,eldim,andim,alphadim),
						fillvalue,dlname,prec="single")

# create netCDF file and put arrays
ncfname <- sprintf('Data/Impacts/finalImpacts.nc')
ncout <- nc_create(ncfname,list(prob_def),force_v4=T)

# put variables
ncvar_put(ncout,prob_def,netCdfArray)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"lat","axis","Y")
ncatt_put(ncout,"time","axis","T")

# add global attributes
ncatt_put(ncout,0,"title",'Lenssen et al. 2020 Results')
ncatt_put(ncout,0,"institution",'IRI')
history <- paste("N. Lenssen", date(), sep=", ")
ncatt_put(ncout,0,"history",history)

# close the file, writing data to disk
nc_close(ncout)



# test

handle <- nc_open('Data/Impacts/finalImpacts.nc')

lon <- ncvar_get(handle, 'lon')
lat <- ncvar_get(handle, 'lat')

dat <- ncvar_get(handle, 'prob')

image.plot(lon,lat,dat[,,1,1,1,1])
world(add=T)
