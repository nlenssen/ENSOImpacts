# Load in the run parameters from the namelist
source(sprintf('Impacts/Namelists/namelist_Cru_0.5_1951_2016.Rnl'))

# Detect enso years
source('Impacts/Code/identifyEnso.R')

# get the number of events
ninoCount <- unlist(lapply(ninoList,nrow))
ninaCount <- unlist(lapply(ninaList,nrow))

pdf('Figures/Impacts/counts11.pdf',10,5)
plot(1:12,ninoCount,lwd=1.5, type='b',col='red',ylim=c(7,16),
	xlab='Month', ylab='Number of Events',
	main='Number of ENSO events (1950-2016)')
points(1:12,ninaCount,lwd=1.5, type='b',col='blue')
grid()
abline(h=11,lwd=1.5)
legend(4.35,8,c('El Nino', 'La Nina'), col=c('red','blue'), lwd=1.5, horiz=T,bty='n')
dev.off()
