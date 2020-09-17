# Seasonal Forecast Skill of ENSO Teleconnection Maps
#### Lenssen, N. J. L., L. Goddard, and S. Mason, Seasonal Forecast Skill of ENSO Teleconnection Maps. Wea. Forecasting, doi: https://doi.org/10.1175/WAF-D-19-0235.1.

## Quick Start Guide
The code is designed to be run from an `R` session with `ENSOImpacts/` as the working directory. This should elimate the need for adjusting paths when running the code on a new machine.


1) Install the required packages in `R` using the command
```
install.packages(c('fields','RColorBrewer','ncdf4','lubridate','rjson'))
```

2) Collect the necessary raw data products and run the data processing procedure
```
source('SharedCode/Scripts/cleanAllData.R')
```

3) Run the full analysis with the script
```
source('Results/Scripts/runEverything.R')
```

## Overview of Code Oranization

The codebase is organized into three major sub-projects:
- `SharedCode/` The data cleaning and initial processing 
- `Impacts/` The ENSO impacts analysis 
- `Forecast/` The generation and verification of the statistical ENSO-based forecasts 

The organization of  each of these three sub-projects is similar with each having subdirectories. 
- `Code/` containing the routines needed to run a specific task
- `Namelists/` adjustable parameter settings that are altered to run specific analyses with the routines in `Code/`
- `Scripts/` collections of `Namelists` and routines from `Code` to run full analyses

## Maproom of Results

Maps of ENSO precipitation impacts have been uploaded to the International Research Institute for Climate and Society (IRI) Data Library and are avaiabe to view at [this link](http://iridl.ldeo.columbia.edu/expert/home/.lenssen/.ensoTeleconnections/.prob/figviewer.html?my.help=more+options&map.T.plotvalue=Dec++-+Feb&map.ENSO.plotvalue=ElNino&map.anomaly.plotvalue=Above_Normal&map.significance.plotvalue=85&map.Y.units=degree_north&map.Y.plotlast=90N&map.here.x=0&map.here.y=0&map.url=X+Y+fig-+colors+coasts+-fig&map.domain=+%7B+%2FT+12.5+12.5+plotrange+%2FENSO+%2FElNino+plotvalue+%2Fanomaly+%2FAbove_Normal+plotvalue+%2Fsignificance+85+plotvalue+Y+-60+90+plotrange+%7D&map.domainparam=+%2Fplotaxislength+432+psdef+%2Fplotborder+72+psdef+%2FXOVY+null+psdef&map.zoom=Zoom&map.Y.plotfirst=60S&map.X.plotfirst=180W&map.X.units=degree_east&map.X.modulus=360&map.X.plotlast=180&map.prob.plotfirst=0&map.prob.units=unitless&map.prob.plotlast=1&map.newurl.grid0=X&map.newurl.grid1=Y&map.newurl.land=draw+coasts&map.newurl.plot=colors&map.plotaxislength=800&map.plotborder=72&map.fnt=NimbusSanLSymbol&map.fntsze=12&map.XOVY=auto&map.color_smoothing=1&map.framelbl=framelabelstart&map.framelabeltext=&map.iftime=25&map.mftime=25&map.fftime=200)
