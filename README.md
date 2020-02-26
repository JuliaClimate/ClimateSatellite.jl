# ClimateSatellite.jl

A Julia Package that downloads climate data from various satellite missions.  You can add it via
```
] add ClimateSatellite
```

### Contributors:
* Nathanael Zhixin Wong (Original Dev): nathanaelwong@fas.harvard.edu


## Overview
`ClimateSatellite.jl` contains scripts and functions to download satellite measurements and
observations of climate data.  It can be configured to download data for specific regions
as defined in the `ClimateEasy.jl` dependency.

Since my research work focuses mainly on precipitation and water vapour research, my efforts
will mostly be on satellite missions involving these variables (i.e. PMM (GPM/TRMM) and MIMIC-TPW2m).  However, if you are interested in adding scripts for additional missions, feel free to submit a pull request.

`ClimateSatellite.jl` currently supports the retrieval of data from the following satellites/
missions:
* Global Precipitation Mission (GPM) Research Product (GPM-FINAL)
* Global Precipitation Mission (GPM) Late NRT Product (GPM-LATE)
* Global Precipitation Mission (GPM) Early NRT Product (GPM-EARLY)
* Morphed Integrated Microwave Imagery at CIMSS (MIMIC)-TPW2m (MIMIC)

The retrieval of the following satellites/missions are in development:
* ~Tropical Rainfall Measuring Mission (TRMM)~
* Remote Sensing Systems TM (RSS)
* MODerate resolution Imaging Spectroradiometer (MODIS)

Note: Since GPM IMERG has essentially replaced TRMM 3B42, I am not planning on releasing a version that adds TRMM functionalities at this stage.


## Valid Satellites / Products
A list of satellite data sources and their products (because some have data for multiple
product types), can be found in the `satellites.txt` file, along with their properties,
units of measurement, etc.  However, only the following combinations of Missions / Products
are currently valid in ClimateSatellite - all other options will throw an error.

|    ID   | Mission | Product |
|  :---:  | :---: | --- |
| gpmimerg |  PMM  | Global Precipitation Mission - IMERGv6 Final Research Version |
| gpmlate  |  PMM  | Global Precipitation Mission - IMERGv6 Near Real-Time Late Run |
| gpmearly |  PMM  | Global Precipitation Mission - IMERGv6 Near Real-Time Early Run |
| mtpw2m   | MIMIC | Total Precipitable Water v2m |


## Workflow
By default, `ClimateSatellite.jl` saves all data into a `data` repository that is user-specified, or else it will otherwise default to
```
datadir="~/research/data/$(satellite_acronym)"
```
