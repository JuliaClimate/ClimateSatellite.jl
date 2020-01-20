# ClimateSatellite.jl

A Julia Package that downloads climate data from various satellite missions.

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
* Tropical Rainfall Measuring Mission (TRMM)
* Remote Sensing Systems TM (RSS)
* MODerate resolution Imaging Spectroradiometer (MODIS)

## Installation
`ClimateSatellite.jl` v2 has not been added to the main JuliaRegistry yet.  You may fork the
repository, or install it using pkg via
```
] add https://github.com/natgeo-wong/ClimateSatellite.jl
```

## Valid Satellites / Products
A list of valid satellite data sources and their products (because some have data for
multiple product types), can be found in the `satellites.txt` file, along with their
properties, units of measurement, etc.

## Workflow
By default, `ClimateSatellite.jl` saves all data into a `data` repository that is user-specified, or else it will otherwise default to
```
datadir="~/research/data/$(satellite_acronym)"
```
