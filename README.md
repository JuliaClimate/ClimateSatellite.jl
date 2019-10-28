# ClimateSatellite

ClimateSatellite contains scripts and functions to download satellite
measurements and observations of climate data.  It can be configured to download
data for specific regions as defined in the ClimateTools dependency.

ClimateSatellite currently supports the retrieval of data from the following
satellites/missions:
* Global Precipitation Mission (GPM) Research Product
* Global Precipitation Mission (GPM) Late NRT Product
* Global Precipitation Mission (GPM) Early NRT Product
* Morphed Integrated Microwave Imagery at CIMSS (MIMIC)-TPW2

The retrieval of the following satellites/missions are in development:
* Tropical Rainfall Measuring Mission (TRMM)
* Remote Sensing Systems TM (RSS)
* MODerate resolution Imaging Spectroradiometer (MODIS)

The ClimateSatellite.jl package requires the following Julia dependencies:
* Dates, Printf
* NetCDF, HDF5, FTPClient
* ClimateEasy

Author(s):
* Nathanael Zhixin Wong: nathanaelwong@fas.harvard.edu
