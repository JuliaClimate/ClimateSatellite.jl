# ClimateSatellite

ClimateSatellite contains scripts and functions to download satellite
measurements and observations of climate data.  It can be configured to download
data for specific regions as defined in the ClimateTools dependency.

ClimateSatellite currently supports the retrieval of data from the following
satellites/missions:
* Global Precipitation Mission (GPM)
* Tropical Rainfall Measuring Mission (TRMM)
* Morphed Integrated Microwave Imagery at CIMSS (MIMIC)-TPW2

The retrieval of the following satellites/missions are forthcoming:
* Remote Sensing Systems TM (RSS)
* MODerate resolution Imaging Spectroradiometer (MODIS)

The ClimateSatellite.jl package requires the following Julia dependencies:
* Dates, Memento, Printf
* NetCDF, HDF5, FTPClient
* Conda, PyCall, PyPlot
* ClimateTools

Author(s):
* Nathanael Zhixin Wong: nathanaelwong@fas.harvard.edu
