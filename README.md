# **<div align="center">ClimateSatellite.jl</div>**

<p align="center">
  <a href="https://www.repostatus.org/#active">
    <img alt="Repo Status" src="https://www.repostatus.org/badges/latest/active.svg?style=flat-square" />
  </a>
  <a href="https://travis-ci.com/github/JuliaClimate/ClimateSatellite.jl">
    <img alt="Travis CI" src="https://travis-ci.com/JuliaClimate/ClimateSatellite.jl.svg?branch=master&style=flat-square">
  </a>
  <a href="https://github.com/JuliaClimate/ClimateSatellite.jl/actions?query=workflow%3ADocumentation">
    <img alt="Documentation Build" src="https://github.com/JuliaClimate/ClimateSatellite.jl/workflows/Documentation/badge.svg">
  </a>
  <br>
  <a href="https://mit-license.org">
    <img alt="MIT License" src="https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square">
  </a>
  <img alt="MIT License" src="https://img.shields.io/github/v/release/JuliaClimate/ClimateSatellite.jl">
  <a href="https://juliaclimate.github.io/ClimateSatellite.jl/stable/">
    <img alt="Latest Documentation" src="https://img.shields.io/badge/docs-stable-blue.svg?style=flat-square">
  </a>
  <a href="https://juliaclimate.github.io/ClimateSatellite.jl/dev/">
    <img alt="Latest Documentation" src="https://img.shields.io/badge/docs-latest-blue.svg?style=flat-square">
  </a>
</p>

**Created By:** Nathanael Wong (nathanaelwong@fas.harvard.edu)

**Note:** I am currently working on a new package, [NASAPrecipitation.jl](https://github.com/natgeo-wong/NASAPrecipitation.jl), that contains most of the current functionality (except MIMIC), and aims to expand to maybe the other datasets available from the NASA Precipitation Measurement Mission.  This repository is not currently being developed further, thought this may change in the future as more satellite datasets become relevant to my research.

## **Introduction**

`ClimateSatellite.jl` is a Julia package that aims to streamline the following processes:
* downloading of climate satellite mission data
* basic analysis (mean, maximum, minimum, standard deviation, etc.) of downloaded data
* extraction of data for a given **GeoRegion** (see `GeoRegions.jl` for more information)

`ClimateSatellite.jl` can be installed via
```
] add ClimateSatellite
```

## Supported Satellite Missions
Since my research work focuses mainly on precipitation and water vapour research, my efforts will mostly be on satellite missions involving these variables (i.e. PMM (GPM/TRMM) and MIMIC-TPW2m).  However, if you are interested in adding scripts for additional missions, feel free to submit a pull request.

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

A list of satellite data sources and their products (because some have data for multiple product types), can be found in the `satellites.txt` file, along with their properties, units of measurement, etc.  However, only the following combinations of Missions / Products are currently valid in ClimateSatellite - all other options will throw an error.

|    ID   | Mission | Product |
|  :---:  | :---: | --- |
| gpmimerg |  PMM  | Global Precipitation Mission - IMERGv6 Final Research Version |
| gpmlate  |  PMM  | Global Precipitation Mission - IMERGv6 Near Real-Time Late Run |
| gpmearly |  PMM  | Global Precipitation Mission - IMERGv6 Near Real-Time Early Run |
| mtpw2m   | MIMIC | Total Precipitable Water v2m |
