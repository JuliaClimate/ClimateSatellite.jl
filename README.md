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

### Directories
By default, `ClimateSatellite.jl` saves all data into a `datadir` repository that is user-specified, or else it will otherwise default to
```
datadir="~/research/data/$(ID)"
```

### Regions
ClimateSatellite utlizes `GeoRegions.jl` to specify domains from which data is to be extracted.  If the option is not specified, then `ClimateSatellite` will assume that the
user wishes to process **global** data (which may not be wise especially for GPM due to the
large file sizes involved and memory required).

### Downloads
ClimateSatellite aims to streamline the downloading process by ensuring that the frontend
functions for the user are standard regardless of mission or product type.  For example, if
I wanted to download **global** GPM IMERG research data for January 2007, I would do it by
```
clisatdownload("gpmimerg",Date(2007,1),email=example@domain.com);
```

If I wanted to specify that the data is to go into a specific path `ddir`, I would do
```
clisatdownload("gpmimerg",Date(2007,1),email=example@domain.com,path=ddir);
```

And if I wanted to specify that I wanted to keep only information within the `TRP` domain
as specified in `GeoRegions.jl `(Tropical Belt, [N,S,W,E] = [30,-30,0,360]), I would do
```
clisatdownload("gpmimerg",Date(2007,1),email=example@domain.com,regions=["TRP"],path=ddir);
```

And if I wanted to do all the above, except for MIMIC-TPW2m data, I would do:
```
clisatdownload("mtpw2m",Date(2007,1),email=example@domain.com);
clisatdownload("mtpw2m",Date(2007,1),email=example@domain.com,path=ddir);
clisatdownload("mtpw2m",Date(2007,1),email=example@domain.com,regions=["TRP"],path=ddir);
```

### Extraction of Subdomains
However, sometimes you don't want to keep redownloading data (especially since due to file
size it can take a long time) and therefore sometimes you might want to extract data from
within a subdomain.

For example, in `GeoRegions.jl`, `SMT` (Sumatra, [N,S,W,E] = [6,-6,95,107]) is within the `TRP` domain.  If you have already downloaded data for the `TRP` domain, you can therefore extraction data for the `SMT` region using:
```
clisatsubregion("gpmimerg",Date(2007,1),region="SMT",path=ddir);
```
