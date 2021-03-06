# Downloading Satellite Datasets

This page goes through the following:
1. `clisatdownload` functionality
2. Basic downloading
3. Email / Login requirements
4. Specifying directories

## 1. Using `clisatdownload`

Satellite datasets are downloaded using the function `clisatdownload`, which is the frontend wrapper for all dataset-specific functions.

!!! info "Note"
    `clisatdownload` is currently configured to download data for each **month**.

We first refer to the current list of supported datasets in `ClimateSatellite.jl`:

|    ID    | Mission | Product |
|  :---:   |  :---:  |  :---   |
| gpmimerg |   PMM   | Global Precipitation Mission - IMERGv6 Research |
| gpmlate  |   PMM   | Global Precipitation Mission - IMERGv6 NRT Late Run |
| gpmearly |   PMM   | Global Precipitation Mission - IMERGv6 NRT Early Run |
| mtpw2m   |  MIMIC  | Total Precipitable Water v2m |

## 2. Simplest Example: MIMIC-TPW2m

The simplest example is to download MIMIC-TPW2m data for one month (in this example, for January 2018).

```
julia> using ClimateSatellite
julia> clisatdownload("mptw2m",Date(2018,1))
[ Info: 2020-05-23T22:27:59.099 - No directory path was given.  Setting to default path: /Users/natgeo-wong/research/CliSat/mtpw2m for ClimateSatellite data downloads.
┌ Warning: 2020-05-23T22:27:59.099 - The path /Users/natgeo-wong/research/CliSat/mtpw2m does not exist.  A new directory will be created here.  Therefore if you already have an existing repository for ClimateSatellite data, make sure that /Users/natgeo-wong/research/CliSat/mtpw2m is the correct location.
└ @ ClimateSatellite ~/.julia/dev/ClimateSatellite/src/frontend.jl:21
[ Info: 2020-05-23T22:27:59.1 - Starting data download of TPW2m data for 2018-01-01 ...
...
```

## 3. Datasets with Email / Login Requirements

MIMIC-TPW2m is unique among all the listed satellite datasets in that data can be directly retrieved without needing email registration.  For other datasets, such as GPM IMERG, you need to register online using your email on their website for access to their datasets at the [Global Precipitation Measurement website](https://gpm.nasa.gov/data-access).

Using `ClimateSatellite.jl`, you can include your email as a keyword input argument.  If the dataset requires your email login, but you did not supply an email in the keyword argument, `ClimateSatellite.jl` will throw an error.

For example, the proper usage of `clisatdownload` for GPM IMERG data is:

```
julia> clisatdownload("gpmimerg",Date(2018,1),email="clisatjl@gmail.com")
[ Info: 2020-05-23T22:11:05.072 - No directory path was given.  Setting to default path: /Users/natgeo-wong/research/CliSat/gpmimerg for ClimateSatellite data downloads.
[ Info: 2020-05-23T22:11:05.073 - The default path /Users/natgeo-wong/research/CliSat/gpmimerg exists and therefore can be used as a directory for ClimateSatellite data downloads.
[ Info: 2020-05-23T22:11:05.074 - Starting data download of GPM IMERG data for 2018-01-01 ...
[ Info: 2020-05-23T22:11:05.074 - Opening FTP request to arthurhou.pps.eosdis.nasa.gov.
[ Info: 2020-05-23T22:11:05.074 - Entering IMERG directory for 2018-01-01 ...
...
```

But, if you did not specify the `email` keyword, the following error is thrown:

```
julia> clisatdownload("gpmimerg",Date(2018,1))
[ Info: 2020-05-23T22:09:25.558 - No directory path was given.  Setting to default path: /Users/natgeo-wong/research/CliSat/gpmimerg for ClimateSatellite data downloads.
[ Info: 2020-05-23T22:09:25.559 - The default path /Users/natgeo-wong/research/CliSat/gpmimerg exists and therefore can be used as a directory for ClimateSatellite data downloads.
ERROR: 2020-05-23T22:09:25.559 - Usage of this dataset requires an email address for login.  However, no email was provided.
```

!!! info "Note"
    The email `clisatjl@gmail.com` does not actually exist.  I just used it as a substitute for my own email.

## 4. Directory Configuration

By default, `ClimateSatellite.jl` will download and analyze data in the path given by `$(homedir())/research/CliSat`.  However, it is often the case that data for different projects need to be saved into different locations.  This is done using the `path` keyword, with the following as an example:

```
julia> ddir = "/n/kuangdss01/lab/clisat/"
julia> clisatdownload("gpmimerg",Date(2018,1),email="clisatjl@gmail.com",path=ddir)
[ Info: 2020-02-13T14:41:23.299 - The path /n/kuangdss01/lab/clisat/ exists and therefore can be used as a directory for ClimateSatellite data downloads.
[ Info: 2020-02-13T14:41:26.677 - Starting data download of GPM IMERG data for 2001-01-01 ...
...
[ Info: 2020-02-13T15:09:27.502 - PMM GPM IMERG data for the Tropical Belt region has been saved into file /n/kuangdss01/lab/clisat/gpmimerg/TRP/raw/2001/gpmimerg-TRP-200101.nc
```

!!! warning "Directory Creation:"
    By default it is assumed that the directory that the keyword `path` points to already exists prior to calling `clisatdownload`.  If the directory does not yet exist, it is necessary to create the directory.  You can do this with `clisatroot` (see the following section) or simply by creating `path` in your bash terminal.

### Managing Data Directories using `clisatroot`

By default, if `path` does not exist, then `clisatdownload` will throw an error:

```@repl
using ClimateSatellite, Dates
ddir = joinpath(homedir(),"test"); if isdir(ddir); rm(ddir,recursive=true) end
clisatdownload("mtpw2m",Date(2018,1),path=ddir)
```

This is to prevent the accidental downloading of huge amounts of data into directories / disks that cannot hold the data.  You can use `clisatroot` to force the creation of the directory:

```@repl
using ClimateSatellite
clisatroot("mtpw2m",joinpath(homedir(),"test"),create=true);
```
