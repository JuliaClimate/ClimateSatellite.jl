module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related
# variables from various satellite instruments

## Modules Used
using Dates, Printf
using NetCDF, HDF5, FTPClient, PyCall
using ClimateEasy

## Exporting the following functions:
export
        clisatroot, pmmfinftpopen, pmmnrtftpopen, pmmftpclose,
        gpmfroot, gpmlroot, gpmeroot, trmmroot, mimicroot,
        gpmfrun, gpmlrun, gpmerun, trmmrun, mimicrun

## Including other files in the module
include("clisatinit.jl")

include("gpmfinal.jl")
include("gpmlate.jl")
include("gpmearly.jl")
include("trmm.jl")

include("mimic.jl")
#include("modis.jl")
#include("rss.jl")

end # module
