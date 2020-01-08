module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related
# variables from various satellite instruments

## Modules Used
using Dates, Printf
using NCDatasets, HDF5, FTPClient, PyCall
using ClimateEasy

## Exporting the following functions:
export
        clisatroot, clisatrun, clisatread

## Including other files in the module
include("startup.jl")

include("gpmfinal.jl")
include("gpmlate.jl")
include("gpmearly.jl")
include("trmm.jl")

include("mimic.jl")
#include("modis.jl")
#include("rss.jl")

end # module
