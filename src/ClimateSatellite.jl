module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related
# variables from various satellite instruments

## Modules Used
using Dates, Printf, DelimitedFiles
using NCDatasets, HDF5, FTPClient, PyCall
using ClimateEasy

## Exporting the following functions:
export
        clisatinfo!, clisatdwn, clisatsave,
        clisatextractall, clisatextractpoint, clisatextractgrid,
        clisatroot, clisatfol, clisatncname, clisatlonlat,
        gpmdwn, mimicdwn

## Including other files in the module
#include("startup.jl")
include("general.jl")
include("gpm.jl")
#include("trmm.jl")
include("mimic.jl")
#include("modis.jl")
#include("rss.jl")

end # module
