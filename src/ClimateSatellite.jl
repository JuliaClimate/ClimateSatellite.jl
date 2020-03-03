module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related
# variables from various satellite instruments

## Modules Used
using Logging, Dates
using DelimitedFiles, Printf, Statistics
using NCDatasets, HDF5, FTPClient
using GeoRegions

## Exporting the following functions:
export
        clisatroot, clisatlonlat, clisatinfo!, clisatdownload, clisatanalysis,
        clisatrawregion, clisatrawpoint, clisatrawgrid, clisatrawfol, clisatrawname,
        clisatanaregion, clisatanapoint, clisatanagrid, clisatanafol, clisatananame,
        clisatsubregion,
        gpmdwn, mimicdwn

## Including other files in the module
include("frontend.jl")
include("backend.jl")
include("raw.jl")
include("analysis.jl")
include("extract.jl")
include("gpm.jl")
include("mimic.jl")
#include("modis.jl")
#include("rss.jl")

end # module
