module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related
# variables from various satellite instruments

## Modules Used
using Dates
using DelimitedFiles
using FTPClient
using GeoRegions
using HDF5
using Logging
using NCDatasets
using Printf
using Statistics

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

end # module
