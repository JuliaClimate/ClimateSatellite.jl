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
        clisatrawfol, clisatrawname, clisatrawread,
        clisatanafol, clisatananame, clisatanaread,
        clisatcmpfol, clisatcmpname, clisatcmpread,
        clisatrawregion, clisatrawpoint, clisatrawgrid,
        clisatanaregion, clisatanapoint, clisatanagrid,
        clisatsubregion,
        gpmdwn, gpmlonlat,
        mimicdwn, mimiclonlat

## Including other files in the module
include("frontend.jl")
include("backend.jl")
include("raw.jl")
include("analysis.jl")
include("extract.jl")
include("gpm.jl")
include("mimic.jl")

end # module
