module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related
# variables from various satellite instruments

## Modules Used
using Dates, Printf, DelimitedFiles
using NCDatasets, HDF5, FTPClient, PyCall
using ClimateEasy

## Exporting the following functions:
export
<<<<<<< HEAD
        clisatroot, clisatrun, clisatread
=======
        clisatroot, pmmftpopen, pmmftpclose,
        gpmfroot, gpmlroot, gpmeroot, trmmroot, mimicroot,
        gpmffol, gpmlfol, gpmefol, trmmfol, mimicfol,
        gpmfncfile, gpmlncfile, gpmencfile, trmmncfile, mimicfile,
        gpmfrun, gpmlrun, gpmerun, trmmrun, mimicrun
>>>>>>> d01dccb762f090d5c50e1f3ee52aafec4d546f83

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
