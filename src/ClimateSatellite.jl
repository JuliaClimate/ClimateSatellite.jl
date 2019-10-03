module ClimateSatellite

# Main file for the ClimateSatellite module that downloads climate-related variables from various satellite instruments

## Modules Used
using Dates, Memento, Printf
using NetCDF, HDF5, FTPClient
using Conda, PyCall, PyPlot
using ClimateTools

## Set Logging Format!
function __init__()
    global logger = Memento.config!("debug")
end

## Exporting the following functions:
export
        clisatroot,
        mimicroot, mimiclonlat, mimicfile, mimicfol,
        mimicdt, mimicget, mimicdwn, mimicextract, mimicsave, mimicrmtmp,
        gpmroot, gpmlonlat, gpmncfile, gpmhdf5, gpmfol,
        gpmftpopen, gpmftpcd, gpmftpclose,
        gpmdt, gpmget, gpmdwn, gpmextract, gpmsave, gpmrmtmp,
        mimicrun, gpmrun

## Including other files in the module
#include("clisatinit.jl")

include("gpm.jl")
#include("trmm.jl")

include("mimic.jl")
#include("modis.jl")
#include("rss.jl")

end # module
