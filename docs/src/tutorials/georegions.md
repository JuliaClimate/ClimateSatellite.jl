# Using `GeoRegions.jl` in `ClimateSatellite.jl`

`ClimateSatellite.jl` uses `GeoRegions.jl` in order to specify regions for extraction during download and analysis.  This aids greatly in helping to save computer memory, especially for datasets with both high spatial and temporal resolution (e.g. GPM IMERGv6).

What follows here is a quick overview of the functionalities of `GeoRegions.jl` in the context of `ClimateSatellite.jl`.  This page goes through the following:
1. A brief introduction of `GeoRegions.jl`
2. How to specify GeoRegions in `ClimateSatellite.jl`

!!! tip "More information on `GeoRegions.jl`"
    More information on `GeoRegions.jl` can be found in the full documentation, which is available here:

## Brief Introduction fo `GeoRegions.jl`
