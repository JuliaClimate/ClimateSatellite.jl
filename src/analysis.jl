"""
This file stores the analysis functions of ClimateSatellite.jl that are important for
climatological analysis of the downloaded data.

Current functionalities stored in this file include:
- Yearly, seasonal and monthly means, maximum, minimum, range and standard deviation:
    - On a general climatological basis
    - For a specified hour
    - Of the diurnal variability
- Saving the results of the abovementioned analysis into yearly NetCDF files
- Retrieval of these data from their respective NetCDF files

Functionalities that are in development include:
- Zonal, meridional and domain averages

Because there is progressively more data every year, analysis will be conducted on a yearly
basis instead of for the entire time-series.  This is to ensure generality and ease-of-use
as the datasets become larger and larger.

"""

function clisatanalysis(
    productID::AbstractString, year::Integer;
    varname::AbstractString,
    path::AbstractString="", region::="GLB"
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end
    if !isdir(clisatfol(info,region))
        error("$(Dates.now()) - No data has been downloaded from $(info["source"]) $(info["product"]) in the $(regionfullname(region))")
    end

    rawfol = clisatrawfol(productID,Date(year),region,path=dataroot);
    if !isdir(rawfol)
        error("$(Dates.now()) - There is no data for $(info["source"]) $(info["product"]) for $(yrmo2dir(dateii)).")
    end

    anafol = clisatanafol(productID,Date(year),region,path=dataroot);

    info = Dict{Any,Any}("root"=>dataroot); clisatinfo!(info,productID);

    @info "$(Dates.now()) - Extracting $(info["source"]) $(info["product"]) data for the entire $(regionfullname(region)) region ..."

    if sum(info["varID"] .== varname) == 0
        error("$(Dates.now()) - There is no varname identifier $(varname) for $(info["source"]) $(info["product"])")
    else; vID = info["varID"] .== varname;
        info["varID"] = info["varID"][vID][1];
        info["units"] = info["units"][vID][1];
        info["variable"] = info["variable"][vID][1];
        info["standard"] = info["standard"][vID][1];
        info["scale"] = info["scale"][vID][1];
        info["offset"] = info["offset"][vID][1];
    end

    cd(rawfol); fnc = glob("*.nc"); lfnc = length(fnc);

end
