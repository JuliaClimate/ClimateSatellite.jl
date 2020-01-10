"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general use regardless of satellite and product type, which includes:
    * downloading data (calls different backend functions based on satellite/product type)
    * saving of data, with details of satellite/product attributes taken from a `Dict`
    * satellite/product attribute information retrieval

"""

function clisatinfo(productID::AbstractString,productattr::Dict)

end

function clisatdwn(
    year::Integer;
    productID::AbstractString, email::AbstractString,
    dataroot::AbstractString="", regions::Array{String,1}=["GLB"]
)

    if dataroot == ""; dataroot = clisatroot(product); end

    info = Dict();
    info = clisatinfo(productID,info);

    if info["satellite"] = "GPM"
        if     info["product"] = "early"; gpmedwn(regions,year,info,email);
        elseif info["product"] = "late";  gpmldwn(regions,year,info,email);
        elseif info["product"] = "final"; gpmfdwn(regions,year,info,email);
        end
    elseif info["satellite"] = "TRMM";    trmmdwn(regions,year,info,email);
    elseif info["satellite"] = "MIMIC";   mtpwdwn(regions,year,info);
    end

end

function clisatsave(
    data::Array{Real,3},
    region::AbstractString, info::Dict, yrmo::AbstractString
)

    fnc = clisatncfile(date,region,info);
    nlon = size(data,1); lon = info["longitude"];
    nlat = size(data,2); lat = info["latitude"];
    nt   = size(data,3);

    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_var = info["variable"]; att_var = Dict("units" => "mm/hr");
    var_lon = "longitude";      att_lon = Dict("units" => "degree");
    var_lat = "lattitude";      att_lat = Dict("units" => "degree");

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    @debug "$(Dates.now()) - Creating $(info["satellite"]) $(info["product"]) $(info["variable"]) netCDF file $(fnc) ..."
    nccreate(fnc,var_var,"nlon",nlon,"nlat",nlat,"ntime",nt,atts=att_prcp,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    @info "$(Dates.now()) - Saving $(info["satellite"]) $(info["product"]) $(info["variable"]) data to netCDF file $(fnc) ..."
    ncwrite(data,fnc,var_var);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);

    fol = gpmffol(date,sroot,reg);
    @debug "$(Dates.now()) - Moving $(fnc) to data directory $(fol)"

    if isfile("$(fol)/$(fnc)"); @info "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,"$(fol)/$(fnc)",force=true);

end
