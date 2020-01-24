"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general use regardless of satellite and product type, which includes:
    * downloading data (calls different backend functions based on satellite/product type)
    * saving of data, with details of satellite/product attributes taken from a `Dict`
    * satellite/product attribute information retrieval
    * creation of folders (both data folders and temporary directories)

"""

function clisatinfo!(productattr::Dict,productID::AbstractString)

    fileinfo = readlines(joinpath(@__DIR__,"../data/info.txt"))
    info = Array{Any,2}(undef,length(fileinfo)-1,7)

    for ii = 1 : length(fileinfo)-1

        row = fileinfo[ii+1];
        if row[1] != '#'
            str = split(row,","); vinfo = reshape(str[5:end],3,:);
            info[ii,1:4] .= str[1:4];  info[ii,5] = [vinfo[1,:]];
            info[ii,6] = [vinfo[2,:]]; info[ii,7] = [vinfo[3,:]];
        else
            info[ii,:] .= "missing";
        end

    end

    ID = (info[:,1] .== productID);
    productattr["source"]   = info[ID,2][1]; productattr["short"] = info[ID,4][1];
    productattr["product"]  = info[ID,3][1]; productattr["varID"] = info[ID,6][1];
    productattr["variable"] = info[ID,5];    productattr["units"] = info[ID,7];

    return

end

function clisatfol(info::Dict,date::TimeType,region::AbstractString)

    fol = joinpath(info["root"],region,yr2str(date));

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region and year $(yr2str(date)) does not exist."
        @info "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisattmp(info::Dict)

    fol = joinpath(info["root"],"tmp");

    if !isdir(fol)
        @debug "$(Dates.now()) - Creating temporary directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatncname(info::Dict,date::TimeType,region::AbstractString);
    return "$(info["short"])-$(region)-$(ymd2str(date)).nc"
end

function clisatdwn(
    date::TimeType;
    productID::AbstractString, email::AbstractString,
    dataroot::AbstractString="", regions::Array{<:AbstractString,1}=["GLB"],
    overwrite::Bool=false
)

    if dataroot == ""; dataroot = clisatroot(productID); end

    info = Dict{Any,Any}("root"=>dataroot,"email"=>replace(email,"@"=>"%40"));
    clisatinfo!(info,productID);

    if info["source"] == "PMM"
        if     isprod(info,"gpm");  gpmdwn(regions,date,info,overwrite=overwrite);
        elseif isprod(info,"3b42"); trmmdwn(regions,date,info,overwrite=overwrite);
        end
    elseif info["source"] == "MIMIC"; mtpwdwn(regions,date,info);
    elseif info["source"] == "RSS"
        if     isprod(info,"trmm"); rtmidwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"gpm");  rgmidwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"smif"); rsmidwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"wind"); rwnddwn(regions,date,info,email,overwrite=overwrite);
        elseif isprod(info,"amsr"); rmsrdwn(regions,date,info,email,overwrite=overwrite);
        end
    end

    return info

end

function clisatsave(
    data::Array{Real,3},
    region::AbstractString, info::Dict, date::TimeType
)

    fnc = clisatncname(info,date,region);
    nlon = size(data,1); nlat = size(data,2); nt = size(data,3);

    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_var = info["variable"]; att_var = Dict("units" => "mm/hr");
    var_lon = "longitude";      att_lon = Dict("units" => "degree");
    var_lat = "lattitude";      att_lat = Dict("units" => "degree");

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    @debug "$(Dates.now()) - Creating $(info["source"]) $(info["product"]) $(info["variable"]) netCDF file $(fnc) ..."
    nccreate(fnc,var_var,"nlon",nlon,"nlat",nlat,"ntime",nt,atts=att_prcp,t=NC_FLOAT);
    nccreate(fnc,var_lon,"nlon",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"nlat",nlat,atts=att_lat,t=NC_FLOAT);

    @info "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) $(info["variable"]) data to netCDF file $(fnc) ..."
    ncwrite(data,fnc,var_var);
    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);

    fol = clisatfol(info,date,region);
    @debug "$(Dates.now()) - Moving $(fnc) to data directory $(fol)"

    if isfile(joinpath(fol,fnc)); @info "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,joinpath(fol,fnc),force=true);

end
