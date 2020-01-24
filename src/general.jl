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
            info[ii,1:4] .= str[1:4]; info[ii,5] = vinfo[1,:];
            info[ii,6] = vinfo[2,:];  info[ii,7] = vinfo[3,:];
        else
            info[ii,:] .= "missing";
        end

    end

    ID = (info[:,1] .== productID);
    productattr["source"]   = info[ID,2][1]; productattr["short"] = info[ID,1][1];
    productattr["product"]  = info[ID,3][1]; productattr["varID"] = info[ID,5][1];
    productattr["variable"] = info[ID,6][1]; productattr["units"] = info[ID,7][1];
    productattr["standard"] = info[ID,5][1]; productattr["scale"] = info[ID,8][1];
    productattr["offset"]   = info[ID,9][1];

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
    nvar = length(info["variable"])

    if nlon != length(lon); error("$(Dates.now()) - nlon is $(nlon) but lon contains $(length(lon)) elements") end
    if nlat != length(lat); error("$(Dates.now()) - nlat is $(nlat) but lat contains $(length(lat)) elements") end

    var_var = AbstractVector{<:AbstractString}(undef,nvar)
    att_var = map(x->Dict(),1:nvar)

    for ii = 1 : nvar
        var_var[ii]                  = info["varID"][ii];
        att_var[ii]["units"]         = info["units"][ii]
        att_var[ii]["standard_name"] = info["standard"][ii]
        att_var[ii]["long_name"]     = info["variable"][ii]
        att_var[ii]["scale_factor"]  = info["scale"][ii]
        att_var[ii]["add_offset"]    = info["offset"][ii]
        att_var[ii]["missing_value"] = -32768
        att_var[ii]["_FillValue"]    = -32768
    end

    var_lon = "longitude"; att_lon = Dict("units"=>"degrees_east","long_name"="longitude");
    var_lat = "latitude";  att_lat = Dict("units"=>"degrees_north","long_name"="latitude");

    var_t = "time"; att_t = Dict("calendar"=>gregorian,"long_name"=>"time",
                                 "units"=>"minutes since $(Date.year(date))-$(Date.year(month))-1 0:0:0")

    if isfile(fnc)
        @info "$(Dates.now()) - Unfinished netCDF file $(fnc) detected.  Deleting."
        rm(fnc);
    end

    @debug "$(Dates.now()) - Creating $(info["source"]) $(info["product"]) $(info["variable"]) netCDF file $(fnc) ..."

    for ii = 1 : nvar
        nccreate(fnc,var_var,"longitude",nlon,"latitude",nlat,"time",nt,
                 atts=att_var[ii],t=NC_SHORT);
    end

    nccreate(fnc,var_lon,"longitude",nlon,atts=att_lon,t=NC_FLOAT);
    nccreate(fnc,var_lat,"latitude",nlat,atts=att_lat,t=NC_FLOAT);
    nccreate(fnc,var_t,"time",nlat,atts=att_t,t=NC_INT);

    @info "$(Dates.now()) - Saving $(info["source"]) $(info["product"]) $(info["variable"]) data to netCDF file $(fnc) ..."

    for ii = 1 : nvar
        ncwrite(data,fnc,var_var[ii]);
    end

    ncwrite(lon,fnc,var_lon);
    ncwrite(lat,fnc,var_lat);
    ncwrite(t,fnc,var_t);

    fol = clisatfol(info,date,region);
    @debug "$(Dates.now()) - Moving $(fnc) to data directory $(fol)"

    if isfile(joinpath(fol,fnc)); @info "$(Dates.now()) - An older version of $(fnc) exists in the $(fol) directory.  Overwriting." end

    mv(fnc,joinpath(fol,fnc),force=true);

end
