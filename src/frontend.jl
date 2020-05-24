"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general
use regardless of satellite and product type, which includes:
- detection of home directory on machine
- opening an FTP request to arthuhou and jsimpson (GPM / TRMM data) servers
- creation of folders (both data folders and temporary directories)
- Extraction of all variable attribute information for a given satellite/product mission
- Extraction of variable attribute information for a given satellite/product mission

"""

# Root Functions
function clisatroot(productID::AbstractString)

    path = joinpath("$(homedir())","research","CliSat",productID);
    @info "$(Dates.now()) - No directory path was given.  Setting to default path: $(path) for ClimateSatellite data downloads."

    if isdir(path)
        @info "$(Dates.now()) - The default path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    else
        @warn "$(Dates.now()) - The path $(path) does not exist.  A new directory will be created here.  Therefore if you already have an existing repository for ClimateSatellite data, make sure that $(path) is the correct location."
        @info "$(Dates.now()) - Creating path $(path) ..."
        mkpath(path);
    end

    return path

end

function clisatroot(
    productID::AbstractString, path::AbstractString;
    create::Bool=false
)

    pdir = joinpath(path,productID);

    if isdir(path)

        @info "$(Dates.now()) - The path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
        if !isdir(pdir);
            @info "$(Dates.now()) - Creating path $(pdir) ..." mkpath(pdir);
        end

    else

        if create
            @warn "$(Dates.now()) - The path $(path) does not exist.  A new directory will be created here.  Therefore if you already have an existing repository for ClimateSatellite data, make sure that $(path) is the correct location."
            @info "$(Dates.now()) - Creating path $(pdir) ..."
            mkpath(pdir);
        else
            error("$(Dates.now()) - The path $(path) does not exist.  If you already have an existing repository for ClimateSatellite data, make sure that $(path) is the correct location.")
        end

    end

    return pdir

end

## Extraction of Information for Satellite Missions and Variables

function clisatinfo!(productattr::Dict,productID::AbstractString)

    fileinfo = readlines(joinpath(@__DIR__,"../data/info.csv"))
    info = Array{Any,2}(undef,length(fileinfo)-1,11)

    for ii = 1 : length(fileinfo)-1

        row = fileinfo[ii+1];
        if row[1] != '#'; str = split(row,",");
              info[ii,1:5] .= str[1:5]; vinfo = reshape(str[6:end],6,:);
              info[ii,6]  = vinfo[1,:][1]; info[ii,7]  = vinfo[2,:][1];
              info[ii,8]  = vinfo[3,:][1]; info[ii,9]  = vinfo[4,:][1];
              info[ii,10] = vinfo[5,:][1]; info[ii,11] = vinfo[6,:][1];
        else; info[ii,:] .= "N/A";
        end

    end

    ID = (info[:,1] .== productID);

    productattr["short"]   = info[ID,1][1];
    productattr["source"]  = info[ID,2][1];
    productattr["product"] = info[ID,3][1];
    productattr["gridres"] = parse.(Float64,info[ID,4][1]);
    productattr["dayfreq"] = parse.(Int8,info[ID,5][1]);

    productattr["varID"] = info[ID,6]; productattr["standard"] = info[ID,7];
    productattr["units"] = info[ID,9]; productattr["variable"] = info[ID,8];

    productattr["offset"] = parse.(Float64,info[ID,10]);
    productattr["scale"]  = parse.(Float64,info[ID,11]);

    return

end

function clisatvarinfo!(
    productattr::Dict, productID::AbstractString;
    varname::AbstractString
)

    if sum(productattr["varID"] .== varname) == 0
        error("$(Dates.now()) - There is no varname identifier $(varname) for $(info["source"]) $(info["product"])")
    else; vID = productattr["varID"] .== varname;
        productattr["varID"] = productattr["varID"][vID][1];
        productattr["units"] = productattr["units"][vID][1];
        productattr["variable"] = productattr["variable"][vID][1];
        productattr["standard"] = productattr["standard"][vID][1];
        productattr["scale"] = productattr["scale"][vID][1];
        productattr["offset"] = productattr["offset"][vID][1];
    end

    return

end

## Functions to define Containing Folder for a particular Satellite Data Source

function clisatregfol(
    productID::AbstractString,region::AbstractString;
    path::AbstractString=""
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    fol = joinpath(dataroot,region);

    if !isdir(fol)
        error("Data directory $(fol) does not exist.")
    end

    return fol

end

function clisatregfol(info::Dict,region::AbstractString)

    fol = joinpath(info["root"],region);

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(gregionfullname(region)) region does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatrawfol(
    productID::AbstractString,date::TimeType,region::AbstractString;
    path::AbstractString=""
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    fol = joinpath(dataroot,region,"raw",yr2str(date));

    if !isdir(fol)
        error("Data directory $(fol) does not exist.")
    end

    return fol

end

function clisatrawfol(info::Dict,date::TimeType,region::AbstractString)

    fol = joinpath(info["root"],region,"raw",yr2str(date));

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(gregionfullname(region)) region and year $(yr2str(date)) does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatanafol(
    productID::AbstractString,region::AbstractString;
    path::AbstractString=""
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    fol = joinpath(dataroot,region,"ana");

    if !isdir(fol)
        error("Data directory $(fol) does not exist.")
    end

    return fol

end

function clisatanafol(info::Dict,region::AbstractString)

    fol = joinpath(info["root"],region,"ana");

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(gregionfullname(region)) region does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end


# NetCDF Filenames
function clisatrawname(productID::AbstractString,date::TimeType,region::AbstractString)
    return "$(productID)-$(region)-$(yrmo2str(date)).nc"
end

function clisatrawname(info::Dict,date::TimeType,region::AbstractString)
    return "$(info["short"])-$(region)-$(yrmo2str(date)).nc"
end

function clisatananame(
    productID::AbstractString, varname::AbstractString,
    date::TimeType, region::AbstractString
)
    return "$(productID)-$(varname)-ana-$(region)-$(yr2str(date)).nc"
end

function clisatananame(
    info::Dict, varname::AbstractString,
    date::TimeType, region::AbstractString
)
    return "$(info["short"])-$(varname)-ana-$(region)-$(yr2str(date)).nc"
end

function clisatlonlat(info::Dict)

    if info["source"] == "PMM"
        if     isprod(info,"gpm");  lon,lat = gpmlonlat()
        elseif isprod(info,"3b42"); lon,lat = trmmlonlat()
        end
    elseif info["source"] == "MIMIC"; lon,lat = mimiclonlat()
    else
        error("$(Dates.now()) - The lonlat function for this particular satellite data has not been defined.")
    end

    return lon,lat

end
