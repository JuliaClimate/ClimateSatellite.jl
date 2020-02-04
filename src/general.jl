"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general use regardless of satellite and product type, which includes:
    * detection of home directory on machine
    * opening an FTP request to arthuhou and jsimpson (GPM / TRMM data) servers
    * creation of folders (both data folders and temporary directories)

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

function clisatroot(productID::AbstractString,path::AbstractString)

    pdir = joinpath(path,product);

    if isdir(path)

        @info "$(Dates.now()) - The path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
        if !isdir(pdir);
            @info "$(Dates.now()) - Creating path $(pdir) ..." mkpath(pdir);
        end

    else

        error("$(Dates.now()) - The path $(path) does not exist.  If you already have an existing repository for ClimateSatellite data, make sure that $(path) is the correct location.")

    end

    return pdir

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
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatregfol(info::Dict,region::AbstractString)

    fol = joinpath(info["root"],region);

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region does not exist."
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
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region and year $(yr2str(date)) does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatrawfol(info::Dict,date::TimeType,region::AbstractString)

    fol = joinpath(info["root"],region,"raw",yr2str(date));

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region and year $(yr2str(date)) does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatanafol(
    productID::AbstractString,date::TimeType,region::AbstractString;
    path::AbstractString=""
)

    if path == ""; dataroot = clisatroot(productID);
    else;          dataroot = clisatroot(productID,path);
    end

    fol = joinpath(dataroot,region,"ana",yr2str(date));

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region and year $(yr2str(date)) does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatanafol(info::Dict,date::TimeType,region::AbstractString)

    fol = joinpath(info["root"],region,"ana",yr2str(date));

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region and year $(yr2str(date)) does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisattmpfol(info::Dict)

    fol = joinpath(info["root"],"tmp");

    if !isdir(fol)
        @debug "$(Dates.now()) - Creating temporary directory $(fol)."; mkpath(fol);
    end

    return fol

end

# NetCDF Filenames
function clisatrawname(productID::AbstractString,date::TimeType,region::AbstractString);
    return "$(productID)-$(region)-$(yrmo2str(date)).nc"
end

function clisatrawname(info::Dict,date::TimeType,region::AbstractString);
    return "$(info["short"])-$(region)-$(yrmo2str(date)).nc"
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
