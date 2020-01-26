"""
This file contains all the front-end scripts in ClimateSatellite.jl that are for general use regardless of satellite and product type, which includes:
    * detection of home directory on machine
    * opening an FTP request to arthuhou and jsimpson (GPM / TRMM data) servers
    * downloading data (calls different backend functions based on satellite/product type)
    * saving of data, with details of satellite/product attributes taken from a `Dict`
    * satellite/product attribute information retrieval
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

# FTP Functions
function pmmftpopen(server::AbstractString,email::AbstractString)
    @info "$(Dates.now()) - Opening FTP request to $(server).pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@$(server).pps.eosdis.nasa.gov")
end

function pmmftpclose(ftp)
    @info "$(Dates.now()) - Closing FTP request."
    close(ftp)
end

function isprod(info::Dict,product::AbstractString)
    occursin(product,info["short"])
end

function clisatfol(info::Dict,date::TimeType,region::AbstractString)

    fol = joinpath(info["root"],region,"raw",yr2str(date));

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region and year $(yr2str(date)) does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
    end

    return fol

end

function clisatfol(info::Dict,region::AbstractString)

    fol = joinpath(info["root"],region);

    if !isdir(fol)
        @info "$(Dates.now()) - $(info["source"]) $(info["product"]) data directory for the $(regionfullname(region)) region does not exist."
        @debug "$(Dates.now()) - Creating data directory $(fol)."; mkpath(fol);
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

function clisatrmtmp(flist::Array{<:AbstractString,2},fol::AbstractString)

    for ii = 1 : length(flist)
        fii = joinpath(fol,"$(flist[ii])"); if isfile(fii); rm(fii) end
    end

end

function clisatncname(info::Dict,date::TimeType,region::AbstractString);
    return "$(info["short"])-$(region)-$(yrmo2str(date)).nc"
end
