"""
This file initializes the ClimateSatellite module by providing the basic
functions that are applicable to all the different satellites.  Current
functionalities include:
    - detection of home directory
    - ftp open for arthuhou

"""

# Root Functions
function clisatroot(product::AbstractString)

    path = joinpath("$(homedir())","research","CliSat",product);
    @info "$(Dates.now()) - No directory path was given.  Setting to default path: $(path) for ClimateSatellite data downloads."

    if isdir(path)
        @info "$(Dates.now()) - The default path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    else
        @info "$(Dates.now()) - The path $(path) does not exist.  Creating now ..."
        mkpath(path);
    end

    return path

end

function clisatroot(product::AbstractString,path::AbstractString)
    pdir = joinpath(path,product);
    if isdir(pdir)
        @info "$(Dates.now()) - The path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    else
        @warn "$(Dates.now()) - The path $(path) does not exist.  A new directory will be created here.  Therefore if you already have an existing repository for ClimateSatellite data, make sure that $(path) is the correct location."
        @info "$(Dates.now()) - Creating path $(path) ..."
        mkpath(pdir);
    end
    return pdir
end

# FTP Functions
function pmmftpopen(server::AbstractString)
    email = "natgeo.wong%40outlook.com"
    @info "$(Dates.now()) - Opening FTP request to $(server).pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@$(server).pps.eosdis.nasa.gov")
end

function pmmftpclose(ftp)
    @info "$(Dates.now()) - Closing FTP request."
    close(ftp)
end

# Run
function clisatrun(year::Integer;
                   dataroot::AbstractString="",
                   regions::AbstractArray=["GLB"],
                   product::AbstractString)

    if dataroot == ""; dataroot = clisatroot(product); end

    data,grid = clisatdwn(product,year,dataroot,regions);
    clisatsave(product,data,grid,year,dataroot,regions);

end

function clisatdwn(product::AbstractString)
