"""
This file initializes the ClimateSatellite module by providing the basic
functions that are applicable to all the different satellites.  Current
functionalities include:
    - detection of home directory
    - ftp open for arthuhou

"""

# Root Functions
function clisatroot()

    path = joinpath("$(homedir)","research","data";
    @info "$(Dates.now()) - No directory path was given.  Setting to default path: $(path) for ClimateSatellite data downloads."

    if isdir(path)
        @info "$(Dates.now()) - The default path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    else
        @warn "$(Dates.now()) - The path $(path) does not exist.  Creating now ..."
        mkpath(path);
    end

    return path

end

function clisatroot(path::AbstractString)
    if isdir(path)
        @info "$(Dates.now()) - The path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    end
    return path
end

# FTP Functions
function pmmfinftpopen()
    email = "natgeo.wong%40outlook.com"
    @info "$(Dates.now()) - Opening FTP request to arthurhou.pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@arthurhou.pps.eosdis.nasa.gov")
end

function pmmnrtftpopen()
    email = "natgeo.wong%40outlook.com"
    @info "$(Dates.now()) - Opening FTP request to jsimpson.pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@jsimpson.pps.eosdis.nasa.gov")
end

function pmmftpclose(ftp)
    @info "$(Dates.now()) - Closing FTP request to arthurhou.pps.eosdis.nasa.gov."
    close(ftp)
end
