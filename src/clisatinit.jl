"""
This file initializes the ClimateSatellite module by providing the basic
functions that are applicable to all the different satellites.  Current
functionalities include:
    - detection of home directory
    - ftp open for arthuhou

"""

# Root Functions
function clisatroot()

    svrstr = "/n/kuangdss01/users/nwong/data/";
    svrrun = "/n/kuangdss01/users/nwong/data/";
    dskdir = "/Volumes/CliNat-Sat";
    docdir = "/Users/natgeo-wong/Documents/research/data/";

    if     isdir(svrstr); return svrstr;
        @info "$(Dates.now()) - The path $(svrdir) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    elseif isdir(svrrun); return svrrun;
        @warn "$(Dates.now()) - The path $(svrdir) is not readable."
        @info "$(Dates.now()) - The path $(svrrun) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    elseif isdir(dskdir); return dskdir;
        @info "$(Dates.now()) - Not running on remote server.  Checking for external disks."
        @info "$(Dates.now()) - External disk $(dskdir) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    elseif isdir(docdir); return docdir;
        @info "$(Dates.now()) - Not running on remote server.  Checking for external disks."
        @info "$(Dates.now()) - External disks not found.  Using local research data directory $(docdir) for ClimateSatellite data downloads."
    else
        @error "$(Dates.now()) - The predefined directories in clisatroot.jl do not exist.  They are user-dependent, so please modify/customize accordingly."
    end

end

function clisatroot(path::AbstractString)
    if isdir(path)
        @info "$(Dates.now()) - The path $(path) exists and therefore can be used as a directory for ClimateSatellite data downloads."
    end
    return path
end

# FTP Functions
function ppmftpopen()
    email = "natgeo.wong%40outlook.com"
    @info "$(Dates.now()) - Opening FTP request to arthurhou.pps.eosdis.nasa.gov."
    return FTP("ftp://$(email):$(email)@arthurhou.pps.eosdis.nasa.gov")
end

function ppmftpclose(ftp)
    @info "$(Dates.now()) - Closing FTP request to arthurhou.pps.eosdis.nasa.gov."
    close(ftp)
end
